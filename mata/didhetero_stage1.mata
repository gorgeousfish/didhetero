mata:

// =============================================================================
// Stage 1 dispatch: parametric estimation (GPS + OR) + KDE
//
// Functions:
//   1. didhetero_parametric_func()   - Unified GPS + OR estimation entry point
//   2. didhetero_stage1_dispatch()   - Full Stage 1 dispatch
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_parametric_func()
// Unified entry point for GPS + OR estimation.
// Calls didhetero_gps_estimate() and didhetero_or_estimate().
//
// Args:
//   data          - DidHeteroData struct
//   gteval        - K x 2 matrix of valid (g,t) pairs
//   geval         - K_g x 1 vector of valid treatment groups
//   control_group - "nevertreated" or "notyettreated"
//   anticipation  - anticipation periods (integer >= 0)
//
// Returns:
//   DidHeteroParamResults struct containing GPS and OR estimates
// -----------------------------------------------------------------------------
struct DidHeteroParamResults scalar didhetero_parametric_func(
    struct DidHeteroData scalar data,
    real matrix gteval,
    real colvector geval,
    string scalar control_group,
    real scalar anticipation)
{
    struct DidHeteroParamResults scalar result
    real matrix gps_coef, or_coef

    result = didhetero_init_param_results()
    result.ctrl_type = control_group

    // GPS estimation
    gps_coef = J(0, 0, .)
    result.gps_mat = didhetero_gps_estimate(data, gteval, geval,
        control_group, anticipation, gps_coef)
    result.gps_coef = gps_coef

    // OR estimation
    or_coef = J(0, 0, .)
    result.or_mat = didhetero_or_estimate(data, gteval,
        control_group, anticipation, or_coef)
    result.or_coef = or_coef

    return(result)
}

// -----------------------------------------------------------------------------
// didhetero_stage1_dispatch()
// Stage 1 dispatch: parametric estimation (GPS + OR) + KDE density/derivative.
//
// Orchestrates the full Stage 1 pipeline:
//   1. Parametric estimation via didhetero_parametric_func()
//   2. Z support grid generation via didhetero_gen_z_supp()
//   3. Kernel density estimation via didhetero_kde_density()
//   4. Density derivative estimation via didhetero_kde_deriv()
//   5. Assembly of DidHeteroStage1Results struct
//
// Args:
//   data          - DidHeteroData struct
//   gteval        - K x 2 matrix of valid (g,t) pairs
//   geval         - K_g x 1 vector of unique treatment groups
//   control_group - "nevertreated" or "notyettreated"
//   anticipation  - anticipation periods (integer >= 0)
//   zeval         - M x 1 vector of evaluation points for density estimation
//
// Returns:
//   DidHeteroStage1Results struct with GPS estimates, OR estimates, and KDE results
// -----------------------------------------------------------------------------
struct DidHeteroStage1Results scalar didhetero_stage1_dispatch(
    struct DidHeteroData scalar data,
    real matrix gteval,
    real colvector geval,
    string scalar control_group,
    real scalar anticipation,
    real colvector zeval)
{
    struct DidHeteroParamResults scalar param_results
    struct DidHeteroStage1Results scalar results
    real colvector Z, Z_supp, kd0_Z, kd1_Z
    real scalar r

    results = didhetero_init_stage1_results()

    // Step 1: Parametric estimation (GPS + OR)
    param_results = didhetero_parametric_func(data, gteval, geval, control_group, anticipation)

    // Step 2: Z vector (already extracted in data preparation)
    Z = data.Z

    // Step 3: Z_supp grid generation
    Z_supp = didhetero_gen_z_supp(Z)

    // Step 4: Kernel density estimation kd0_Z using Epanechnikov kernel
    kd0_Z = didhetero_kde_density(Z, zeval)
    // Positive value protection: truncate non-positive to 1e-12
    for (r = 1; r <= rows(kd0_Z); r++) {
        if (kd0_Z[r] <= 0) {
            printf("{txt}Warning: density estimate non-positive at z=%g, truncated to 1e-12\n", zeval[r])
            kd0_Z[r] = 1e-12
        }
    }

    // Step 5: Density derivative estimation kd1_Z using local polynomial (p=3, v=2)
    kd1_Z = didhetero_kde_deriv(Z, zeval)
    // Missing value check (warn only, no truncation)
    for (r = 1; r <= rows(kd1_Z); r++) {
        if (kd1_Z[r] == . | kd1_Z[r] >= .) {
            printf("{txt}Warning: density derivative estimate missing at z=%g, WLS matrix may be singular\n", zeval[r])
        }
    }

    // Step 6: Assemble return struct
    results.gps_mat  = param_results.gps_mat
    results.gps_coef = param_results.gps_coef
    results.or_mat   = param_results.or_mat
    results.or_coef  = param_results.or_coef
    results.ctrl_type = control_group
    results.Z_supp   = Z_supp
    results.kd0_Z    = kd0_Z
    results.kd1_Z    = kd1_Z

    return(results)
}

end
