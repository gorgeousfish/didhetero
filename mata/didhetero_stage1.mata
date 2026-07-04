mata:

// =============================================================================
// didhetero_stage1.mata
// Stage 1 dispatch: parametric estimation (GPS + OR) + KDE
//
// Provides:
//   - didhetero_parametric_func()  // Unified GPS + OR estimation entry point
//   - didhetero_stage1_dispatch()  // Full Stage 1 pipeline dispatch
//
// Requires:
//   - didhetero_types.mata         (DidHeteroData, DidHeteroParamResults,
//                                   DidHeteroStage1Results)
//   - didhetero_gps.mata           (didhetero_gps_estimate)
//   - didhetero_or.mata            (didhetero_or_estimate)
//   - didhetero_kde.mata           (didhetero_kde_density, didhetero_kde_deriv)
//   - didhetero_utils_init.mata    (didhetero_gen_z_supp,
//                                   didhetero_init_param_results)
//
// Paper reference: Section 2.3-2.4, parametric first stage
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
    // Data-driven density truncation with optional trim logic
    // Replaces hardcoded 1e-12 with max(density) * 1e-6 (paper Assumption 3.1)
    {
        real scalar _kde_trunc_count, _kde_max_density, _kde_threshold
        real colvector _kde_trimmed_vec

        _kde_max_density = max(kd0_Z)
        // Fallback: if all density values are non-positive, use 1e-12
        if (_kde_max_density <= 0 | _kde_max_density >= .) {
            _kde_threshold = 1e-12
        }
        else {
            _kde_threshold = _kde_max_density * 1e-6
        }

        _kde_trunc_count = 0
        _kde_trimmed_vec = J(rows(kd0_Z), 1, 0)

        for (r = 1; r <= rows(kd0_Z); r++) {
            if (kd0_Z[r] <= 0 | kd0_Z[r] < _kde_threshold) {
                _kde_trimmed_vec[r] = 1
                kd0_Z[r] = _kde_threshold
                _kde_trunc_count++
            }
        }
        if (_kde_trunc_count > 0) {
            printf("{txt}Note: %g/%g evaluation points below density threshold\n",
                _kde_trunc_count, rows(kd0_Z))
            if (data.verbose) {
                printf("{txt}      Threshold: %.2e, low-density z-values:", _kde_threshold)
                for (r = 1; r <= rows(kd0_Z); r++) {
                    if (_kde_trimmed_vec[r] == 1) printf(" %.4g", zeval[r])
                }
                printf("\n")
                printf("{txt}      These points may have unreliable variance estimates (1/f_Z inflation)\n")
            }
        }

        results.kde_trimmed = _kde_trimmed_vec
    }

    // Step 5: Density derivative estimation kd1_Z using local polynomial (p=3, v=2)
    kd1_Z = didhetero_kde_deriv(Z, zeval)
    // Missing value check (warn only, no truncation) — aggregated
    {
        real scalar _kd1_miss_count
        _kd1_miss_count = 0
        for (r = 1; r <= rows(kd1_Z); r++) {
            if (kd1_Z[r] == . | kd1_Z[r] >= .) _kd1_miss_count++
        }
        if (_kd1_miss_count > 0) {
            printf("{txt}Note: %g/%g density derivative estimates missing (WLS may be singular)\n",
                _kd1_miss_count, rows(kd1_Z))
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
