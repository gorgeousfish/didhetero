mata:

// =============================================================================
// Second and Third Stage Estimation: Influence Function Construction
//
// This module implements the doubly robust estimation pipeline for conditional
// average treatment effects on the treated (CATT). The estimation proceeds
// in two phases:
//
// Phase I (via didhetero_catt_core()):
//   - Construct intermediate variables
//   - Estimate nuisance parameters mu_E, mu_F via local linear regression
//   - Estimate conditional expectations mu_G, mu_R via local polynomial regression
//   - Construct A_{i,g,t}(z_r) and B_{i,g,t}(z_r) influence functions
//   - Compute DR point estimates via local polynomial regression on A
//
// Phase II (local computation):
//   - Standard error and analytical uniform confidence bands
//   - Bootstrap uniform confidence bands
//
// Functions:
//   - _didhetero_construct_A(): Construct DR core term A_{i,g,t}(z_r)
//   - _didhetero_construct_B(): Construct influence function B_{i,g,t}(z_r)
//   - didhetero_stage23(): Main estimation orchestrator
// =============================================================================


// -----------------------------------------------------------------------------
// _didhetero_construct_A()
// Construct A_{i,g,t}(z_r) for a single (g,t) pair and all evaluation points.
//
// The DR core term is defined as:
//   A_i(z_r) = (G_{i,g} / mu_G(z_r) - R_{i,g,t} / mu_R(z_r)) * Y_tilde_i
//
// where:
//   - Treatment group (G=1, R=0): A_i = Y_tilde_i / mu_G(z_r)
//   - Comparison group (G=0, R>0): A_i = -R_i * Y_tilde_i / mu_R(z_r)
//   - Other units (G=0, R=0): A_i = 0
//
// Arguments:
//   G_ig    - n x 1, treatment group indicator
//   R_g     - n x 1, IPW weight
//   Y_diff  - n x 1, outcome adjustment Y_tilde
//   mu_G    - R x 1, conditional mean of G at evaluation points
//   mu_R    - R x 1, conditional mean of R at evaluation points
//
// Returns:
//   n x R matrix of A_{i,g,t}(z_r) values
// -----------------------------------------------------------------------------
real matrix _didhetero_construct_A(
    real colvector G_ig,
    real colvector R_g,
    real colvector Y_diff,
    real colvector mu_G,
    real colvector mu_R)
{
    real scalar n, num_zeval, r
    real matrix A_mat

    n = rows(G_ig)
    num_zeval = rows(mu_G)
    A_mat = J(n, num_zeval, 0)

    for (r = 1; r <= num_zeval; r++) {
        A_mat[., r] = (G_ig / mu_G[r] - R_g / mu_R[r]) :* Y_diff
    }

    return(A_mat)
}


// -----------------------------------------------------------------------------
// _didhetero_construct_B()
// Construct B_{i,g,t}(z_r) for a single (g,t) pair and all evaluation points.
//
// The influence function with Hajek bias correction is defined as:
//   B_i(z_r) = A_i(z_r) + (mu_E(z_r) / mu_R(z_r)^2) * R_{i,g,t}
//                      - (mu_F(z_r) / mu_G(z_r)^2) * G_{i,g}
//
// Arguments:
//   A_mat   - n x R matrix of A_{i,g,t}(z_r)
//   G_ig    - n x 1, treatment group indicator
//   R_g     - n x 1, IPW weight
//   mu_G    - R x 1, conditional mean of G
//   mu_R    - R x 1, conditional mean of R
//   mu_E    - R x 1, conditional mean of E
//   mu_F    - R x 1, conditional mean of F
//
// Returns:
//   n x R matrix of B_{i,g,t}(z_r) values
// -----------------------------------------------------------------------------
real matrix _didhetero_construct_B(
    real matrix A_mat,
    real colvector G_ig,
    real colvector R_g,
    real colvector mu_G,
    real colvector mu_R,
    real colvector mu_E,
    real colvector mu_F)
{
    real scalar n, num_zeval, r
    real matrix B_mat

    n = rows(G_ig)
    num_zeval = rows(mu_G)
    B_mat = J(n, num_zeval, 0)

    for (r = 1; r <= num_zeval; r++) {
        B_mat[., r] = A_mat[., r] ///
            + (mu_E[r] / mu_R[r]^2) * R_g ///
            - (mu_F[r] / mu_G[r]^2) * G_ig
    }

    return(B_mat)
}


// -----------------------------------------------------------------------------
// didhetero_stage23()
// Main estimation orchestrator for the second and third stages.
//
// This function coordinates the doubly robust estimation pipeline by delegating
// core computations to didhetero_catt_core() and then computing inference
// quantities locally. The workflow comprises:
//
// Core estimation (via didhetero_catt_core()):
//   - Construct intermediate variables
//   - Estimate nuisance parameters mu_E, mu_F via local linear regression
//   - Estimate conditional expectations mu_G, mu_R via local polynomial regression
//   - Construct A_{i,g,t}(z_r) and B_{i,g,t}(z_r) influence functions
//   - Compute DR point estimates via local polynomial regression on A
//
// Inference computation (local):
//   - Standard errors and analytical uniform confidence bands (per group-time pair)
//   - Bootstrap uniform confidence bands (across all group-time pairs)
//
// Arguments:
//   data          - DidHeteroData struct (modified in place)
//   gps_mat       - Generalized propensity score matrix from Stage 1
//   or_mat        - Outcome regression matrix from Stage 1
//   bw            - K x 1 common bandwidth vector (pre-resolved)
//   bwselect      - bandwidth selection method string (for interface compatibility)
//   seed          - bootstrap RNG seed (applied when bootstrap is requested)
//
// Modified fields (via data struct):
//   data.A_g_t    - 1 x K pointer array, each n x num_zeval
//   data.B_g_t    - 1 x K pointer array, each n x num_zeval
//   data.G_g      - n x K group indicator matrix
//   data.mu_G_g   - num_zeval x K conditional mean of G
//
// Returns:
//   num_zeval x K matrix of DR point estimates
// -----------------------------------------------------------------------------
real matrix didhetero_stage23(
    struct DidHeteroData scalar data,
    real matrix gps_mat,
    real matrix or_mat,
    real colvector bw,
    string scalar bwselect,
    real scalar seed)
{
    real scalar n, K, num_zeval, id_gt
    real scalar h_gt, se_c_hat_gt
    real matrix est
    real colvector se_gt, ci1_l_gt, ci1_u_gt
    real colvector mathcal_V_gt
    real scalar has_kd0_Z
    real scalar r_pos, n_pos_missing, eps_pos
    struct DidHeteroCattResult scalar catt_result

    // --- Dimensions ---
    n = data.n
    K = data.num_gteval
    num_zeval = data.num_zeval

    // =====================================================================
    // Core CATT estimation via didhetero_catt_core()
    // Delegates intermediate variable construction, nuisance parameter
    // estimation, influence function construction, and DR point estimation
    // to the reentrant core function. Bandwidth is pre-resolved, hence
    // bwselect="manual" is specified.
    // =====================================================================
    catt_result = didhetero_catt_core(data, gps_mat, or_mat,
                      data.gteval, bw, "manual", data.porder, data.kernel)

    // --- Copy core results into data struct ---
    data.A_g_t = catt_result.A_g_t
    data.B_g_t = catt_result.B_g_t
    data.G_g   = catt_result.G_g
    data.mu_G_g = catt_result.mu_G_g
    data.mu_E_g_t = catt_result.mu_E_g_t
    data.mu_F_g_t = catt_result.mu_F_g_t
    est = catt_result.catt_est

    // Initialize SE and analytical UCB output matrices
    data.se = J(num_zeval, K, .)
    data.mathcal_V = J(num_zeval, K, .)
    data.ci1_lower = J(num_zeval, K, .)
    data.ci1_upper = J(num_zeval, K, .)
    data.c_hat = J(K, 1, .)

    // Check if kd0_Z is available for SE computation
    has_kd0_Z = (rows(data.kd0_Z) == num_zeval)

    // =====================================================================
    // Standard errors and analytical uniform confidence bands
    // Loop over group-time pairs using B_g_t and estimates from core results.
    // Computed only when kd0_Z is available from Stage 1.
    // =====================================================================
    for (id_gt = 1; id_gt <= K; id_gt++) {
        if (has_kd0_Z) {
            h_gt = bw[id_gt]

            didhetero_se_analytical_ucb(
                *data.B_g_t[id_gt], data.Z, data.zeval, data.Z_supp,
                data.kd0_Z,
                h_gt, n, data.porder, data.kernel,
                data.const_V, data.lambda, data.alp,
                est[., id_gt],
                se_gt, ci1_l_gt, ci1_u_gt, se_c_hat_gt,
                mathcal_V_gt)

            data.se[., id_gt] = se_gt
            data.mathcal_V[., id_gt] = mathcal_V_gt
            data.ci1_lower[., id_gt] = ci1_l_gt
            data.ci1_upper[., id_gt] = ci1_u_gt
            data.c_hat[id_gt] = se_c_hat_gt

            // =================================================================
            // Positivity guard for standard errors and confidence intervals
            //
            // When the generalized propensity score mu_G = E[G_ig | Z=z_r] falls
            // below a threshold, the overlap assumption is violated and the
            // CATT is not point-identified at that evaluation point. The
            // influence function B_i is inversely proportional to mu_G, so
            // small values cause variance inflation.
            //
            // The threshold eps_pos = 1e-3 identifies GPS values below 0.1%,
            // indicating evaluation points outside the treated cohort's support.
            // Standard errors and confidence intervals are set to missing to
            // indicate lack of point identification.
            // =================================================================
            eps_pos = 1e-3
            n_pos_missing = 0
            for (r_pos = 1; r_pos <= num_zeval; r_pos++) {
                if (data.mu_G_g[r_pos, id_gt] < eps_pos) {
                    data.se[r_pos, id_gt] = .
                    data.ci1_lower[r_pos, id_gt] = .
                    data.ci1_upper[r_pos, id_gt] = .
                    n_pos_missing++
                }
            }
            if (n_pos_missing > 0) {
                printf("Warning: SE=. at %g of %g evaluation point(s): " +
                       "GPS positivity violation (E[G=g|Z=z] < 1e-3)\n",
                       n_pos_missing, num_zeval)
                printf("  (zeval outside treated cohort support; " +
                       "SE reported as missing)\n")
            }
        }
    }

    // =================================================================
    // Bootstrap uniform confidence bands
    // When bootstrap iterations (biters) are specified, compute uniform
    // confidence bands via weighted bootstrap. Otherwise, set CI2 fields
    // to missing values.
    // =================================================================
    if (data.biters > 0) {
        if (seed >= 0 & seed < .) {
            rseed(seed)
        }

        didhetero_boot_ucb_optimized(
            data.A_g_t, est, data.se, bw, data.Z, data.zeval,
            n, data.porder, data.kernel, data.alp, data.biters,
            data.uniformall, K, num_zeval,
            data.ci2_lower, data.ci2_upper, data.c_check_bs)
    }
    else {
        data.ci2_lower = J(num_zeval, K, .)
        data.ci2_upper = J(num_zeval, K, .)
        data.c_check_bs = J(K, 1, .)
    }

    return(est)
}

end
