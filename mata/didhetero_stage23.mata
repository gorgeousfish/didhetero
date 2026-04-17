mata:

// =============================================================================
// didhetero_stage23.mata
// Second/third stage estimation and influence function construction
//
// Delegates the core estimation loop (Steps 1-5) to didhetero_catt_core(),
// then runs SE computation (Step 6) and bootstrap UCB (Step 7) locally.
//
// For each (g,t) pair (via catt_core):
//   1. Construct intermediate variables
//   2. Estimate mu_E, mu_F via LLR with MSE-DPI bandwidth
//   3. Estimate mu_G, mu_R (porder, common bandwidth)
//   4. Construct A_{i,g,t}(z_r) and B_{i,g,t}(z_r)
//   5. Compute DR point estimates via LPR on A
//
// Locally (after catt_core returns):
//   6. SE and analytical UCB (per g,t pair)
//   7. Bootstrap UCB (across all g,t pairs)
//
// Functions:
//   1. _didhetero_construct_A()  - DR core term A_{i,g,t}(z_r)
//   2. _didhetero_construct_B()  - Influence function B_{i,g,t}(z_r)
//   3. didhetero_stage23()       - Main estimation orchestrator
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Section 4.2.1 (DR estimation, Stages 2 and 3)
// =============================================================================


// -----------------------------------------------------------------------------
// _didhetero_construct_A()
// Construct A_{i,g,t}(z_r) for a single (g,t) pair and all zeval points.
//
// A_i(z_r) = (G_{i,g} / mu_G(z_r) - R_{i,g,t} / mu_R(z_r)) * Y_tilde_i
//
// Properties:
//   - Treatment group (G=1, R=0): A_i = Y_tilde_i / mu_G(z_r)
//   - Comparison group (G=0, R>0): A_i = -R_i * Y_tilde_i / mu_R(z_r)
//   - Other (G=0, R=0): A_i = 0
//
// Args:
//   G_ig    - n x 1, treatment group indicator
//   R_g     - n x 1, IPW weight
//   Y_diff  - n x 1, outcome adjustment Y_tilde
//   mu_G    - R x 1, conditional mean of G at zeval
//   mu_R    - R x 1, conditional mean of R at zeval
//
// Returns:
//   n x R matrix, A_{i,g,t}(z_r)
//
// Paper ref: Section 4.2.1, Eq. 14 (A_{i,g,t} construction)
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
        // A_i(z_r) = (G_ig / mu_G[r] - R_g / mu_R[r]) * Y_diff
        A_mat[., r] = (G_ig / mu_G[r] - R_g / mu_R[r]) :* Y_diff
    }

    return(A_mat)
}


// -----------------------------------------------------------------------------
// _didhetero_construct_B()
// Construct B_{i,g,t}(z_r) for a single (g,t) pair and all zeval points.
//
// B_i(z_r) = A_i(z_r) + (mu_E(z_r) / mu_R(z_r)^2) * R_{i,g,t}
//                      - (mu_F(z_r) / mu_G(z_r)^2) * G_{i,g}
//
// This is the full influence function (Paper Eq. 14) with Hajek bias correction.
//
// Args:
//   A_mat   - n x R matrix, A_{i,g,t}(z_r)
//   G_ig    - n x 1, treatment group indicator
//   R_g     - n x 1, IPW weight
//   mu_G    - R x 1, conditional mean of G
//   mu_R    - R x 1, conditional mean of R
//   mu_E    - R x 1, conditional mean of E
//   mu_F    - R x 1, conditional mean of F
//
// Returns:
//   n x R matrix, B_{i,g,t}(z_r)
//
// Paper ref: Section 4.2.2 (B_{i,g,t} residual construction)
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
        // B_i(z_r) = A_i(z_r) + (mu_E[r] / mu_R[r]^2) * R_g
        //                      - (mu_F[r] / mu_G[r]^2) * G_ig
        B_mat[., r] = A_mat[., r] ///
            + (mu_E[r] / mu_R[r]^2) * R_g ///
            - (mu_F[r] / mu_G[r]^2) * G_ig
    }

    return(B_mat)
}


// -----------------------------------------------------------------------------
// didhetero_stage23()
// Execute the second/third stage estimation loop.
//
// Delegates the core estimation (Steps 1-5) to didhetero_catt_core():
//   1. Construct intermediate variables via didhetero_intermediate_vars()
//   2. Estimate mu_E, mu_F via LLR with MSE-DPI bandwidth
//   3. Estimate mu_G, mu_R using porder and common bandwidth bw[id_gt]
//   4. Construct A and B influence functions
//   5. Compute DR point estimates via LPR on A (per zeval point)
//
// Then runs SE and bootstrap locally:
//   6. SE and analytical UCB (per g,t pair)
//   7. Bootstrap UCB (across all g,t pairs)
//
// Args:
//   data          - DidHeteroData struct (modified in place)
//   gps_mat       - GPS result matrix from Stage 1
//   or_mat        - OR result matrix from Stage 1
//   bw            - K x 1 common bandwidth vector (already resolved)
//   bwselect      - bandwidth selection method string (unused after refactor;
//                   catt_core is always called with "manual" since BW is resolved)
//   seed          - bootstrap RNG seed; applied only when bootstrap runs
//
// Modifies in place (via data struct):
//   data.A_g_t    - 1 x K pointer array, each n x num_zeval
//   data.B_g_t    - 1 x K pointer array, each n x num_zeval
//   data.G_g      - n x K group indicator matrix
//   data.mu_G_g   - num_zeval x K conditional mean of G
//
// Returns:
//   num_zeval x K matrix of DR point estimates
//
// Paper ref: Section 4.2.1, Stages 2-3 (DR estimation pipeline)
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
    // Steps 1-5: Core CATT estimation via didhetero_catt_core()
    // Delegates intermediate vars, mu_E/F/G/R estimation, A/B construction,
    // and DR point estimates to the reentrant core function.
    // BW is already resolved (by didhetero_bw_preloop or user), so we
    // always call with bwselect="manual".
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
    // Step 6: SE and analytical UCB (Paper Section 4.2.2-4.2.3)
    // Loop over (g,t) pairs; uses B_g_t and est from catt_core results.
    // Only computed when kd0_Z is available from Stage 1.
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
            // Positivity post-guard: zero out SE and CI for zeval points where
            // the GPS (mu_G = E[G_ig | Z=z_r]) is below a minimum threshold.
            //
            // When zeval[r] lies far outside the treated cohort's support,
            // mu_G[r] \u2248 0 causes the DR influence function (B_i \u221d 1/mu_G) to
            // blow up. The sigma2_raw<0 guard in didhetero_se catches downward
            // extrapolation (SE=.) but NOT the upward case (SE=huge).
            //
            // Direct positivity check: if mu_G[r, id_gt] < eps_pos, the overlap
            // assumption fails at z_r and the CATT is not point-identified there.
            // Setting SE=. (and CI=.) correctly communicates this to the user.
            //
            // Threshold eps_pos = 1e-3: flags GPS below 0.1%, safe margin
            // above typical minimum GPS at interior zeval (~5-30%).
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
    // Step 7: Bootstrap UCB (Paper Section 4.2.4)
    // When biters > 0, run weighted bootstrap for uniform confidence bands.
    // Otherwise, set CI2 fields to missing.
    // =================================================================
    if (data.biters > 0) {
        if (seed >= 0 & seed < .) {
            rseed(seed)
        }

        // Use optimized bootstrap with precompute + batched weights
        // Semantics are identical to didhetero_bootstrap_ucb().
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
