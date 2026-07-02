mata:

// =============================================================================
// didhetero_catt_core.mata
// CATT core estimation (doubly robust point estimates)
//
// Provides:
//   - didhetero_catt_core()  // Core CATT estimator for specified (g,t) pairs
//
// Requires:
//   - didhetero_types.mata             (DidHeteroData, DidHeteroCattResult)
//   - didhetero_intermediate.mata      (didhetero_intermediate_vars)
//   - didhetero_lpr.mata               (didhetero_lpr)
//   - didhetero_bwselect_lp.mata       (_didhetero_lpbwselect_mse)
//   - didhetero_stage23.mata           (_didhetero_construct_A, _didhetero_construct_B)
//
// Paper reference: Section 3.1, three-step estimation procedure
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_catt_core()
//
// Core estimation function for CATT. Executes the full estimation loop for
// specified (g,t) pairs with given bandwidth.
//
// Parameters:
//   data      - DidHeteroData struct containing:
//               n (sample size), Z (covariate), zeval (evaluation points),
//               num_zeval (number of evaluation points), and other data arrays
//   gps_mat   - Generalized propensity score estimates from first-stage estimation
//   or_mat    - Outcome regression estimates from first-stage estimation
//   gteval    - K x 2 matrix of (group, time) pairs to estimate
//   bw_in     - scalar or K x 1 vector of bandwidth values
//   bwselect  - bandwidth selection method ("manual" or other)
//   porder    - polynomial order for local polynomial regression (1 or 2)
//   kernel    - kernel type ("gau" for Gaussian, "epa" for Epanechnikov)
//
// Returns:
//   DidHeteroCattResult struct containing:
//   - catt_est: num_zeval x K matrix of CATT point estimates
//   - A_g_t, B_g_t: influence function matrices for each (g,t)
//   - G_g: group indicator matrix
//   - mu_G_g, mu_E_g_t, mu_F_g_t: conditional mean estimates
//   - bw_vec: K x 1 vector of bandwidths used
// -----------------------------------------------------------------------------
struct DidHeteroCattResult scalar didhetero_catt_core(
    struct DidHeteroData scalar data,
    real matrix gps_mat,
    real matrix or_mat,
    real matrix gteval,
    real colvector bw_in,
    string scalar bwselect,
    real scalar porder,
    string scalar kernel)
{
    struct DidHeteroCattResult scalar result
    struct DidHeteroIntermediate scalar intermed
    real scalar n, K, num_zeval, id_gt, g1, t1, r, h_gt
    real matrix A_mat, B_mat
    real colvector mu_G_col, mu_R_col, bw_vec
    real colvector E_g_t_vec, F_g_t_vec, mu_E_bw, mu_F_bw
    real colvector mu_E_col, mu_F_col
    real colvector R_g, Y_diff, G_ig

    // Dimensions
    n = data.n
    K = rows(gteval)
    num_zeval = data.num_zeval

    // Bandwidth handling: expand scalar to vector or validate vector length
    if (bwselect == "manual") {
        if (rows(bw_in) == 1) {
            bw_vec = J(K, 1, bw_in[1])
        }
        else if (rows(bw_in) == K) {
            bw_vec = bw_in
        }
        else {
            _error("bw must be scalar or vector of length " + strofreal(K))
        }
    }
    else {
        _error("didhetero_catt_core: non-manual bwselect not supported. " +
               "Use manual bandwidth selection.")
    }

    // Initialize result struct
    result.n = n
    result.num_gteval = K
    result.num_zeval = num_zeval
    result.bw_vec = bw_vec
    result.A_g_t = J(1, K, NULL)
    result.B_g_t = J(1, K, NULL)
    result.G_g = J(n, K, 0)
    result.mu_G_g = J(num_zeval, K, .)
    result.mu_E_g_t = J(num_zeval, K, .)
    result.mu_F_g_t = J(num_zeval, K, .)
    result.catt_est = J(num_zeval, K, .)

    // Main loop over (g,t) pairs
    for (id_gt = 1; id_gt <= K; id_gt++) {
        g1 = gteval[id_gt, 1]
        t1 = gteval[id_gt, 2]
        h_gt = bw_vec[id_gt]

        // Step 1: Compute intermediate variables
        // Per Equation 2.6: R_g, Y_diff, E_g_t, F_g_t depend only on
        // Stage 1 estimates (GPS + OR), not on the LPR bandwidth h.
        {
            real matrix _tmp_G_g
            _tmp_G_g = J(n, 1, 0)
            intermed = didhetero_intermediate_vars(data, gps_mat, or_mat,
                           g1, t1, 1, data.control_group,
                           data.anticipation, _tmp_G_g)
            R_g       = intermed.R_g
            Y_diff    = intermed.Y_diff
            G_ig      = (data.G :== g1)
            E_g_t_vec = R_g :* Y_diff
            F_g_t_vec = G_ig :* Y_diff
        }

        // Step 2: Estimate mu_E and mu_F via local linear regression
        // These are auxiliary conditional means for bias correction
        mu_E_bw = _didhetero_lpbwselect_mse(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel)
        mu_E_col = didhetero_lpr(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel, mu_E_bw)

        mu_F_bw = _didhetero_lpbwselect_mse(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel)
        mu_F_col = didhetero_lpr(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel, mu_F_bw)

        result.mu_E_g_t[., id_gt] = mu_E_col
        result.mu_F_g_t[., id_gt] = mu_F_col

        // Step 3: Estimate mu_G and mu_R via local polynomial regression
        mu_G_col = didhetero_lpr(G_ig, data.Z,
                       data.zeval, porder, 0, kernel, h_gt)
        result.mu_G_g[., id_gt] = mu_G_col

        mu_R_col = didhetero_lpr(R_g, data.Z,
                       data.zeval, porder, 0, kernel, h_gt)

        // Step 4: Construct influence functions A and B
        A_mat = _didhetero_construct_A(G_ig, R_g,
                    Y_diff, mu_G_col, mu_R_col)

        B_mat = _didhetero_construct_B(A_mat, G_ig, R_g,
                    mu_G_col, mu_R_col, mu_E_col, mu_F_col)

        // Store A and B matrices (force independent copies)
        result.A_g_t[id_gt] = &(1 * A_mat)
        result.B_g_t[id_gt] = &(1 * B_mat)

        // Store group indicators
        result.G_g[., id_gt] = G_ig

        // Step 5: Compute doubly robust point estimates via LPR on A
        for (r = 1; r <= num_zeval; r++) {
            result.catt_est[r, id_gt] = didhetero_lpr(A_mat[., r], data.Z,
                                data.zeval[r], porder, 0, kernel, h_gt)
        }
    }

    return(result)
}

end
