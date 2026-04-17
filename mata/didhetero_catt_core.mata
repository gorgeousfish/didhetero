mata:

// =============================================================================
// didhetero_catt_core.mata
// Reentrant CATT core estimation function
//
// Extracts the core CATT estimation loop from didhetero_stage23() into a
// standalone function that can be called by both the main pipeline and
// aggte Pass 2.
//
// The function runs the per-(g,t) estimation loop:
//   1. Construct intermediate variables
//   2. Estimate mu_G, mu_R via LPR
//   3. Construct A and B influence functions
//   4. Compute DR point estimates
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Section 4 (CATT estimation core)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_catt_core()
//
// Reentrant CATT estimation core. Runs the full estimation loop for a given
// set of (g,t) pairs with specified bandwidth.
//
// This function is called by:
//   1. didhetero_stage23() — main pipeline (original behavior, future refactor)
//   2. didhetero_aggte_pass2() — Pass 2 re-estimation with h_agg
//
// Parameters:
//   data      - DidHeteroData struct (read-only for core arrays;
//               intermediate vars are constructed internally)
//   gps_mat   - GPS result matrix from Stage 1
//   or_mat    - OR result matrix from Stage 1
//   gteval    - K x 2 matrix of (g,t) pairs to estimate
//   bw_in     - scalar or K x 1 vector of bandwidths
//   bwselect  - bandwidth selection method ("manual", "IMSE2", etc.)
//   porder    - polynomial order (1 or 2)
//   kernel    - kernel type ("gau" or "epa")
//
// Returns:
//   DidHeteroCattResult struct with B_g_t, G_g, mu_G_g, catt_est, bw_vec
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
    real scalar n, K, num_zeval, id_gt, g1, t1, r, h_gt
    real matrix A_mat, B_mat, G_g_local
    real colvector mu_G_col, mu_R_col, bw_vec
    real colvector E_g_t_vec, F_g_t_vec, mu_E_bw, mu_F_bw
    real colvector mu_E_col, mu_F_col
    struct DidHeteroIntermediate scalar intermed

    // --- Dimensions ---
    n = data.n
    K = rows(gteval)
    num_zeval = data.num_zeval

    // --- Bandwidth handling ---
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
        // For non-manual BW selection, we need to run the BW pre-loop.
        // This path is used by the main pipeline (didhetero_stage23).
        // For Pass 2, bwselect is always "manual".
        _error("didhetero_catt_core: non-manual bwselect not yet supported " +
               "in core function. Use didhetero_bw_preloop() + manual BW.")
    }

    // --- Initialize result struct ---
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

    // Local G_g matrix — NOT data.G_g.
    // didhetero_intermediate_vars() writes to the G_g matrix passed to it,
    // and we must not modify the original data struct's G_g.
    G_g_local = J(n, K, 0)

    // =====================================================================
    // Main loop over (g,t) pairs
    // =====================================================================
    for (id_gt = 1; id_gt <= K; id_gt++) {
        g1 = gteval[id_gt, 1]
        t1 = gteval[id_gt, 2]
        h_gt = bw_vec[id_gt]

        // =================================================================
        // Step 1: Intermediate variables
        // =================================================================
        intermed = didhetero_intermediate_vars(data, gps_mat, or_mat,
                       g1, t1, id_gt, data.control_group,
                       data.anticipation, G_g_local)

        // =================================================================
        // Step 2: mu_E, mu_F conditional estimation
        // For manual BW (Pass 2), always estimate these.
        // Uses p=1 LLR with MSE-DPI bandwidth (same as manual branch
        // in didhetero_stage23).
        // =================================================================
        E_g_t_vec = intermed.R_g :* intermed.Y_diff
        F_g_t_vec = intermed.G_ig :* intermed.Y_diff

        mu_E_bw = _didhetero_lpbwselect_mse(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel)
        mu_E_col = didhetero_lpr(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel, mu_E_bw)

        mu_F_bw = _didhetero_lpbwselect_mse(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel)
        mu_F_col = didhetero_lpr(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, kernel, mu_F_bw)

        // Store mu_E and mu_F in result struct
        result.mu_E_g_t[., id_gt] = mu_E_col
        result.mu_F_g_t[., id_gt] = mu_F_col

        // =================================================================
        // Step 3: mu_G, mu_R estimation with porder and common bandwidth
        // =================================================================
        mu_G_col = didhetero_lpr(G_g_local[., id_gt], data.Z,
                       data.zeval, porder, 0, kernel, h_gt)
        result.mu_G_g[., id_gt] = mu_G_col

        mu_R_col = didhetero_lpr(intermed.R_g, data.Z,
                       data.zeval, porder, 0, kernel, h_gt)

        // =================================================================
        // Step 4: Construct A and B
        // =================================================================
        A_mat = _didhetero_construct_A(intermed.G_ig, intermed.R_g,
                    intermed.Y_diff, mu_G_col, mu_R_col)

        B_mat = _didhetero_construct_B(A_mat, intermed.G_ig, intermed.R_g,
                    mu_G_col, mu_R_col, mu_E_col, mu_F_col)

        // Store A and B — force independent copies via expression to avoid
        // pointer aliasing across loop iterations.
        result.A_g_t[id_gt] = &(1 * A_mat)
        result.B_g_t[id_gt] = &(1 * B_mat)

        // Store G_g
        result.G_g[., id_gt] = G_g_local[., id_gt]

        // =================================================================
        // Step 5: DR point estimate via LPR on A
        // Each z_r uses a different A column (A depends on z_r)
        // =================================================================
        for (r = 1; r <= num_zeval; r++) {
            result.catt_est[r, id_gt] = didhetero_lpr(A_mat[., r], data.Z,
                                data.zeval[r], porder, 0, kernel, h_gt)
        }
    }

    return(result)
}

end
