// =============================================================================
// didhetero_bootstrap_opt.mata
// Optimized bootstrap main loop for the didhetero-stata package
//
// Uses batch weight generation and pre-computed invariants to accelerate
// the bootstrap UCB procedure. Replaces the naive per-iteration LPR calls
// with pre-computed R'KR / R'KA and batched Mammen weight draws.
//
// Dependencies:
//   - didhetero_bootstrap.mata (#13): didhetero_bootstrap_precompute(),
//     didhetero_bootstrap_iter_gt(), _didhetero_bs_quantile(),
//     didhetero_max_nonmissing()
//   - didhetero_boot.mata (#14): _didhetero_mammen_weights_batch()
//
// Paper ref: Section 4.2.4, multiplier bootstrap UCB
// =============================================================================

mata:

// -----------------------------------------------------------------------------
// _dh_catt_boot_optimized()
// Optimized bootstrap UCB with batch weight generation and progress reporting.
//
// Same semantics as didhetero_bootstrap_ucb() but:
//   1. Uses didhetero_bootstrap_precompute() for O(1) per-iteration kernel work
//   2. Uses _didhetero_mammen_weights_batch() for vectorized weight generation
//   3. Allocates only mb_sup_t (no mb_est/mb_t arrays)
//   4. Reports ETA-aware progress every 100 iterations
//
// Args:
//   A_g_t         - pointer(real matrix) rowvector, 1 x num_gteval
//   est           - real matrix, num_zeval x num_gteval point estimates
//   se            - real matrix, num_zeval x num_gteval standard errors
//   bw            - real colvector, num_gteval x 1 bandwidths
//   Z             - real colvector, n x 1 covariate values
//   zeval         - real colvector, num_zeval x 1 evaluation points
//   n             - real scalar, sample size
//   porder        - real scalar, polynomial order for LPR
//   kernel        - string scalar, kernel type ("epa" or "gau")
//   alp           - real scalar, significance level
//   biters        - real scalar, number of bootstrap iterations
//   uniformall    - real scalar, 1=uniform across all (g,t), 0=per (g,t)
//   num_gteval    - real scalar, number of (g,t) pairs
//   num_zeval     - real scalar, number of evaluation points
//   ci2_lower     - real matrix, output num_zeval x num_gteval lower CI bounds
//   ci2_upper     - real matrix, output num_zeval x num_gteval upper CI bounds
//   c_check       - real colvector, output num_gteval x 1 critical values
//   batch_size    - (optional) real scalar, weight batch size (default 100)
//   show_progress - (optional) real scalar, 1=show progress (default 1)
//
// Paper ref: Section 4.2.4, multiplier bootstrap UCB
// -----------------------------------------------------------------------------
void didhetero_boot_ucb_optimized(
    pointer(real matrix) rowvector A_g_t,
    real matrix est,
    real matrix se,
    real colvector bw,
    real colvector Z,
    real colvector zeval,
    real scalar n,
    real scalar porder,
    string scalar kernel,
    real scalar alp,
    real scalar biters,
    real scalar uniformall,
    real scalar num_gteval,
    real scalar num_zeval,
    real matrix ci2_lower,
    real matrix ci2_upper,
    real colvector c_check,
    | real scalar batch_size,
    real scalar show_progress)
{
    real scalar b, id_gt, effective_uniformall
    real scalar batch_start, batch_end, actual_batch, b_in_batch
    real scalar c_val
    real matrix mb_sup_t, weight_batch
    real colvector mb_weight, mb_sup_t_all
    struct BootPrecomp scalar precomp

    // === Default optional parameters ===
    if (args() < 18) {
        batch_size = 100
    }
    if (args() < 19) {
        show_progress = 1
    }

    // === Auto-override: uniformall forced to FALSE when num_gteval == 1 ===
    // When there is only one (g,t) pair, row-max across columns is identity,
    // so uniformall=TRUE degenerates to uniformall=FALSE. Force it explicitly.
    if (num_gteval == 1) {
        effective_uniformall = 0
    }
    else {
        effective_uniformall = uniformall
    }

    // === Pre-compute bootstrap invariants ===
    didhetero_bootstrap_precompute(A_g_t, Z, zeval, bw, porder, kernel, ///
        n, num_gteval, num_zeval, precomp)

    // === Allocate only sup-t array (no mb_est/mb_t) ===
    mb_sup_t = J(biters, num_gteval, .)

    // === Bootstrap batch loop ===
    for (batch_start = 1; batch_start <= biters; batch_start = batch_start + batch_size) {

        // Determine actual batch size (last batch may be smaller)
        batch_end = min((batch_start + batch_size - 1, biters))
        actual_batch = batch_end - batch_start + 1

        // Generate weights for entire batch at once
        weight_batch = _didhetero_mammen_weights_batch(actual_batch, n)

        // Inner loop over iterations within this batch
        for (b_in_batch = 1; b_in_batch <= actual_batch; b_in_batch++) {

            b = batch_start + b_in_batch - 1

            // Extract weight vector for this iteration (row -> column)
            mb_weight = weight_batch[b_in_batch, .]'

            // Loop over (g,t) pairs
            for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
                mb_sup_t[b, id_gt] = didhetero_bootstrap_iter_gt( ///
                    mb_weight, precomp, A_g_t[id_gt], ///
                    est[., id_gt], se[., id_gt], id_gt)
            }

            // Progress display every 100 iterations
            if (show_progress) {
                if (mod(b, 100) == 0 | b == biters) {
                    printf("Bootstrap: %g/%g (%5.1f%%)\n", ///
                        b, biters, 100 * b / biters)
                }
            }
        }
    }

    // === Critical value extraction ===
    // (Copied from didhetero_bootstrap_ucb)
    c_check = J(num_gteval, 1, .)

    if (effective_uniformall) {
        // uniformall=TRUE: row-wise max across (g,t), then single quantile
        mb_sup_t_all = J(biters, 1, .)
        for (b = 1; b <= biters; b++) {
            mb_sup_t_all[b] = didhetero_max_nonmissing(mb_sup_t[b, .]')
        }
        c_val = _didhetero_bs_quantile(mb_sup_t_all, 1 - alp)
        c_check = J(num_gteval, 1, c_val)
    }
    else {
        // uniformall=FALSE: column-wise quantile, one per (g,t)
        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
            c_check[id_gt] = _didhetero_bs_quantile(mb_sup_t[., id_gt], 1 - alp)
        }
    }

    // === Construct CI2 ===
    ci2_lower = J(num_zeval, num_gteval, .)
    ci2_upper = J(num_zeval, num_gteval, .)

    for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
        if (c_check[id_gt] < .) {
            ci2_lower[., id_gt] = est[., id_gt] - c_check[id_gt] * se[., id_gt]
            ci2_upper[., id_gt] = est[., id_gt] + c_check[id_gt] * se[., id_gt]
        }
    }
}

end
