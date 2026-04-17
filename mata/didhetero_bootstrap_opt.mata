// =============================================================================
// Optimized Bootstrap Implementation for Uniform Confidence Bands
//
// This module implements a computationally efficient multiplier bootstrap
// procedure for constructing uniform confidence bands. The algorithm employs
// batch weight generation and pre-computed kernel invariants to reduce
// per-iteration computational overhead.
//
// Key optimizations:
//   1. Pre-computation of R'KR and R'KA matrices eliminates repeated kernel
//      evaluations across bootstrap iterations
//   2. Vectorized generation of Mammen multipliers in batches reduces
//      function call overhead
//   3. Memory-efficient storage using only the sup-t statistic matrix
// =============================================================================

mata:

// -----------------------------------------------------------------------------
// didhetero_boot_ucb_optimized()
//
// Computes uniform confidence bands via multiplier bootstrap with optimized
// memory allocation and batch processing. This implementation reduces
// computational overhead through pre-computed invariants and vectorized
// weight generation.
//
// Parameters:
//   A_g_t         - pointer rowvector to influence matrices for each (g,t)
//   est           - point estimates (num_zeval x num_gteval)
//   se            - standard errors (num_zeval x num_gteval)
//   bw            - bandwidth vector (num_gteval x 1)
//   Z             - covariate values (n x 1)
//   zeval         - evaluation points (num_zeval x 1)
//   n             - sample size
//   porder        - polynomial order for local polynomial regression
//   kernel        - kernel type ("epa" or "gau")
//   alp           - significance level
//   biters        - number of bootstrap iterations
//   uniformall    - flag for uniform band type (1=all (g,t), 0=per (g,t))
//   num_gteval    - number of (g,t) pairs
//   num_zeval     - number of evaluation points
//   ci2_lower     - output: lower confidence bounds (num_zeval x num_gteval)
//   ci2_upper     - output: upper confidence bounds (num_zeval x num_gteval)
//   c_check       - output: critical values (num_gteval x 1)
//   batch_size    - optional: batch size for weight generation (default 100)
//   show_progress - optional: display progress indicator (default 1)
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

    // Set default values for optional parameters
    if (args() < 18) {
        batch_size = 100
    }
    if (args() < 19) {
        show_progress = 1
    }

    // Override uniformall when only one (g,t) pair exists
    if (num_gteval == 1) {
        effective_uniformall = 0
    }
    else {
        effective_uniformall = uniformall
    }

    // Pre-compute kernel invariants for efficient iteration
    didhetero_bootstrap_precompute(A_g_t, Z, zeval, bw, porder, kernel, ///
        n, num_gteval, num_zeval, precomp)

    // Allocate storage for supremum t-statistics
    mb_sup_t = J(biters, num_gteval, .)

    // Main bootstrap loop with batch processing
    for (batch_start = 1; batch_start <= biters; batch_start = batch_start + batch_size) {

        // Adjust batch size for final iteration
        batch_end = min((batch_start + batch_size - 1, biters))
        actual_batch = batch_end - batch_start + 1

        // Generate Mammen multipliers for current batch
        weight_batch = _didhetero_mammen_weights_batch(actual_batch, n)

        // Process each bootstrap iteration in batch
        for (b_in_batch = 1; b_in_batch <= actual_batch; b_in_batch++) {

            b = batch_start + b_in_batch - 1

            // Extract weight vector for current iteration
            mb_weight = weight_batch[b_in_batch, .]'

            // Compute sup-t statistic for each (g,t) pair
            for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
                mb_sup_t[b, id_gt] = didhetero_bootstrap_iter_gt( ///
                    mb_weight, precomp, A_g_t[id_gt], ///
                    est[., id_gt], se[., id_gt], id_gt)
            }

            // Display progress at specified intervals
            if (show_progress) {
                if (mod(b, 100) == 0 | b == biters) {
                    printf("Bootstrap: %g/%g (%5.1f%%)\n", ///
                        b, biters, 100 * b / biters)
                }
            }
        }
    }

    // Extract critical values from bootstrap distribution
    c_check = J(num_gteval, 1, .)

    if (effective_uniformall) {
        // Uniform band: max across all (g,t), then single quantile
        mb_sup_t_all = J(biters, 1, .)
        for (b = 1; b <= biters; b++) {
            mb_sup_t_all[b] = didhetero_max_nonmissing(mb_sup_t[b, .]')
        }
        c_val = _didhetero_bs_quantile(mb_sup_t_all, 1 - alp)
        c_check = J(num_gteval, 1, c_val)
    }
    else {
        // Pointwise band: separate quantile for each (g,t)
        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
            c_check[id_gt] = _didhetero_bs_quantile(mb_sup_t[., id_gt], 1 - alp)
        }
    }

    // Construct uniform confidence intervals
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
