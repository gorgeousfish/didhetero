mata:

// =============================================================================
// didhetero_bootstrap.mata
// Bootstrap UCB helper functions for the didhetero-stata package
//
// Functions:
//   1. didhetero_max_nonmissing()   - Max ignoring missing values
//   2. didhetero_lpr_point()        - Single-point weighted LPR wrapper
//   3. _didhetero_bs_mammen_weights() - Bootstrap-specific Mammen weights
//   4. _didhetero_bs_quantile()     - Bootstrap-specific quantile
//   5. didhetero_bootstrap_ucb()    - Main bootstrap UCB procedure
//
//
// Paper ref: Section 4.2.4, multiplier bootstrap UCB
// =============================================================================

// NOTE: Mammen weights generator unified in didhetero_boot.mata
// All bootstrap routines should call _didhetero_mammen_weights(n).

// Bootstrap-specific quantile with missing value handling
real scalar _didhetero_bs_quantile(real colvector x, real scalar prob)
{
    real colvector x_valid
    real scalar B, h, lo, hi

    x_valid = select(x, x :< .)
    if (rows(x_valid) == 0) return(.)

    _sort(x_valid, 1)
    B = rows(x_valid)

    if (B == 1) return(x_valid[1])

    h  = (B - 1) * prob + 1
    lo = floor(h)
    hi = ceil(h)

    if (lo < 1) lo = 1
    if (hi > B) hi = B

    if (lo == hi) return(x_valid[lo])
    return(x_valid[lo] + (h - lo) * (x_valid[hi] - x_valid[lo]))
}

// -----------------------------------------------------------------------------
// didhetero_max_nonmissing
// Compute maximum of non-missing values in a vector.
// Returns . if all values are missing.
//
// Maximum ignoring missing values.
// Returns . (Stata missing) if all values are missing.
// -----------------------------------------------------------------------------
real scalar didhetero_max_nonmissing(real colvector x)
{
    real colvector x_valid

    x_valid = select(x, x :< .)
    if (rows(x_valid) == 0) return(.)
    return(max(x_valid))
}

// -----------------------------------------------------------------------------
// didhetero_lpr_point
// Compute weighted LPR estimate at a single evaluation point.
// Thin wrapper around didhetero_lpr for single-point evaluation
// in the bootstrap inner loop.
//
// Paper ref: Section 4.2.1, LPR at single evaluation point
// -----------------------------------------------------------------------------
real scalar didhetero_lpr_point(
    real colvector y,
    real colvector x,
    real scalar eval_pt,
    real scalar p,
    real scalar deriv,
    string scalar kernel,
    real scalar h,
    real colvector weight)
{
    return(didhetero_lpr(y, x, eval_pt, p, deriv, kernel, h, weight))
}


// -----------------------------------------------------------------------------
// didhetero_bootstrap_precompute()
// Pre-compute invariants for bootstrap loop optimization.
// Called once before bootstrap loop. Stores kernel values, design matrices,
// R'diag(K)R, and R'diag(K)A for all (g,t,r) combinations.
//
// These quantities do not change across bootstrap iterations b, so computing
// them once avoids O(B × K × R) redundant kernel/design matrix evaluations.
//
// Args:
//   A_g_t       - pointer(real matrix) rowvector, 1 x num_gteval
//                 Each pointer -> n x num_zeval matrix of DR scores
//   Z           - real colvector, n x 1, continuous covariate values
//   zeval       - real colvector, num_zeval x 1, evaluation points
//   bw          - real colvector, num_gteval x 1, bandwidths per (g,t)
//   porder      - real scalar, polynomial order (typically 2)
//   kernel      - string scalar, kernel type ("epa" or "gau")
//   n           - real scalar, sample size
//   num_gteval  - real scalar, number of (g,t) pairs
//   num_zeval   - real scalar, number of evaluation points
//   precomp     - struct BootPrecomp scalar, output (modified in place)
//
// Returns: void (modifies precomp in place)
//
// Paper ref: Section 4.2.4 (Bootstrap UCB construction)
// -----------------------------------------------------------------------------
void didhetero_bootstrap_precompute(
    pointer(real matrix) rowvector A_g_t,
    real colvector Z,
    real colvector zeval,
    real colvector bw,
    real scalar porder,
    string scalar kernel,
    real scalar n,
    real scalar num_gteval,
    real scalar num_zeval,
    struct BootPrecomp scalar precomp)
{
    real scalar id_gt, r, j
    real scalar h
    real colvector diff, K_vec, A_r
    real matrix R_mat

    // Store dimensions for validation
    precomp.num_gteval = num_gteval
    precomp.num_zeval  = num_zeval
    precomp.n          = n
    precomp.porder     = porder

    // Allocate pointer matrices (num_gteval × num_zeval)
    precomp.kernel_vals = J(num_gteval, num_zeval, NULL)
    precomp.design_mats = J(num_gteval, num_zeval, NULL)
    // Note: RtKR/RtKA are not precomputed anymore because bootstrap uses
    // iteration-specific weights V* that enter as diag(V*) in R' diag(V*) diag(K) R
    // and R' diag(V*) diag(K) A. Since V* varies by iteration, the products
    // cannot be reused across iterations. Keep fields as NULL to signal unused.
    precomp.RtKR        = J(num_gteval, num_zeval, NULL)
    precomp.RtKA        = J(num_gteval, num_zeval, NULL)

    for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
        h = bw[id_gt]

        for (r = 1; r <= num_zeval; r++) {

            // Differences: diff_i = (Z_i - z_r)
            diff = (Z :- zeval[r])

            // Kernel values: K((Z_i - z_r)/h), n × 1
            // Use standardized distance only for kernel, match didhetero_lpr()
            K_vec = didhetero_kernel_eval(diff / h, kernel)

            // Design matrix R_{z_r}: n × (p+1) built on unscaled diff
            // Column j+1 = (Z_i - z_r)^j, j = 0, ..., porder
            R_mat = J(n, porder + 1, .)
            for (j = 0; j <= porder; j++) {
                R_mat[., j + 1] = diff :^ j
            }

            // Store copies via pointer to expression (creates new objects)
            precomp.kernel_vals[id_gt, r] = &(K_vec[., .])
            precomp.design_mats[id_gt, r] = &(R_mat[., .])

            // Do NOT precompute RtKR / RtKA; weights vary at each bootstrap
            // iteration, so these matrices would need to be recomputed anyway.
            // Leave precomp.RtKR and precomp.RtKA entries as NULL.
        }
    }
}


// -----------------------------------------------------------------------------
// didhetero_bootstrap_iter_gt()
// Optimized bootstrap iteration for one (g,t) pair.
// Uses pre-computed invariants; only updates weights per iteration.
//
// For each evaluation point r:
//   1. Compute effective weights: w_eff = mb_weight .* K_vec
//   2. Compute R'W*R = cross(R_mat, w_eff, R_mat)
//   3. Compute R'W*A = cross(R_mat, w_eff, A_r)
//   4. Solve WLS: beta = lusolve(R'W*R, R'W*A)
//   5. Extract intercept: mb_est_r = beta[1]
//   6. Compute t-stat: |mb_est_r - est[r]| / se[r]
//   7. Update sup-t (max over r)
//
// Args:
//   mb_weight  - real colvector, n × 1, Mammen weights for this iteration
//   precomp    - struct BootPrecomp scalar, pre-computed invariants
//   A_gt_ptr   - pointer(real matrix) scalar, -> n × num_zeval DR scores
//   est_col    - real colvector, num_zeval × 1, point estimates DR(z_r)
//   se_col     - real colvector, num_zeval × 1, standard errors SE(z_r)
//   id_gt      - real scalar, index of current (g,t) pair (1-based)
//
// Returns: real scalar, sup-t statistic for this (g,t) pair
//
// Paper ref: Section 4.2.4 (Bootstrap iteration)
// -----------------------------------------------------------------------------
real scalar didhetero_bootstrap_iter_gt(
    real colvector mb_weight,
    struct BootPrecomp scalar precomp,
    pointer(real matrix) scalar A_gt_ptr,
    real colvector est_col,
    real colvector se_col,
    real scalar id_gt)
{
    real scalar r, sup_t, mb_est_r, t_r
    real colvector w_eff, A_r, A_sub, RtWA, beta, w_sub
    real matrix R_mat, RtWR, R_sub, RtWR_inv
    real rowvector eff_index

    sup_t = 0

    for (r = 1; r <= precomp.num_zeval; r++) {

        // Effective weights: V*_i × K_h(Z_i - z_r), element-wise
        w_eff = mb_weight :* (*precomp.kernel_vals[id_gt, r])

        // R'W*R: (p+1) × (p+1)
        R_mat = *precomp.design_mats[id_gt, r]
        RtWR = cross(R_mat, w_eff, R_mat)

        // R'W*A: (p+1) × 1
        A_r = (*A_gt_ptr)[., r]
        eff_index = didhetero_selectindex(w_eff :!= 0)

        if (length(eff_index) == 0) {
            mb_est_r = .
        }
        else if (length(eff_index) == 1) {
            // Mirror didhetero_lpr(): one effective observation implies the
            // local fit collapses to the unique observed score value.
            mb_est_r = A_r[eff_index]
        }
        else {
            R_sub = R_mat[eff_index', .]
            A_sub = A_r[eff_index']
            w_sub = w_eff[eff_index']

            RtWR = cross(R_sub, w_sub, R_sub)
            RtWR_inv = cholinv(RtWR)

            if (hasmissing(RtWR_inv)) {
                mb_est_r = .
            }
            else {
                RtWA = cross(R_sub, w_sub, A_sub)
                beta = RtWR_inv * RtWA
                mb_est_r = beta[1]
            }
        }

        // t statistic: |DR*_b(z_r) - DR(z_r)| / SE(z_r)
        // Paper Eq. 23: sup-t statistic for bootstrap UCB
        if (mb_est_r != . & se_col[r] != . & se_col[r] > 0) {
            t_r = abs(mb_est_r - est_col[r]) / se_col[r]
            if (t_r > sup_t) {
                sup_t = t_r
            }
        }
    }

    return(sup_t)
}


// -----------------------------------------------------------------------------
// didhetero_bootstrap_ucb()
// Main bootstrap uniform confidence band (UCB) procedure.
//
// Computes bootstrap critical values and constructs simultaneous confidence
// bands for heterogeneous treatment effect estimates across evaluation points
// and (g,t) pairs. Uses Mammen wild bootstrap with studentized sup-t statistic.
//
// Args:
//   A_g_t       - pointer(real matrix) rowvector, 1 x num_gteval pointer array,
//                 each pointing to an n x num_zeval matrix of influence functions
//   est         - real matrix, num_zeval x num_gteval point estimates
//   se          - real matrix, num_zeval x num_gteval standard errors
//   bw          - real colvector, num_gteval x 1 bandwidths
//   Z           - real colvector, n x 1 covariate values
//   zeval       - real colvector, num_zeval x 1 evaluation points
//   n           - real scalar, sample size
//   porder      - real scalar, polynomial order for LPR
//   ci2_lower   - real matrix, output num_zeval x num_gteval lower CI bounds
//   ci2_upper   - real matrix, output num_zeval x num_gteval upper CI bounds
//   c_check     - real colvector, output num_gteval x 1 critical values
//
// Paper ref: Section 4.2.4 Eq. 23 (Bootstrap UCB construction)
// -----------------------------------------------------------------------------
void didhetero_bootstrap_ucb(
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
    real colvector c_check)
{
    real scalar b, id_gt, r, col_idx, c_val, effective_uniformall
    real matrix mb_est, mb_t, mb_sup_t
    real colvector mb_weight, est_temp, se_temp, t_vals, mb_sup_t_all

    // === Auto-override: uniformall forced to FALSE when num_gteval == 1 ===
    // When there is only one (g,t) pair, row-max across columns is identity,
    // so uniformall=TRUE degenerates to uniformall=FALSE. Force it explicitly.
    if (num_gteval == 1) {
        effective_uniformall = 0
    }
    else {
        effective_uniformall = uniformall
    }

    // === Initialization ===
    // Flattened 2D arrays: biters x (num_zeval * num_gteval)
    mb_est   = J(biters, num_zeval * num_gteval, .)
    mb_t     = J(biters, num_zeval * num_gteval, .)
    mb_sup_t = J(biters, num_gteval, .)

    // === Bootstrap loop ===
    for (b = 1; b <= biters; b++) {

        // Generate Mammen weights - ALL (g,t) pairs share same weights per iteration
        // Unified Mammen weights generator
        mb_weight = _didhetero_mammen_weights(n)

        // Loop over (g,t) pairs
        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {

            est_temp = est[., id_gt]
            se_temp  = se[., id_gt]

            // Loop over evaluation points
            for (r = 1; r <= num_zeval; r++) {

                col_idx = (id_gt - 1) * num_zeval + r

                // Weighted LPR on A_{i,g,t}(z_r)
                mb_est[b, col_idx] = didhetero_lpr_point(
                    (*A_g_t[id_gt])[., r], Z, zeval[r], porder, 0, kernel,
                    bw[id_gt], mb_weight
                )

                // Studentized t-statistic
                if (se_temp[r] != . & se_temp[r] > 0) {
                    mb_t[b, col_idx] = abs(mb_est[b, col_idx] - est_temp[r]) / se_temp[r]
                }
            }

            // Sup-t statistic: max over evaluation points for this (g,t)
            t_vals = J(num_zeval, 1, .)
            for (r = 1; r <= num_zeval; r++) {
                col_idx = (id_gt - 1) * num_zeval + r
                t_vals[r] = mb_t[b, col_idx]
            }
            mb_sup_t[b, id_gt] = didhetero_max_nonmissing(t_vals)
        }

        // Progress display every 50 iterations or at the end
        if (mod(b, 50) == 0 | b == biters) {
            printf("Bootstrap iteration %g of %g\n", b, biters)
        }
    }

    // === Critical value extraction ===
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
