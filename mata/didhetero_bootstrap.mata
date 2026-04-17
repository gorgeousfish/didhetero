mata:

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
// didhetero_max_nonmissing()
// Compute maximum of non-missing values in a vector.
// Returns missing (.) if all values are missing.
// -----------------------------------------------------------------------------
real scalar didhetero_max_nonmissing(real colvector x)
{
    real colvector x_valid

    x_valid = select(x, x :< .)
    if (rows(x_valid) == 0) return(.)
    return(max(x_valid))
}

// -----------------------------------------------------------------------------
// didhetero_lpr_point()
// Compute weighted local polynomial regression estimate at a single evaluation
// point. Thin wrapper for bootstrap inner loop efficiency.
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
// Pre-compute kernel values and design matrices for bootstrap loop efficiency.
// Quantities invariant across bootstrap iterations are computed once to avoid
// redundant evaluations.
//
// Arguments:
//   A_g_t       - pointer(real matrix) rowvector, influence function matrices
//   Z           - real colvector, continuous covariate values
//   zeval       - real colvector, evaluation points
//   bw          - real colvector, bandwidths per (g,t) pair
//   porder      - real scalar, polynomial order
//   kernel      - string scalar, kernel type ("epa" or "gau")
//   n           - real scalar, sample size
//   num_gteval  - real scalar, number of (g,t) pairs
//   num_zeval   - real scalar, number of evaluation points
//   precomp     - struct BootPrecomp scalar, output structure (modified in place)
//
// Returns: void
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
    // RtKR and RtKA are not precomputed because bootstrap weights vary by
    // iteration, making these products iteration-dependent.
    precomp.RtKR        = J(num_gteval, num_zeval, NULL)
    precomp.RtKA        = J(num_gteval, num_zeval, NULL)

    for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
        h = bw[id_gt]

        for (r = 1; r <= num_zeval; r++) {

            // Compute differences and kernel weights
            diff = (Z :- zeval[r])
            K_vec = didhetero_kernel_eval(diff / h, kernel)

            // Build design matrix: column j contains (Z_i - z_r)^(j-1)
            R_mat = J(n, porder + 1, .)
            for (j = 0; j <= porder; j++) {
                R_mat[., j + 1] = diff :^ j
            }

            // Store precomputed values
            precomp.kernel_vals[id_gt, r] = &(K_vec[., .])
            precomp.design_mats[id_gt, r] = &(R_mat[., .])
        }
    }
}


// -----------------------------------------------------------------------------
// didhetero_bootstrap_iter_gt()
// Execute one bootstrap iteration for a single (g,t) pair.
// Computes studentized supremum t-statistic across evaluation points.
//
// Arguments:
//   mb_weight  - real colvector, Mammen weights for this iteration
//   precomp    - struct BootPrecomp scalar, pre-computed kernel and design matrices
//   A_gt_ptr   - pointer(real matrix) scalar, influence function matrix
//   est_col    - real colvector, point estimates
//   se_col     - real colvector, standard errors
//   id_gt      - real scalar, index of current (g,t) pair
//
// Returns: real scalar, supremum t-statistic for this (g,t) pair
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

        // Compute effective weights and weighted least squares
        w_eff = mb_weight :* (*precomp.kernel_vals[id_gt, r])
        R_mat = *precomp.design_mats[id_gt, r]
        RtWR = cross(R_mat, w_eff, R_mat)
        A_r = (*A_gt_ptr)[., r]
        eff_index = didhetero_selectindex(w_eff :!= 0)

        if (length(eff_index) == 0) {
            mb_est_r = .
        }
        else if (length(eff_index) == 1) {
            // Single effective observation: local fit equals observed value
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

        // Compute studentized t-statistic
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
// Construct uniform confidence bands via wild bootstrap.
// Computes bootstrap critical values and simultaneous confidence bands
// for heterogeneous treatment effect estimates.
//
// Arguments:
//   A_g_t       - pointer(real matrix) rowvector, influence function matrices
//   est         - real matrix, point estimates
//   se          - real matrix, standard errors
//   bw          - real colvector, bandwidths
//   Z           - real colvector, covariate values
//   zeval       - real colvector, evaluation points
//   n           - real scalar, sample size
//   porder      - real scalar, polynomial order
//   kernel      - string scalar, kernel type
//   alp         - real scalar, significance level
//   biters      - real scalar, bootstrap iterations
//   uniformall  - real scalar, uniform across all (g,t) pairs flag
//   num_gteval  - real scalar, number of (g,t) pairs
//   num_zeval   - real scalar, number of evaluation points
//   ci2_lower   - real matrix, output lower confidence bounds
//   ci2_upper   - real matrix, output upper confidence bounds
//   c_check     - real colvector, output critical values
//
// Returns: void
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

    // Override uniformall when only one (g,t) pair exists
    if (num_gteval == 1) {
        effective_uniformall = 0
    }
    else {
        effective_uniformall = uniformall
    }

    // Initialize storage matrices
    mb_est   = J(biters, num_zeval * num_gteval, .)
    mb_t     = J(biters, num_zeval * num_gteval, .)
    mb_sup_t = J(biters, num_gteval, .)

    // Bootstrap iterations
    for (b = 1; b <= biters; b++) {

        // Generate Mammen weights (shared across all (g,t) pairs)
        mb_weight = _didhetero_mammen_weights(n)

        // Loop over (g,t) pairs
        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {

            est_temp = est[., id_gt]
            se_temp  = se[., id_gt]

            for (r = 1; r <= num_zeval; r++) {

                col_idx = (id_gt - 1) * num_zeval + r

                // Compute bootstrap estimate via weighted local polynomial regression
                mb_est[b, col_idx] = didhetero_lpr_point(
                    (*A_g_t[id_gt])[., r], Z, zeval[r], porder, 0, kernel,
                    bw[id_gt], mb_weight
                )

                // Compute studentized t-statistic
                if (se_temp[r] != . & se_temp[r] > 0) {
                    mb_t[b, col_idx] = abs(mb_est[b, col_idx] - est_temp[r]) / se_temp[r]
                }
            }

            // Compute supremum t-statistic across evaluation points
            t_vals = J(num_zeval, 1, .)
            for (r = 1; r <= num_zeval; r++) {
                col_idx = (id_gt - 1) * num_zeval + r
                t_vals[r] = mb_t[b, col_idx]
            }
            mb_sup_t[b, id_gt] = didhetero_max_nonmissing(t_vals)
        }

        // Display progress
        if (mod(b, 50) == 0 | b == biters) {
            printf("Bootstrap iteration %g of %g\n", b, biters)
        }
    }

    // Extract critical values
    c_check = J(num_gteval, 1, .)

    if (effective_uniformall) {
        // Uniform across all (g,t) pairs: max over pairs then quantile
        mb_sup_t_all = J(biters, 1, .)
        for (b = 1; b <= biters; b++) {
            mb_sup_t_all[b] = didhetero_max_nonmissing(mb_sup_t[b, .]')
        }
        c_val = _didhetero_bs_quantile(mb_sup_t_all, 1 - alp)
        c_check = J(num_gteval, 1, c_val)
    }
    else {
        // Pair-specific critical values
        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
            c_check[id_gt] = _didhetero_bs_quantile(mb_sup_t[., id_gt], 1 - alp)
        }
    }

    // Construct confidence bands
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
