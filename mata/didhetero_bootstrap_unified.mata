// =============================================================================
// didhetero_bootstrap_unified.mata
// Unified bootstrap inference module -- Application Layer
//
// Provides:
//   - _didhetero_mammen_weights()       // Single draw of Mammen weights
//   - _didhetero_mammen_weights_batch() // Batch Mammen weight generation
//   - _didhetero_bs_quantile()          // Bootstrap quantile (delegates to engine)
//   - didhetero_max_nonmissing()        // Max ignoring missing (delegates)
//   - didhetero_lpr_point()             // Single-point LPR for bootstrap iter
//   - didhetero_bootstrap_precompute()  // Precompute kernel/design matrices
//   - didhetero_bootstrap_iter_gt()     // Single bootstrap iteration per (g,t)
//   - didhetero_bootstrap_ucb()         // Main bootstrap UCB interface
//   - didhetero_boot_ucb_optimized()    // Optimized vectorized bootstrap
//
// Requires:
//   - didhetero_bootstrap_engine.mata   (dh_boot_generate_weights, dh_boot_quantile,
//                                        dh_boot_compute_sup, dh_boot_max_nonmissing)
//   - didhetero_types.mata              (DidHeteroData, BootPrecomp)
//   - didhetero_kernel.mata             (didhetero_kernel_eval)
//   - didhetero_lpr.mata                (didhetero_polynomial)
//
// Paper reference: Section 3.1.4, multiplier bootstrap for uniform bands
// =============================================================================

mata:

// =============================================================================
// Section 1: Mammen Weight Generation
// Delegates to dh_boot_generate_weights() from bootstrap_engine.
// Retains original function signatures for backward compatibility.
// =============================================================================

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights()
//
// Generate Mammen (1993) two-point distribution weights.
// Returns an n x 1 vector of i.i.d. draws satisfying E[V*] = 1, Var[V*] = 1.
//
// Delegates to: dh_boot_generate_weights(n, 1, "mammen")
// -----------------------------------------------------------------------------
real colvector _didhetero_mammen_weights(real scalar n)
{
    // Delegate to engine: generate n x 1 weight matrix, return as colvector
    return(dh_boot_generate_weights(n, 1, "mammen")[., 1])
}

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights_batch()
//
// Generate Mammen weights for all bootstrap iterations.
// Returns a B x n matrix where each row is an independent draw of n weights.
//
// Delegates to: dh_boot_generate_weights(n, B, "mammen")' (transposed)
// -----------------------------------------------------------------------------
real matrix _didhetero_mammen_weights_batch(real scalar B, real scalar n)
{
    // Engine returns n x B; transpose to B x n for legacy API compatibility
    return(dh_boot_generate_weights(n, B, "mammen")')
}


// =============================================================================
// Section 2: Internal Utilities
// Delegates to engine primitives; retains original function names.
// =============================================================================

// -----------------------------------------------------------------------------
// _didhetero_bs_quantile()
// Bootstrap-specific quantile with missing value handling.
// Delegates to: dh_boot_quantile()
// -----------------------------------------------------------------------------
real scalar _didhetero_bs_quantile(real colvector x, real scalar prob)
{
    return(dh_boot_quantile(x, prob))
}

// -----------------------------------------------------------------------------
// didhetero_max_nonmissing()
// Compute maximum of non-missing values in a vector.
// Delegates to: dh_boot_max_nonmissing()
// -----------------------------------------------------------------------------
real scalar didhetero_max_nonmissing(real colvector x)
{
    return(dh_boot_max_nonmissing(x))
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


// =============================================================================
// Section 3: Precomputation
// =============================================================================

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

    // Allocate pointer matrices (num_gteval x num_zeval)
    precomp.kernel_vals = J(num_gteval, num_zeval, NULL)
    precomp.design_mats = J(num_gteval, num_zeval, NULL)
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


// =============================================================================
// Section 4: Core Iteration
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_bootstrap_iter_gt()
// Execute one bootstrap iteration for a single (g,t) pair.
// Computes studentized supremum t-statistic across evaluation points.
//
// Uses precomputed R'KR Gram matrix to avoid redundant O(n*(p+1)^2) computation
// per iteration. The Gram matrix R'*diag(K)*R is iteration-invariant because
// kernel weights K_h((Z_i - z_r)/h) depend only on fixed quantities.
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
    real colvector w_eff, A_r, RtWA, beta
    real matrix R_mat, RtKR, RtKR_inv
    real rowvector eff_index

    sup_t = 0

    for (r = 1; r <= precomp.num_zeval; r++) {

        // Effective weights: Mammen * kernel (kernel is precomputed)
        w_eff = mb_weight :* (*precomp.kernel_vals[id_gt, r])
        A_r = (*A_gt_ptr)[., r]

        // Effective index: since Mammen weights are always positive (v1>0, v2>0),
        // eff_index depends solely on kernel support (invariant across iterations)
        eff_index = didhetero_selectindex((*precomp.kernel_vals[id_gt, r]) :!= 0)

        if (length(eff_index) == 0) {
            mb_est_r = .
        }
        else if (length(eff_index) == 1) {
            // Single effective observation: local fit equals observed value
            mb_est_r = A_r[eff_index]
        }
        else {
            // Bootstrap WLS: beta = (R'*diag(V*K)*R)^{-1} * R'*diag(V*K)*A
            // V (Mammen weights) MUST appear in both Gram matrix and moment vector.
            // Only R_mat and kernel_vals are precomputed (iteration-invariant);
            // the Gram matrix R'*diag(w_eff)*R must be recomputed each iteration
            // because w_eff = V * K_h changes with V.
            R_mat = *precomp.design_mats[id_gt, r]
            RtKR = cross(R_mat, w_eff, R_mat)
            RtKR_inv = cholinv(RtKR)

            if (hasmissing(RtKR_inv)) {
                mb_est_r = .
            }
            else {
                RtWA = cross(R_mat, w_eff, A_r)
                beta = RtKR_inv * RtWA
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


// =============================================================================
// Section 5: Shared Post-processing
// Delegates critical value extraction to dh_boot_compute_sup() from engine.
// =============================================================================

// -----------------------------------------------------------------------------
// _didhetero_boot_extract_crit()
// Extract critical values from bootstrap supremum t-statistic distribution.
// Delegates to: dh_boot_compute_sup()
//
// Arguments:
//   mb_sup_t            - real matrix (biters x num_gteval), sup-t statistics
//   biters              - real scalar, number of bootstrap iterations
//   num_gteval          - real scalar, number of (g,t) pairs
//   effective_uniformall - real scalar, flag for uniform band type
//   alp                 - real scalar, significance level
//   c_check             - real colvector, output critical values (modified)
//
// Returns: void
// -----------------------------------------------------------------------------
void _didhetero_boot_extract_crit(
    real matrix mb_sup_t,
    real scalar biters,
    real scalar num_gteval,
    real scalar effective_uniformall,
    real scalar alp,
    real colvector c_check)
{
    // Delegate entirely to engine's strategy-based sup computation
    c_check = dh_boot_compute_sup(mb_sup_t, effective_uniformall, alp)
}

// -----------------------------------------------------------------------------
// _didhetero_bootstrap_build_ci()
// Construct uniform confidence intervals from estimates, standard errors,
// and bootstrap critical values.
//
// Arguments:
//   est         - real matrix (num_zeval x num_gteval), point estimates
//   se          - real matrix (num_zeval x num_gteval), standard errors
//   c_check     - real colvector (num_gteval x 1), critical values
//   num_zeval   - real scalar, number of evaluation points
//   num_gteval  - real scalar, number of (g,t) pairs
//   ci2_lower   - real matrix, output lower bounds (modified)
//   ci2_upper   - real matrix, output upper bounds (modified)
//
// Returns: void
// -----------------------------------------------------------------------------
void _didhetero_bootstrap_build_ci(
    real matrix est,
    real matrix se,
    real colvector c_check,
    real scalar num_zeval,
    real scalar num_gteval,
    real matrix ci2_lower,
    real matrix ci2_upper)
{
    real scalar id_gt

    ci2_lower = J(num_zeval, num_gteval, .)
    ci2_upper = J(num_zeval, num_gteval, .)

    for (id_gt = 1; id_gt <= num_gteval; id_gt++) {
        if (c_check[id_gt] < .) {
            ci2_lower[., id_gt] = est[., id_gt] - c_check[id_gt] * se[., id_gt]
            ci2_upper[., id_gt] = est[., id_gt] + c_check[id_gt] * se[., id_gt]
        }
    }
}


// =============================================================================
// Section 5b: Vectorized Bootstrap Core (Matrix Acceleration)
// =============================================================================

// -----------------------------------------------------------------------------
// _dh_boot_vectorized_p1()
//
// Fully vectorized bootstrap for porder=1 (local linear regression).
// Eliminates the b-loop by using the analytic 2x2 WLS solution:
//   beta_0 = (s2*m0 - s1*m1) / (s0*s2 - s1^2)
//
// where s0, s1, s2, m0, m1 are computed as matrix products across all B
// iterations simultaneously.
//
// Arguments:
//   V_mat       - n x B matrix of Mammen weights (pre-generated)
//   precomp     - struct BootPrecomp scalar, pre-computed kernel/design
//   A_g_t       - pointer rowvector to influence matrices
//   est         - point estimates (num_zeval x num_gteval)
//   se          - standard errors (num_zeval x num_gteval)
//   biters      - number of bootstrap iterations B
//   num_gteval  - number of (g,t) pairs
//   num_zeval   - number of evaluation points
//   mb_sup_t    - output: biters x num_gteval sup-t statistics (modified)
//
// Returns: void
// -----------------------------------------------------------------------------
void _dh_boot_vectorized_p1(
    real matrix V_mat,
    struct BootPrecomp scalar precomp,
    pointer(real matrix) rowvector A_g_t,
    real matrix est,
    real matrix se,
    real scalar biters,
    real scalar num_gteval,
    real scalar num_zeval,
    real matrix mb_sup_t)
{
    real scalar id_gt, r, n
    real colvector K_vec, diff, A_r, KA, Kdiff, KdiffA, Kdiff2
    real rowvector s0, s1, s2, m0, m1, det, beta0_row, t_row
    real matrix t_mat
    real scalar est_r, se_r

    n = precomp.n

    for (id_gt = 1; id_gt <= num_gteval; id_gt++) {

        // Allocate t-statistic matrix: num_zeval x B
        t_mat = J(num_zeval, biters, 0)

        for (r = 1; r <= num_zeval; r++) {

            K_vec = *precomp.kernel_vals[id_gt, r]
            A_r = (*A_g_t[id_gt])[., r]
            est_r = est[r, id_gt]
            se_r = se[r, id_gt]

            // Skip if SE is missing or zero
            if (se_r == . | se_r <= 0) continue

            // For Epa kernel, check if enough support points exist
            diff = (*precomp.design_mats[id_gt, r])[., 2]

            // Precompute kernel-weighted vectors (n x 1)
            KA     = K_vec :* A_r
            Kdiff  = K_vec :* diff
            Kdiff2 = Kdiff :* diff
            KdiffA = Kdiff :* A_r

            // Matrix products: (n x 1)' * (n x B) -> (1 x B)
            s0 = K_vec' * V_mat
            s1 = Kdiff' * V_mat
            s2 = Kdiff2' * V_mat
            m0 = KA' * V_mat
            m1 = KdiffA' * V_mat

            // Determinant: s0*s2 - s1^2 (1 x B)
            det = s0 :* s2 - s1 :* s1

            // Analytic intercept: beta0 = (s2*m0 - s1*m1) / det
            beta0_row = (s2 :* m0 - s1 :* m1) :/ det

            // Handle singular cases: set to missing where |det| is tiny
            beta0_row = beta0_row :* (abs(det) :> 1e-300)

            // Studentized t-statistic: |beta0 - est_r| / se_r
            t_row = abs(beta0_row :- est_r) / se_r

            // Store in t_mat (row r)
            t_mat[r, .] = t_row
        }

        // Supremum across evaluation points: colmax(t_mat) -> 1 x B
        mb_sup_t[., id_gt] = colmax(t_mat)'
    }
}


// -----------------------------------------------------------------------------
// _dh_boot_vectorized_general()
//
// Partially vectorized bootstrap for general porder (p >= 2).
// Weight matrix is pre-generated; inner WLS per evaluation point uses
// the cross() approach but avoids redundant weight regeneration.
//
// Arguments: same as _dh_boot_vectorized_p1()
// Returns: void
// -----------------------------------------------------------------------------
void _dh_boot_vectorized_general(
    real matrix V_mat,
    struct BootPrecomp scalar precomp,
    pointer(real matrix) rowvector A_g_t,
    real matrix est,
    real matrix se,
    real scalar biters,
    real scalar num_gteval,
    real scalar num_zeval,
    real matrix mb_sup_t)
{
    real scalar id_gt, r, b, n
    real scalar mb_est_r, t_r, sup_t
    real colvector K_vec, A_r, w_eff, RtWA, beta
    real matrix R_mat, RtKR, RtKR_inv
    real rowvector eff_index

    n = precomp.n

    for (b = 1; b <= biters; b++) {

        for (id_gt = 1; id_gt <= num_gteval; id_gt++) {

            sup_t = 0

            for (r = 1; r <= num_zeval; r++) {

                // Effective weights: Mammen_b * kernel
                K_vec = *precomp.kernel_vals[id_gt, r]
                w_eff = V_mat[., b] :* K_vec
                A_r = (*A_g_t[id_gt])[., r]

                eff_index = didhetero_selectindex(K_vec :!= 0)

                if (length(eff_index) == 0) {
                    mb_est_r = .
                }
                else if (length(eff_index) == 1) {
                    mb_est_r = A_r[eff_index]
                }
                else {
                    R_mat = *precomp.design_mats[id_gt, r]
                    RtKR = cross(R_mat, w_eff, R_mat)
                    RtKR_inv = cholinv(RtKR)

                    if (hasmissing(RtKR_inv)) {
                        mb_est_r = .
                    }
                    else {
                        RtWA = cross(R_mat, w_eff, A_r)
                        beta = RtKR_inv * RtWA
                        mb_est_r = beta[1]
                    }
                }

                // Studentized t-statistic
                if (mb_est_r != . & se[r, id_gt] != . & se[r, id_gt] > 0) {
                    t_r = abs(mb_est_r - est[r, id_gt]) / se[r, id_gt]
                    if (t_r > sup_t) sup_t = t_r
                }
            }

            mb_sup_t[b, id_gt] = sup_t
        }
    }
}


// =============================================================================
// Section 6: Main Public Interface
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_bootstrap_ucb()
//
// Unified public entry point for uniform confidence bands via multiplier
// bootstrap. Uses optimized batch processing with pre-computed invariants.
// Weight generation delegates to engine via _didhetero_mammen_weights_batch().
// Critical value extraction delegates to engine via _didhetero_boot_extract_crit().
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
    real colvector c_check,
    | real scalar batch_size,
    real scalar show_progress)
{
    real scalar effective_uniformall
    real matrix mb_sup_t, V_mat
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

    // Generate full n x B Mammen weight matrix via engine
    if (show_progress) {
        printf("Bootstrap: generating %g x %g weight matrix...\n", n, biters)
    }
    V_mat = dh_boot_generate_weights(n, biters, "mammen")

    // Dispatch to vectorized implementation based on polynomial order
    if (porder == 1) {
        // Fully vectorized: analytic 2x2 WLS solution, no b-loop
        if (show_progress) {
            printf("Bootstrap: using vectorized p=1 path (%g iters)...\n", biters)
        }
        _dh_boot_vectorized_p1(V_mat, precomp, A_g_t, est, se, ///
            biters, num_gteval, num_zeval, mb_sup_t)
    }
    else {
        // General porder: loop over b with precomputed weights
        if (show_progress) {
            printf("Bootstrap: using general p=%g path (%g iters)...\n", ///
                porder, biters)
        }
        _dh_boot_vectorized_general(V_mat, precomp, A_g_t, est, se, ///
            biters, num_gteval, num_zeval, mb_sup_t)
    }

    if (show_progress) {
        printf("Bootstrap: %g/%g (100.0%s) done.\n", biters, biters, "%")
    }

    // Extract critical values via engine and build confidence intervals
    _didhetero_boot_extract_crit(mb_sup_t, biters, num_gteval, ///
        effective_uniformall, alp, c_check)
    _didhetero_bootstrap_build_ci(est, se, c_check, num_zeval, num_gteval, ///
        ci2_lower, ci2_upper)
}


// -----------------------------------------------------------------------------
// _didhetero_bootstrap_ucb_basic()
// Legacy basic implementation (retained for debugging/comparison).
// Uses per-iteration weight generation and direct LPR calls without
// pre-computed invariants.
// -----------------------------------------------------------------------------
void _didhetero_bootstrap_ucb_basic(
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
    real scalar b, id_gt, r, col_idx, effective_uniformall
    real matrix mb_est, mb_t, mb_sup_t
    real colvector mb_weight, est_temp, se_temp, t_vals

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

    // Extract critical values and build confidence intervals
    _didhetero_boot_extract_crit(mb_sup_t, biters, num_gteval, ///
        effective_uniformall, alp, c_check)
    _didhetero_bootstrap_build_ci(est, se, c_check, num_zeval, num_gteval, ///
        ci2_lower, ci2_upper)
}


// -----------------------------------------------------------------------------
// didhetero_boot_ucb_optimized()
// Backward-compatibility alias. Forwards to didhetero_bootstrap_ucb().
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
    didhetero_bootstrap_ucb(A_g_t, est, se, bw, Z, zeval, n, porder, ///
        kernel, alp, biters, uniformall, num_gteval, num_zeval, ///
        ci2_lower, ci2_upper, c_check, batch_size, show_progress)
}

end
