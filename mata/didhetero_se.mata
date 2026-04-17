mata:

// =============================================================================
// Standard error estimation and analytical uniform confidence bands
//
// Implements the 4-step variance estimation pipeline for each evaluation
// point z_r:
//   Step 1: Bandwidth selection for conditional mean of influence function
//   Step 2: Local polynomial regression of influence function on covariate
//   Step 3: Residuals used to estimate conditional variance
//   Step 4: Variance assembly and standard error computation
//
// Analytical uniform confidence bands computed via extreme value approximation.
// =============================================================================


// -----------------------------------------------------------------------------
// didhetero_analytical_crit()
// Compute analytical critical value for uniform confidence bands.
// Based on extreme value distribution approximation.
//
// Formula:
//   a^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
//   c_hat = sqrt(a^2 - 2*log(log(1/sqrt(1-alpha))))
//
// Arguments:
//   zeval   - num_zeval x 1 evaluation points
//   bw_gt   - scalar, common bandwidth for current (g,t)
//   lambda  - scalar, kernel constant (epa: 2.5, gau: 0.5)
//   alp     - scalar, significance level (e.g. 0.05)
//
// Returns:
//   scalar critical value (. if computation invalid)
// -----------------------------------------------------------------------------
real scalar didhetero_analytical_crit(
    real colvector zeval,
    real scalar bw_gt,
    real scalar lambda,
    real scalar alp)
{
    real scalar a_val, b_val, a_sq, inner, log_inner, arg, c_hat

    a_val = min(zeval)
    b_val = max(zeval)

    // Guard: need valid range and bandwidth
    if (b_val <= a_val | bw_gt <= 0) return(.)

    // Compute a^2 term
    a_sq = 2 * log((b_val - a_val) / bw_gt) + 2 * log(sqrt(lambda) / (2 * c("pi")))

    // Guard: valid significance level
    if (alp <= 0 | alp >= 1) return(.)
    inner = log(1 / sqrt(1 - alp))

    if (inner <= 0) {
        printf("Warning: log(1/sqrt(1-alp)) <= 0, invalid alp\n")
        return(.)
    }

    log_inner = log(inner)
    arg = a_sq - 2 * log_inner

    if (arg <= 0) {
        printf("Warning: argument to square root non-positive\n")
        return(.)
    }

    c_hat = sqrt(arg)

    return(c_hat)
}



// -----------------------------------------------------------------------------
// didhetero_se_analytical_ucb()
// Compute standard errors and analytical uniform confidence bands for a
// single (g,t) pair. Implements the 4-step variance estimation pipeline.
//
// For each evaluation point z_r:
//   Step 1: Bandwidth selection for conditional mean of influence function
//   Step 2: Local polynomial regression of influence function on covariate
//   Step 3a: Compute residuals
//   Step 3b: Bandwidth selection for conditional variance
//   Step 3c: Local polynomial regression of squared residuals
//   Step 4: Assemble variance and compute standard error
//
// Uniform confidence bands constructed using extreme value critical value.
//
// Arguments:
//   B_gt      - n x num_zeval influence function matrix
//   Z         - n x 1 covariate vector
//   zeval     - num_zeval x 1 evaluation points
//   Z_supp    - 100 x 1 support grid for bandwidth selection
//   kd0_Z     - num_zeval x 1 kernel density estimates
//   bw_gt     - scalar, common bandwidth for current (g,t)
//   n         - scalar, number of individuals
//   porder    - scalar, polynomial order (1 or 2)
//   kernel    - string, kernel type ("epa" or "gau")
//   const_V   - scalar, variance constant
//   lambda    - scalar, kernel constant for critical value
//   alp       - scalar, significance level
//   est_gt    - num_zeval x 1 point estimates for current (g,t)
//
// Returns (via pointer arguments):
//   se_gt       - num_zeval x 1 standard errors
//   ci1_l_gt    - num_zeval x 1 lower confidence bounds
//   ci1_u_gt    - num_zeval x 1 upper confidence bounds
//   c_hat_gt    - scalar, critical value
//   mathcal_V_gt - num_zeval x 1 variance estimates
// -----------------------------------------------------------------------------
void didhetero_se_analytical_ucb(
    real matrix B_gt,
    real colvector Z,
    real colvector zeval,
    real colvector Z_supp,
    real colvector kd0_Z,
    real scalar bw_gt,
    real scalar n,
    real scalar porder,
    string scalar kernel,
    real scalar const_V,
    real scalar lambda,
    real scalar alp,
    real colvector est_gt,
    real colvector se_gt,
    real colvector ci1_l_gt,
    real colvector ci1_u_gt,
    real scalar c_hat_gt,
    real colvector mathcal_V_gt)
{
    real scalar R, r, mu_B_bw, sigma2_bw_scalar, sigma2, V_hat_r
    real scalar n_missing, est_scale
    real colvector B_r, mu_B_hat, U_hat, sigma2_bw_vec

    R = rows(zeval)

    // Initialize output vectors
    se_gt = J(R, 1, .)
    mathcal_V_gt = J(R, 1, .)
    ci1_l_gt = J(R, 1, .)
    ci1_u_gt = J(R, 1, .)

    // Guard against invalid inputs
    if (n <= 0 | bw_gt <= 0) {
        c_hat_gt = .
        return
    }

    // Main loop over evaluation points
    for (r = 1; r <= R; r++) {

        // Extract influence function column for current evaluation point
        B_r = B_gt[., r]

        // Bandwidth selection for conditional mean of influence function
        mu_B_bw = _didhetero_lpbwselect_imse(B_r, Z, Z_supp,
                      porder, 0, kernel)

        if (mu_B_bw == . | mu_B_bw <= 0) {
            se_gt[r] = .
            continue
        }

        // Local polynomial regression of influence function on covariate
        mu_B_hat = didhetero_lpr(B_r, Z, Z, porder, 0, kernel, mu_B_bw)

        // Compute residuals
        U_hat = B_r - mu_B_hat

        // Bandwidth selection for conditional variance
        sigma2_bw_vec = _didhetero_lpbwselect_mse(U_hat :^ 2, Z,
                            zeval[r], 1, 0, kernel)
        sigma2_bw_scalar = sigma2_bw_vec[1]

        if (sigma2_bw_scalar == . | sigma2_bw_scalar <= 0) {
            se_gt[r] = .
            continue
        }

        // Local polynomial regression of squared residuals
        real scalar sigma2_raw
        sigma2_raw = didhetero_lpr(U_hat :^ 2, Z, zeval[r],
                         1, 0, kernel, sigma2_bw_scalar)

        // Truncate negative variance estimates to zero
        sigma2 = sigma2_raw
        if (sigma2 < 0) sigma2 = 0

        // Assemble variance estimate
        if (kd0_Z[r] > 0) {
            V_hat_r = const_V * sigma2 / kd0_Z[r]
        }
        else {
            V_hat_r = .
        }

        // Store variance estimate
        mathcal_V_gt[r] = V_hat_r

        // Compute standard error
        if (sigma2_raw < 0) {
            // Negative variance indicates insufficient local data
            se_gt[r] = .
        }
        else if (V_hat_r != . & V_hat_r >= 0) {
            se_gt[r] = sqrt(V_hat_r / (n * bw_gt))
            // Scale-invariant guard against pathologically large SEs
            est_scale = median(abs(est_gt)) + 1
            if (se_gt[r] > 1e3 * est_scale) {
                se_gt[r] = .
            }
        }
        else {
            se_gt[r] = .
        }
    }

    // Compute analytical critical value via extreme value approximation
    c_hat_gt = didhetero_analytical_crit(zeval, bw_gt, lambda, alp)

    // Handle invalid critical value
    if (c_hat_gt == .) {
        printf("{txt}Warning: analytical UCB critical value cannot be computed; leaving analytical CI missing\n")
    }

    // Construct uniform confidence bands
    ci1_l_gt = est_gt - c_hat_gt * se_gt
    ci1_u_gt = est_gt + c_hat_gt * se_gt

    // Report missing standard errors
    n_missing = sum(se_gt :== .)
    if (n_missing > 0) {
        printf("Warning: SE=. at %g of %g evaluation point(s): insufficient local data\n",
               n_missing, R)
        printf("  (boundary z-values far from treated units; CI reported as missing)\n")
    }
}

end
