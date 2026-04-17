mata:

// =============================================================================
// didhetero_se.mata
// Standard error estimation and analytical uniform confidence band (UCB)
//
// Implements the 4-step variance estimation pipeline and analytical UCB
// as described in Imai, Qin, and Yanagi (2025), Section 4.2.
//
// For each evaluation point z_r:
//   Step 1: IMSE-DPI bandwidth for residual mean of B_{i,g,t}(z_r)
//   Step 2: LPR estimate of E[B_i | Z_i] at all observations
//   Step 3: Residuals → conditional variance σ²_B(z_r)
//   Step 4: Assemble V_hat(z_r) = const_V * σ²_B(z_r) / f_Z(z_r)
//
// Then: SE(z_r) = sqrt(V_hat(z_r) / (n * h_{g,t}))
// Analytical UCB: θ(z_r) ± c_hat * SE(z_r)
//
// Functions:
//   1. didhetero_analytical_crit()        - Analytical critical value (CCK 2014)
//   2. didhetero_se_analytical_ucb()      - Main SE + UCB function
//
// Dependencies:
//   _didhetero_lpbwselect_imse()  - IMSE-DPI bandwidth (in didhetero_bwselect.mata)
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Sections 4.2.3 (SE estimation) and 4.2.4 (analytical UCB)
//   CCK:   Chernozhukov, Chetverikov, Kato (2014) — extreme value UCB
// =============================================================================


// -----------------------------------------------------------------------------
// didhetero_analytical_crit()
// Compute analytical critical value for uniform confidence band.
// Based on extreme value distribution approximation (CCK 2014).
//
// Formula (Paper Eq. 22):
//   a^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
//   c_hat = sqrt(a^2 - 2*log(log(1/sqrt(1-alpha))))
//
// Args:
//   zeval   - num_zeval x 1 evaluation points
//   bw_gt   - scalar, common bandwidth for current (g,t)
//   lambda  - scalar, kernel constant (epa: 2.5, gau: 0.5)
//   alp     - scalar, significance level (e.g. 0.05)
//
// Returns:
//   scalar critical value (. if computation invalid)
//
// Paper ref: Section 4.2.4, analytical critical value
// -----------------------------------------------------------------------------
real scalar didhetero_analytical_crit(
    real colvector zeval,
    real scalar bw_gt,
    real scalar lambda,
    real scalar alp)
{
    real scalar a_val, b_val, a_sq, inner, log_inner, arg, c_hat

    a_val = min(zeval)   // z_1
    b_val = max(zeval)   // z_R

    // Guard: need b > a and bw > 0
    if (b_val <= a_val | bw_gt <= 0) return(.)

    // a^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
    // Note: a_sq CAN be negative — this is fine because the second term
    // 2*log(log(1/sqrt(1-alp))) is typically more negative, so the
    // subtraction a_sq - 2*log(log(...)) yields a positive result.
    a_sq = 2 * log((b_val - a_val) / bw_gt) + 2 * log(sqrt(lambda) / (2 * c("pi")))

    // inner = log(1/sqrt(1-alp))
    if (alp <= 0 | alp >= 1) return(.)
    inner = log(1 / sqrt(1 - alp))

    if (inner <= 0) {
        printf("Warning: log(1/sqrt(1-alp)) <= 0, invalid alp\n")
        return(.)
    }

    // log_inner = log(log(1/sqrt(1-alp)))
    log_inner = log(inner)

    // c_hat = sqrt(a^2 - 2*log(log(1/sqrt(1-alp))))
    arg = a_sq - 2 * log_inner

    if (arg <= 0) {
        printf("Warning: a^2 - 2*log(log(1/sqrt(1-alp))) <= 0\n")
        return(.)
    }

    c_hat = sqrt(arg)

    return(c_hat)
}



// -----------------------------------------------------------------------------
// didhetero_se_analytical_ucb()
// Compute standard errors and analytical UCB for a single (g,t) pair.
// Implements the 4-step variance estimation pipeline.
//
// For each evaluation point z_r:
//   Step 1: IMSE-DPI BW for residual mean (p=porder, eval=Z_supp)
//   Step 2: LPR of B_r on Z at ALL Z (eval=Z, p=porder, h=Step1_BW)
//   Step 3a: Residuals U = B_r - mu_B_hat
//   Step 3b: MSE-DPI BW for U^2 (p=1, eval=zeval[r])
//   Step 3c: LPR of U^2 at zeval[r] (p=1, h=Step3b_BW)
//   Step 4: V_hat[r] = const_V * sigma2 / kd0_Z[r]
//   SE[r] = sqrt(V_hat[r] / (n * bw_gt))
//
// Then analytical UCB:
//   c_hat = analytical critical value
//   CI1_lower = est - c_hat * SE
//   CI1_upper = est + c_hat * SE
//
// Args:
//   B_gt      - n x num_zeval influence function matrix
//   Z         - n x 1 covariate vector
//   zeval     - num_zeval x 1 evaluation points
//   Z_supp    - 100 x 1 support grid for IMSE-DPI BW selection
//   kd0_Z     - num_zeval x 1 kernel density estimates f_hat(z_r)
//   bw_gt     - scalar, common bandwidth for current (g,t)
//   n         - scalar, number of individuals
//   porder    - scalar, polynomial order (1 or 2)
//   kernel    - string, kernel type ("epa" or "gau")
//   const_V   - scalar, variance constant selected by porder
//   lambda    - scalar, kernel constant for analytical critical value
//   alp       - scalar, significance level
//   est_gt    - num_zeval x 1 DR point estimates for current (g,t)
//
// Returns (via pointer arguments):
//   se_gt       - num_zeval x 1 standard errors
//   ci1_l_gt    - num_zeval x 1 analytical UCB lower bounds
//   ci1_u_gt    - num_zeval x 1 analytical UCB upper bounds
//   c_hat_gt    - scalar, analytical critical value for current (g,t)
//   mathcal_V_gt - num_zeval x 1 variance estimates V_hat(z_r)
//
// Paper ref: Section 4.2, variance estimation and analytical UCB
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

    // === Initialize output vectors ===
    se_gt = J(R, 1, .)
    mathcal_V_gt = J(R, 1, .)
    ci1_l_gt = J(R, 1, .)
    ci1_u_gt = J(R, 1, .)

    // === Guard: invalid inputs ===
    if (n <= 0 | bw_gt <= 0) {
        c_hat_gt = .
        return
    }

    // =================================================================
    // Main loop over evaluation points
    // Paper ref: Section 4.2.3, per-evaluation-point variance estimation
    // =================================================================
    for (r = 1; r <= R; r++) {

        // Extract B_{i,g,t}(z_r) column
        B_r = B_gt[., r]

        // =============================================================
        // Step 1: IMSE-DPI bandwidth for residual mean
        // Paper ref: Section 4.2.3, bandwidth for conditional mean of B
        // =============================================================
        mu_B_bw = _didhetero_lpbwselect_imse(B_r, Z, Z_supp,
                      porder, 0, kernel)

        if (mu_B_bw == . | mu_B_bw <= 0) {
            se_gt[r] = .
            continue
        }

        // =============================================================
        // Step 2: LPR estimate of residual mean at ALL observations
        // Paper ref: Section 4.2.3, E[B_i | Z_i] estimation
        // =============================================================
        mu_B_hat = didhetero_lpr(B_r, Z, Z, porder, 0, kernel, mu_B_bw)

        // =============================================================
        // Step 3a: Compute residuals
        // Paper ref: Section 4.2.3, U_hat = B - E[B|Z]
        // =============================================================
        U_hat = B_r - mu_B_hat

        // =============================================================
        // Step 3b: MSE-DPI bandwidth for sigma^2
        // Paper ref: Section 4.2.3, bandwidth for conditional variance
        // =============================================================
        sigma2_bw_vec = _didhetero_lpbwselect_mse(U_hat :^ 2, Z,
                            zeval[r], 1, 0, kernel)
        sigma2_bw_scalar = sigma2_bw_vec[1]

        if (sigma2_bw_scalar == . | sigma2_bw_scalar <= 0) {
            se_gt[r] = .
            continue
        }

        // =============================================================
        // Step 3c: LPR estimate of sigma^2
        // Paper ref: Section 4.2.3, sigma^2_B(z) estimation
        // =============================================================
        real scalar sigma2_raw
        sigma2_raw = didhetero_lpr(U_hat :^ 2, Z, zeval[r],
                         1, 0, kernel, sigma2_bw_scalar)

        // Truncate negative variance to 0.
        // A negative LPR estimate of U^2 at boundary z-points is a known
        // finite-sample artifact (p=1 linear extrapolation near data edge).
        // Flag these cases so SE is set to missing rather than 0.
        sigma2 = sigma2_raw
        if (sigma2 < 0) sigma2 = 0

        // =============================================================
        // Step 4: Variance term V_hat(z_r) = const_V * sigma^2 / f_Z(z_r)
        // Paper ref: Section 4.2.3, Eq. 19 variance assembly
        // =============================================================
        if (kd0_Z[r] > 0) {
            V_hat_r = const_V * sigma2 / kd0_Z[r]
        }
        else {
            V_hat_r = .
        }

        // Store V_hat(z_r) into output vector
        mathcal_V_gt[r] = V_hat_r

        // =============================================================
        // SE computation: SE(z_r) = sqrt(V_hat(z_r) / (n * h))
        // Paper ref: Section 4.2.3
        // =============================================================
        if (sigma2_raw < 0) {
            // sigma2_raw < 0: LPR of U^2 extrapolates below zero at this z_r.
            // This reliably indicates insufficient local data near z_r
            // (typically a boundary evaluation point far from treated units).
            // Setting SE=. is correct: a nonparametric estimator cannot have
            // zero variance, so SE=0 would produce a degenerate [est,est] CI.
            se_gt[r] = .
        }
        else if (V_hat_r != . & V_hat_r >= 0) {
            se_gt[r] = sqrt(V_hat_r / (n * bw_gt))
            // Scale-invariant SE guard: flag SE when it is implausibly large
            // relative to the outcome scale. When zeval[r] falls far outside
            // the treated cohort's support, GPS≈0 inflates B_r—and hence
            // Var(B_r)—so a Var(B_r)-based upper bound also scales up and
            // cannot catch the pathological case.
            //
            // Instead, compare SE to a scale reference derived from the
            // point estimates themselves: est_scale = median_abs(est_gt) + 1.
            // Guard: flag as missing when SE > 1e3 * est_scale.
            // 1e3 allows SEs up to 1000× the outcome scale while catching
            // blow-ups that are 10^20+ times larger.
            est_scale = median(abs(est_gt)) + 1
            if (se_gt[r] > 1e3 * est_scale) {
                se_gt[r] = .
            }
        }
        else {
            se_gt[r] = .
        }
    }

    // =================================================================
    // Analytical critical value via extreme value approximation
    // Paper ref: Section 4.2.4, Eq. 22 (CCK 2014)
    // =================================================================
    c_hat_gt = didhetero_analytical_crit(zeval, bw_gt, lambda, alp)

    // =================================================================
    // Analytical UCB invalid-state handling
    // When the analytical joint critical value cannot be computed, keep the
    // analytical band unavailable instead of silently switching to a
    // pointwise normal interval. This matches the paper's specification.
    // =================================================================
    if (c_hat_gt == .) {
        printf("{txt}Warning: analytical UCB critical value cannot be computed; leaving analytical CI missing\n")
    }

    // =================================================================
    // Analytical UCB construction: est +/- c_hat * SE
    // Paper ref: Section 2, Eq. 3 (UCB structure)
    // =================================================================
    ci1_l_gt = est_gt - c_hat_gt * se_gt
    ci1_u_gt = est_gt + c_hat_gt * se_gt

    // === Diagnostic: count missing SEs ===
    n_missing = sum(se_gt :== .)
    if (n_missing > 0) {
        printf("Warning: SE=. at %g of %g evaluation point(s): insufficient local data\n",
               n_missing, R)
        printf("  (boundary z-values far from treated units; CI reported as missing)\n")
    }
}

end
