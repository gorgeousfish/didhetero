mata:

// =============================================================================
// didhetero_kde.mata
// Kernel density estimation (KDE) with given bandwidth
//
// Implements the standard KDE estimator:
//   f_hat(z) = (1/nh) * sum(K((Z_i - z)/h))
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//     - Section 4.2.2 (density estimation role in bandwidth selection)
//   Section 4.2.2 (kernel density estimation for bandwidth selection)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_kde()
// Kernel density estimator at given evaluation points.
//
// Computes f_hat(z) = (1/nh) * sum_{i=1}^{n} K((x_i - z) / h)
// for each evaluation point z.
//
// Args:
//   x      - n x 1 sample data (covariate Z observations)
//   eval   - num_eval x 1 evaluation points
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//   h      - scalar or num_eval x 1 bandwidth vector
//
// Returns:
//   num_eval x 1 colvector of density estimates
//
// PaperRef: Section 4.2 (kernel density estimation)
// -----------------------------------------------------------------------------
real colvector didhetero_kde(
    real colvector x,
    real colvector eval,
    string scalar kernel,
    real vector h)
{
    real scalar n, num_eval, j
    real colvector Estimate, u, K_vals, h_local


    // --- Dimensions ---
    n = rows(x)
    num_eval = rows(eval)

    // --- Input validation: h vector length ---
    // Check before scalar expansion: must be scalar or match eval length
    if (length(h) != 1 & length(h) != num_eval) {
        _error("h must be scalar or length(eval) vector")
    }

    // --- Bandwidth scalar expansion ---
    // Use local copy to avoid modifying caller's h (Mata passes by reference)
    if (length(h) == 1) {
        h_local = J(num_eval, 1, h)
    }
    else {
        h_local = h
    }

    // --- Input validation: h > 0 ---
    if (any(h_local :<= 0)) {
        _error("h must be positive")
    }

    // --- Initialize output ---
    Estimate = J(num_eval, 1, .)

    // --- Main loop over evaluation points ---
    for (j = 1; j <= num_eval; j++) {

        // Step 1: Standardized distances
        u = (x :- eval[j]) / h_local[j]

        // Step 2: Kernel function values
        K_vals = didhetero_kernel_eval(u, kernel)

        // Step 3: Sum and normalize
        // KDE: f_hat(z) = (1/nh) * sum(K((Z_i - z)/h))
        // For Epanechnikov: sum = 0 when all |u_i| > 1
        // For Gaussian: sum > 0 always
        // Non-negativity: K(u) >= 0, h > 0 => f_hat >= 0
        Estimate[j] = sum(K_vals) / (n * h_local[j])
    }

    return(Estimate)
}

// -----------------------------------------------------------------------------
// didhetero_kde_density()
// Wrapper for kernel density estimation at evaluation points.
// Hardcoded: kernel = "epa", bwselect = "mse-dpi"
//
// Uses didhetero_kdrobust() which implements local polynomial density
// estimation with MSE-DPI automatic bandwidth selection.
// Kernel: Epanechnikov, bandwidth selection: MSE-DPI.
//
// Args:
//   Z     - n x 1 sample data
//   zeval -  x 1 evaluation points
//
// Returns:
//   x 1 colvector of density estimates f_hat(z_r)
//
// PaperRef: Section 4.2.1, density estimation via kernel smoothing
// -----------------------------------------------------------------------------
real colvector didhetero_kde_density(real colvector Z, real colvector zeval)
{
    return(didhetero_kdrobust(Z, zeval, "epa"))
}

// -----------------------------------------------------------------------------
// didhetero_kde_deriv()
// Wrapper for density derivative estimation at evaluation points.
// Hardcoded: kernel = "epa", p = 3, v = 2, bwselect = "mse-dpi"
//
// Uses didhetero_lpdensity() which implements local polynomial density
// derivative estimation with MSE-DPI bandwidth.
// Hardcoded: kernel=Epanechnikov, p=3, v=2, bandwidth selection: MSE-DPI.
//
// Args:
//   Z     - n x 1 sample data
//   zeval -  x 1 evaluation points
//
// Returns:
//   x 1 colvector of density derivative estimates f'_hat(z_r)
//
// PaperRef: Section 4.2.1, density derivative via local polynomial (p=3, v=2)
// -----------------------------------------------------------------------------
real colvector didhetero_kde_deriv(real colvector Z, real colvector zeval)
{
    return(didhetero_lpdensity(Z, zeval, 3, 2, "epa"))
}

end
