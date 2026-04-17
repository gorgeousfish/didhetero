mata:

// =============================================================================
// Kernel density estimation with given bandwidth
//
// Implements the standard KDE estimator:
//   f_hat(z) = (1/nh) * sum(K((x_i - z)/h))
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_kde()
// Kernel density estimator at given evaluation points.
//
// Computes f_hat(z) = (1/nh) * sum_{i=1}^{n} K((x_i - z) / h)
// for each evaluation point z.
//
// Parameters:
//   x      - n x 1 sample data
//   eval   - num_eval x 1 evaluation points
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//   h      - scalar or num_eval x 1 bandwidth vector
//
// Returns:
//   num_eval x 1 colvector of density estimates
// -----------------------------------------------------------------------------
real colvector didhetero_kde(
    real colvector x,
    real colvector eval,
    string scalar kernel,
    real vector h)
{
    real scalar n, num_eval, j
    real colvector Estimate, u, K_vals, h_local


    // Dimensions
    n = rows(x)
    num_eval = rows(eval)

    // Validate bandwidth vector length
    if (length(h) != 1 & length(h) != num_eval) {
        _error("h must be scalar or length(eval) vector")
    }

    // Bandwidth scalar expansion (use local copy to preserve input)
    if (length(h) == 1) {
        h_local = J(num_eval, 1, h)
    }
    else {
        h_local = h
    }

    // Validate bandwidth is positive
    if (any(h_local :<= 0)) {
        _error("h must be positive")
    }

    // Initialize output vector
    Estimate = J(num_eval, 1, .)

    // Compute density at each evaluation point
    for (j = 1; j <= num_eval; j++) {

        // Standardized distances
        u = (x :- eval[j]) / h_local[j]

        // Kernel function values
        K_vals = didhetero_kernel_eval(u, kernel)

        // Sum and normalize: f_hat(z) = (1/nh) * sum(K((x_i - z)/h))
        Estimate[j] = sum(K_vals) / (n * h_local[j])
    }

    return(Estimate)
}

// -----------------------------------------------------------------------------
// didhetero_kde_density()
// Kernel density estimation at evaluation points.
//
// Uses Epanechnikov kernel with MSE-DPI bandwidth selection via
// local polynomial density estimation.
//
// Parameters:
//   Z     - n x 1 sample data
//   zeval - m x 1 evaluation points
//
// Returns:
//   m x 1 colvector of density estimates
// -----------------------------------------------------------------------------
real colvector didhetero_kde_density(real colvector Z, real colvector zeval)
{
    return(didhetero_kdrobust(Z, zeval, "epa"))
}

// -----------------------------------------------------------------------------
// didhetero_kde_deriv()
// Density derivative estimation at evaluation points.
//
// Computes second-order density derivatives using local polynomial
// regression with Epanechnikov kernel (p=3, v=2) and MSE-DPI bandwidth.
//
// Parameters:
//   Z     - n x 1 sample data
//   zeval - m x 1 evaluation points
//
// Returns:
//   m x 1 colvector of density derivative estimates
// -----------------------------------------------------------------------------
real colvector didhetero_kde_deriv(real colvector Z, real colvector zeval)
{
    return(didhetero_lpdensity(Z, zeval, 3, 2, "epa"))
}

end
