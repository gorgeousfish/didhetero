// =============================================================================
// didhetero_boot.mata
// Bootstrap weight generation for didhetero-stata
//
// Implements Mammen (1993) two-point distribution weights.
// Satisfies Assumption 10: E[V*]=1, Var[V*]=1, sub-exponential tails.
//
// Paper ref: Section 4.2.4 (Mammen two-point distribution for multiplier bootstrap)
// =============================================================================

mata:

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights()
//
// Generate Mammen (1993) two-point distribution weights.
// Each element is i.i.d. from the Mammen two-point distribution.
// Satisfies E[V*] = 1, Var[V*] = 1 (Assumption 10).
//
// Paper ref: Assumption 10, Mammen two-point distribution for multiplier bootstrap.
//
// Parameters:
//   n - sample size (positive integer)
//
// Returns:
//   n x 1 vector of Mammen weights
// -----------------------------------------------------------------------------
real colvector _didhetero_mammen_weights(real scalar n)
{
    real scalar kappa, v1, v2, p1
    real colvector u, weights

    // Input validation
    if (missing(n)) {
        _error("n is missing")
    }
    if (n <= 0) {
        _error("n must be positive")
    }
    if (n != floor(n)) {
        errprintf("warning: n is not integer, truncating to %g\n", floor(n))
        n = floor(n)
    }

    // Constants: golden ratio kappa = (sqrt(5)+1)/2
    kappa = (sqrt(5) + 1) / 2
    v1 = 2 - kappa
    v2 = 1 + kappa
    p1 = kappa / sqrt(5)

    // Generate uniform random numbers
    u = runiform(n, 1)

    // Vectorized two-point distribution via inverse CDF
    // Align threshold convention with bootstrap use sites: u < p1 -> v1; u >= p1 -> v2
    weights = (u :< p1) :* v1 + (u :>= p1) :* v2

    return(weights)
}

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights_batch()
//
// Batch generate Mammen weights for all bootstrap iterations.
// Each row is an independent draw of n i.i.d. Mammen weights.
// Used by aggte multiplier bootstrap.
//
// Paper ref: Section 5, multiplier bootstrap for aggregation.
//
// Parameters:
//   B - number of bootstrap iterations (positive integer)
//   n - sample size (positive integer)
//
// Returns:
//   B x n matrix of Mammen weights
// -----------------------------------------------------------------------------
real matrix _didhetero_mammen_weights_batch(real scalar B, real scalar n)
{
    real scalar kappa, v1, v2, p1
    real matrix u, weights_mat

    // Input validation
    if (missing(B)) {
        _error("B is missing")
    }
    if (missing(n)) {
        _error("n is missing")
    }
    if (B <= 0) {
        _error("B must be positive")
    }
    if (n <= 0) {
        _error("n must be positive")
    }
    if (B != floor(B)) {
        errprintf("warning: B is not integer, truncating to %g\n", floor(B))
        B = floor(B)
    }
    if (n != floor(n)) {
        errprintf("warning: n is not integer, truncating to %g\n", floor(n))
        n = floor(n)
    }

    // Constants (same as single generation)
    kappa = (sqrt(5) + 1) / 2
    v1 = 2 - kappa
    v2 = 1 + kappa
    p1 = kappa / sqrt(5)

    // Generate B x n uniform random matrix
    u = runiform(B, n)

    // Vectorized two-point distribution
    weights_mat = (u :< p1) :* v1 + (u :>= p1) :* v2

    return(weights_mat)
}

end
