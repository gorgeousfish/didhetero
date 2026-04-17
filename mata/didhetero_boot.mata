// =============================================================================
// didhetero_boot.mata
// Bootstrap weight generation
//
// Implements Mammen (1993) two-point distribution for multiplier bootstrap.
// Weights satisfy E[V*] = 1, Var[V*] = 1, and sub-exponential tails.
// =============================================================================

mata:

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights()
//
// Generate Mammen (1993) two-point distribution weights.
// Returns an n x 1 vector of i.i.d. draws satisfying E[V*] = 1, Var[V*] = 1.
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

    // Vectorized two-point distribution: u < p1 yields v1, u >= p1 yields v2
    weights = (u :< p1) :* v1 + (u :>= p1) :* v2

    return(weights)
}

// -----------------------------------------------------------------------------
// _didhetero_mammen_weights_batch()
//
// Generate Mammen weights for all bootstrap iterations.
// Returns a B x n matrix where each row is an independent draw of n weights.
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

    // Constants for two-point distribution
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
