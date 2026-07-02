mata:

// =============================================================================
// didhetero_utils_numerical.mata
// Numerical computation utility functions
//
// Provides:
//   - didhetero_seq()             // Generate equispaced sequence
//   - didhetero_trapz()           // Trapezoidal rule integration
//   - didhetero_quantile()        // Type-7 quantile function
//   - didhetero_mammen_weights()  // Mammen two-point bootstrap weights
//   - didhetero_analytical_cv()   // Analytical critical value (legacy)
//   - _didhetero_intersect()      // Set intersection for colvectors
//   - _didhetero_setdiff()        // Set difference for colvectors
//
// Requires:
//   - (none -- standalone numerical utilities)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_seq()
// Generate equispaced sequence from `from` to `to` with `length` points
//
// Args:
//   from   - start value
//   to     - end value
//   length - number of points (must be >= 1)
//
// Returns:
//   real colvector of equispaced values
//
// -----------------------------------------------------------------------------
real colvector didhetero_seq(real scalar from, real scalar to, real scalar length)
{
    // Input validation
    if (length < 1) _error("length must be at least 1")
    
    // Special case: single point
    if (length == 1) return(from)
    
    return(rangen(from, to, length))
}

// -----------------------------------------------------------------------------
// didhetero_trapz()
// Trapezoidal rule numerical integration
//
// Args:
//   x - real colvector of grid points (abscissae)
//   y - real colvector of function values at grid points
//
// Returns:
//   real scalar approximation of integral of y over x
//
// -----------------------------------------------------------------------------
real scalar didhetero_trapz(real colvector x, real colvector y)
{
    real scalar n, integral, i
    
    // Input validation: x and y must have same length
    n = length(x)
    if (n != length(y)) _error("x and y must have same length")
    
    // Edge case: fewer than 2 points, no interval to integrate
    if (n < 2) return(0)
    
    // Trapezoidal rule: sum of (x[i+1]-x[i]) * (y[i]+y[i+1]) / 2
    integral = 0
    for (i = 1; i < n; i++) {
        integral = integral + (x[i+1] - x[i]) * (y[i] + y[i+1]) / 2
    }
    return(integral)
}

// -----------------------------------------------------------------------------
// didhetero_quantile()
// Quantile function using type-7 algorithm
//
// Args:
//   x     - real colvector of data values (must be non-empty)
//   probs - real vector of probabilities in [0, 1]
//
// Returns:
//   real vector of quantile values corresponding to probs
//
// Type 7 continuous index: h = (n-1)*p + 1, with linear interpolation
// -----------------------------------------------------------------------------
real vector didhetero_quantile(real colvector x, real vector probs)
{
    real colvector sorted_x
    real vector result
    real scalar n, i, h_idx, lo, hi
    
    // Input validation: x must be non-empty
    if (length(x) < 1) _error("x must be non-empty")
    
    // Input validation: probs must be in [0, 1]
    if (min(probs) < 0 | max(probs) > 1) _error("probs must be between 0 and 1")
    
    // Sort data ascending
    sorted_x = sort(x, 1)
    n = length(sorted_x)
    result = J(length(probs), 1, .)
    
    for (i = 1; i <= length(probs); i++) {
        // Linear interpolation quantile: continuous index h = (n-1)*p + 1
        h_idx = (n - 1) * probs[i] + 1
        lo = floor(h_idx)
        hi = ceil(h_idx)
        
        // Boundary protection
        if (lo < 1) lo = 1
        if (hi > n) hi = n
        
        if (lo == hi) {
            result[i] = sorted_x[lo]
        }
        else {
            // Linear interpolation
            result[i] = sorted_x[lo] + (h_idx - lo) * (sorted_x[hi] - sorted_x[lo])
        }
    }
    return(result)
}

// -----------------------------------------------------------------------------
// didhetero_mammen_weights()
// Generate Mammen two-point wild bootstrap weights
//
// Args:
//   n - number of weights to generate (must be >= 1)
//
// Returns:
//   real colvector of bootstrap weights
//
// kappa = (sqrt(5)+1)/2 (golden ratio)
// V* = 2-kappa w.p. kappa/sqrt(5), 1+kappa w.p. 1-kappa/sqrt(5)
// Moment properties: E[V]=1, Var[V]=1, E[(V-1)^3]=1
// -----------------------------------------------------------------------------
real colvector didhetero_mammen_weights(real scalar n)
{
    real scalar kappa, p_low
    real colvector u, weights
    
    // Input validation
    if (n < 1) _error("n must be a positive integer")
    
    // Golden ratio
    kappa = (sqrt(5) + 1) / 2
    p_low = kappa / sqrt(5)
    
    // Two-point distribution via uniform threshold
    // V* = 2-kappa with prob kappa/sqrt(5), 1+kappa with prob 1-kappa/sqrt(5)
    u = runiform(n, 1)
    weights = (u :<= p_low) :* (2 - kappa) + (u :> p_low) :* (1 + kappa)
    
    return(weights)
}

// -----------------------------------------------------------------------------
// didhetero_analytical_cv()
// Analytical critical value from extreme value distribution
//
// Args:
//   a      - lower bound of the interval
//   b      - upper bound of the interval (must be > a)
//   h      - bandwidth (must be > 0)
//   lambda - kernel-dependent constant
//   alpha  - significance level in (0, 1)
//
// Returns:
//   real scalar critical value c_hat
//
// a_n^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
// c_hat = sqrt(a_n^2 - 2*log(log(1/sqrt(1-alpha))))
// -----------------------------------------------------------------------------
real scalar didhetero_analytical_cv(real scalar a, real scalar b,
                                     real scalar h, real scalar lambda,
                                     real scalar alpha)
{
    real scalar a_n_sq, c_hat, inner_log
    
    // Input validation
    if (h <= 0) _error("h must be positive")
    if (b <= a) _error("b must be greater than a")
    if (alpha <= 0 | alpha >= 1) _error("alpha must be between 0 and 1 (exclusive)")
    
    // a_n^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
    a_n_sq = 2 * log((b - a) / h) + 2 * log(sqrt(lambda) / (2 * c("pi")))
    
    // inner_log = log(log(1/sqrt(1-alpha)))
    inner_log = log(log(1 / sqrt(1 - alpha)))
    
    // Boundary protection: fallback to normal quantile when a_n^2 is too small
    if (a_n_sq - 2 * inner_log <= 0) {
        c_hat = invnormal(1 - alpha/2)
    }
    else {
        // c_hat = sqrt(a_n_sq - 2*inner_log)
        c_hat = sqrt(a_n_sq - 2 * inner_log)
    }
    
    return(c_hat)
}

// -----------------------------------------------------------------------------
// _didhetero_intersect()
// Return sorted intersection of two vectors
//
// Args:
//   a - real colvector
//   b - real colvector
//
// Returns:
//   real colvector of elements present in both a and b, sorted ascending
//   Returns J(0, 1, .) if intersection is empty
// -----------------------------------------------------------------------------
real colvector _didhetero_intersect(real colvector a, real colvector b)
{
    real colvector result
    real scalar i, n_a
    
    // Handle empty inputs
    n_a = rows(a)
    if (n_a == 0 | rows(b) == 0) return(J(0, 1, .))
    
    result = J(0, 1, .)
    for (i = 1; i <= n_a; i++) {
        if (any(b :== a[i])) {
            result = result \ a[i]
        }
    }
    
    if (rows(result) == 0) return(J(0, 1, .))
    return(sort(result, 1))
}

// -----------------------------------------------------------------------------
// _didhetero_setdiff()
// Return sorted set difference: elements in a but not in b
//
// Args:
//   a - real colvector
//   b - real colvector
//
// Returns:
//   real colvector of elements in a but not in b, sorted ascending
//   Returns J(0, 1, .) if result is empty
// -----------------------------------------------------------------------------
real colvector _didhetero_setdiff(real colvector a, real colvector b)
{
    real colvector result
    real scalar i, n_a
    
    // Handle empty inputs
    n_a = rows(a)
    if (n_a == 0) return(J(0, 1, .))
    if (rows(b) == 0) return(sort(a, 1))
    
    result = J(0, 1, .)
    for (i = 1; i <= n_a; i++) {
        if (!any(b :== a[i])) {
            result = result \ a[i]
        }
    }
    
    if (rows(result) == 0) return(J(0, 1, .))
    return(sort(result, 1))
}

end
