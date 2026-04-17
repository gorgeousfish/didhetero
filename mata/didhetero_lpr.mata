mata:

// =============================================================================
// Local polynomial regression (LPR) estimator
//
// Implements p-th order local polynomial regression for:
//   - Conditional mean estimation (p=1 local linear, p=2 local quadratic)
//   - Derivative estimation (deriv=0,1,2,...)
//   - Weighted estimation with user-specified observation weights
//
// Two functions:
//   1. didhetero_polynomial() - Polynomial design matrix construction
//   2. didhetero_lpr()        - Main LPR estimator
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_polynomial()
// Construct polynomial design matrix [1, u, u^2, ..., u^p]
//
// Parameters:
//   u - n x 1 colvector of centered values (x - eval_point)
//   p - polynomial order (non-negative integer)
//
// Returns:
//   n x (p+1) matrix where column j+1 contains u^j for j=0,...,p
// -----------------------------------------------------------------------------
real matrix didhetero_polynomial(real colvector u, real scalar p)
{
    real matrix result
    real scalar j, n
    
    n = rows(u)
    result = J(n, p + 1, .)
    
    // Construct column j+1 as u raised to power j
    for (j = 0; j <= p; j++) {
        result[., j + 1] = u :^ j
    }
    
    return(result)
}


// -----------------------------------------------------------------------------
// didhetero_lpr()
// Local polynomial regression estimator
//
// Solves weighted least squares at each evaluation point:
//   beta(z) = (Z' K_h Z)^{-1} Z' K_h y
// Returns the nu-th derivative estimate: nu! * e_nu' * beta(z)
//
// Parameters:
//   y      - n x 1 dependent variable
//   x      - n x 1 independent variable
//   eval   - num_eval x 1 evaluation points
//   p      - polynomial order (0=Nadaraya-Watson, 1=local linear, 2=local quadratic)
//   deriv  - derivative order (0 <= deriv <= p)
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//   h      - scalar or num_eval x 1 bandwidth vector
//   weight - n x 1 observation weights (optional, defaults to unit weights)
//
// Returns:
//   num_eval x 1 colvector of LPR estimates
// -----------------------------------------------------------------------------
real colvector didhetero_lpr(
    real colvector y,
    real colvector x,
    real colvector eval,
    real scalar p,
    real scalar deriv,
    string scalar kernel,
    real vector h,
    | real colvector weight)
{
    real scalar n, num_eval, j
    real colvector Estimate, K_h, y_sub, K_h_sub, e_vec, beta, h_local
    real matrix Z_mat, Z_sub, Gamma, Gamma_inv
    real rowvector index
    
    // Dimensions
    n = rows(y)
    num_eval = rows(eval)
    
    // Input validation
    if (rows(x) != n) {
        _error("y and x must have the same length")
    }
    if (p < 0) {
        _error("p must be non-negative")
    }
    if (deriv > p) {
        _error("deriv must be <= p")
    }
    
    // Default weights if not provided
    if (args() < 8) weight = J(n, 1, 1)
    
    // Bandwidth handling: replicate scalar bandwidth for all evaluation points
    if (length(h) == 1) {
        h_local = J(num_eval, 1, h)
    }
    else {
        h_local = h
    }
    if (length(h_local) != num_eval) {
        _error("h must be scalar or length(eval) vector")
    }
    
    // Unit vector to extract the (deriv+1)-th coefficient from beta
    e_vec = J(p + 1, 1, 0)
    e_vec[deriv + 1] = 1
    
    // Initialize output vector
    Estimate = J(num_eval, 1, .)
    
    // Main loop over evaluation points
    for (j = 1; j <= num_eval; j++) {
        
        // Construct polynomial design matrix Z_mat at evaluation point eval[j]
        Z_mat = didhetero_polynomial(x :- eval[j], p)
        
        // Compute kernel weights multiplied by user-specified weights
        K_h = weight :* didhetero_kernel_eval((x :- eval[j]) / h_local[j], kernel)
        
        // Select observations with positive kernel weights
        index = didhetero_selectindex(K_h :!= 0)
        
        // Compute estimate based on number of effective observations
        if (length(index) == 0) {
            // No effective observations: return missing value
            Estimate[j] = .
        }
        else if (length(index) == 1) {
            // Single observation: return observed value
            Estimate[j] = y[index]
        }
        else {
            // Multiple observations: solve weighted least squares
            
            // Extract subsets with positive weights
            Z_sub   = Z_mat[index', .]
            y_sub   = y[index']
            K_h_sub = K_h[index']
            
            // Compute weighted Gram matrix
            Gamma = cross(Z_sub, K_h_sub, Z_sub)
            
            // Compute inverse via Cholesky decomposition
            Gamma_inv = cholinv(Gamma)
            
            // Handle singular or ill-conditioned Gram matrix
            if (hasmissing(Gamma_inv)) {
                Estimate[j] = .
            }
            else {
                // Compute weighted least squares coefficients
                beta = Gamma_inv * cross(Z_sub, K_h_sub, y_sub)
                
                // Compute derivative estimate by extracting relevant coefficient
                Estimate[j] = factorial(deriv) * (e_vec' * beta)
            }
        }
    }
    
    return(Estimate)
}

end
