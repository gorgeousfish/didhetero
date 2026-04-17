mata:

// =============================================================================
// didhetero_lpr.mata
// Local polynomial regression (LPR) estimator
//
// Implements p-th order local polynomial regression for:
//   - Conditional mean estimation (p=1 LLR, p=2 LQR)
//   - Derivative estimation (deriv=0,1,2,...)
//   - Bootstrap weighted estimation (Mammen wild bootstrap)
//
// Two functions:
//   1. didhetero_polynomial() - Polynomial design matrix construction
//   2. didhetero_lpr()        - Main LPR estimator
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//     - Section 4.2.1 Eq. 12-13 (LPR estimator definition)
//     - Section 4.2.4 Eq. 23 (Bootstrap weighted LPR)
//   Section 4.2.1 Eq. 12-13 (LPR estimator, polynomial design matrix)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_polynomial()
// Construct polynomial design matrix [1, u, u^2, ..., u^p]
//
// Args:
//   u - n x 1 colvector of centered values (x - eval_point)
//   p - polynomial order (non-negative integer)
//
// Returns:
//   n x (p+1) matrix where column j+1 = u^j, j=0,...,p
//
// Paper ref: Eq. 12, polynomial basis r_p(u) = (1, u, ..., u^p)'
// -----------------------------------------------------------------------------
real matrix didhetero_polynomial(real colvector u, real scalar p)
{
    real matrix result
    real scalar j, n
    
    n = rows(u)
    result = J(n, p + 1, .)
    
    // Column j+1 = u^j (j=0 gives all-ones column)
    for (j = 0; j <= p; j++) {
        result[., j + 1] = u :^ j
    }
    
    return(result)
}


// -----------------------------------------------------------------------------
// didhetero_lpr()
// Local polynomial regression estimator.
//
// Solves weighted least squares at each evaluation point:
//   beta(z) = (Z' K_h Z)^{-1} Z' K_h y
// Returns nu-th derivative estimate: nu! * e_nu' * beta(z)
//
// Args:
//   y      - n x 1 dependent variable
//   x      - n x 1 independent variable (covariate Z)
//   eval   - num_eval x 1 evaluation points
//   p      - polynomial order (0=NW, 1=LLR, 2=LQR, ...)
//   deriv  - derivative order (0 <= deriv <= p)
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//   h      - scalar or num_eval x 1 bandwidth vector
//   weight - n x 1 observation weights (optional, default all 1s)
//
// Returns:
//   num_eval x 1 colvector of LPR estimates
//
// Paper ref: Section 4.2.1 Eq. 12-13, local polynomial estimation
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
    
    // --- Dimensions ---
    n = rows(y)
    num_eval = rows(eval)
    
    // --- Input validation ---
    if (rows(x) != n) {
        _error("y and x must have the same length")
    }
    if (p < 0) {
        _error("p must be non-negative")
    }
    if (deriv > p) {
        _error("deriv must be <= p")
    }
    
    // --- Optional weight default ---
    if (args() < 8) weight = J(n, 1, 1)
    
    // --- Bandwidth handling ---
    // Use local copy to avoid modifying caller's h (Mata passes by reference)
    if (length(h) == 1) {
        h_local = J(num_eval, 1, h)
    }
    else {
        h_local = h
    }
    if (length(h_local) != num_eval) {
        _error("h must be scalar or length(eval) vector")
    }
    
    // --- Unit vector for derivative extraction ---
    // e_vec selects the (deriv+1)-th coefficient from beta
    e_vec = J(p + 1, 1, 0)
    e_vec[deriv + 1] = 1
    
    // --- Initialize output ---
    Estimate = J(num_eval, 1, .)
    
    // --- Main loop over evaluation points ---
    for (j = 1; j <= num_eval; j++) {
        
        // Step 1: Polynomial design matrix
        // Z_mat[i,.] = (1, x_i-z_j, (x_i-z_j)^2, ..., (x_i-z_j)^p)
        Z_mat = didhetero_polynomial(x :- eval[j], p)
        
        // Step 2: Kernel weights * user weights
        // K_h[i] = weight[i] * K((x_i - z_j) / h_j)
        K_h = weight :* didhetero_kernel_eval((x :- eval[j]) / h_local[j], kernel)
        
        // Step 3: Select non-zero weight observations
        index = didhetero_selectindex(K_h :!= 0)
        
        // Step 4: Branch on effective observation count
        if (length(index) == 0) {
            // Case A: No effective observations
            Estimate[j] = .
        }
        else if (length(index) == 1) {
            // Case B: Single observation
            Estimate[j] = y[index]
        }
        else {
            // Case C: Multiple observations - solve WLS
            
            // Extract effective subsets
            Z_sub   = Z_mat[index', .]
            y_sub   = y[index']
            K_h_sub = K_h[index']
            
            // Weighted Gram matrix: Z' diag(K_h) Z
            Gamma = cross(Z_sub, K_h_sub, Z_sub)
            
            // Cholesky inverse
            Gamma_inv = cholinv(Gamma)
            
            // Cholesky failure detection
            if (hasmissing(Gamma_inv)) {
                Estimate[j] = .
            }
            else {
                // WLS coefficient vector: (Z'KZ)^{-1} Z'Ky
                beta = Gamma_inv * cross(Z_sub, K_h_sub, y_sub)
                
                // Derivative estimate: nu! * e_nu' * beta (Eq. 12)
                Estimate[j] = factorial(deriv) * (e_vec' * beta)
            }
        }
    }
    
    return(Estimate)
}

end
