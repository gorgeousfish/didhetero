mata:

// =============================================================================
// didhetero_lpr.mata
// Local polynomial regression (LPR) estimator
//
// Provides:
//   - didhetero_polynomial()         // Polynomial design matrix construction
//   - didhetero_polynomial_inplace() // In-place polynomial basis construction
//   - _dh_bsearch_ge()              // Binary search: first index >= target
//   - _dh_bsearch_le()              // Binary search: last index <= target
//   - _dh_lpr_analytic_p1()         // Analytic 2x2 Cramer solver (p=1)
//   - didhetero_lpr()               // Main LPR estimator
//
// Requires:
//   - didhetero_kernel.mata         (didhetero_kernel_eval)
//   - didhetero_utils_formula.mata  (didhetero_selectindex)
//
// Paper reference: Section 3.1, local polynomial smoothing
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
// didhetero_polynomial_inplace()
// In-place polynomial basis construction (avoids repeated memory allocation)
// Overwrites all elements of pre-allocated result matrix
//
// Parameters:
//   u      - n x 1 colvector of centered values (x - eval_point)
//   p      - polynomial order (non-negative integer)
//   result - pre-allocated n x (p+1) matrix, modified in place
// -----------------------------------------------------------------------------
void didhetero_polynomial_inplace(real colvector u, real scalar p, real matrix result)
{
    real scalar j
    
    for (j = 0; j <= p; j++) {
        result[., j + 1] = u :^ j
    }
}


// -----------------------------------------------------------------------------
// _dh_bsearch_ge()
// Binary search: return the first index i such that sorted_x[i] >= target
// If no such index exists (all elements < target), returns rows(sorted_x)+1
// -----------------------------------------------------------------------------
real scalar _dh_bsearch_ge(real colvector sorted_x, real scalar target)
{
    real scalar lo, hi, mid
    
    lo = 1
    hi = rows(sorted_x) + 1
    while (lo < hi) {
        mid = trunc((lo + hi) / 2)
        if (sorted_x[mid] < target) lo = mid + 1
        else hi = mid
    }
    return(lo)
}


// -----------------------------------------------------------------------------
// _dh_bsearch_le()
// Binary search: return the last index i such that sorted_x[i] <= target
// If no such index exists (all elements > target), returns 0
// -----------------------------------------------------------------------------
real scalar _dh_bsearch_le(real colvector sorted_x, real scalar target)
{
    real scalar lo, hi, mid
    
    lo = 0
    hi = rows(sorted_x)
    while (lo < hi) {
        mid = trunc((lo + hi + 1) / 2)
        if (sorted_x[mid] <= target) lo = mid
        else hi = mid - 1
    }
    return(lo)
}


// -----------------------------------------------------------------------------
// _dh_lpr_analytic_p1()
// Analytic 2x2 Cramer's rule solver for local linear regression (p=1)
//
// For p=1, the normal equations form a 2x2 system:
//   [S0, S1] [beta0]   [T0]
//   [S1, S2] [beta1] = [T1]
//
// Cramer's rule gives:
//   beta0 = (S2*T0 - S1*T1) / (S0*S2 - S1^2)   (deriv==0)
//   beta1 = (S0*T1 - S1*T0) / (S0*S2 - S1^2)   (deriv==1)
//
// Falls back to missing when determinant is near-zero (degenerate case).
//
// Parameters:
//   y          - n x 1 dependent variable
//   x          - n x 1 independent variable
//   eval       - num_eval x 1 evaluation points
//   deriv      - derivative order (0 or 1)
//   kernel     - "epa" (Epanechnikov) or "gau" (Gaussian)
//   h          - scalar or num_eval x 1 bandwidth vector
//   weight     - n x 1 observation weights
//
// Returns:
//   num_eval x 1 colvector of LPR estimates
// -----------------------------------------------------------------------------
real colvector _dh_lpr_analytic_p1(
    real colvector y,
    real colvector x,
    real colvector eval,
    real scalar deriv,
    string scalar kernel,
    real colvector h_local,
    real colvector weight)
{
    real scalar n, num_eval, j, k1, k2
    real scalar S0, S1, S2, T0, T1, det
    real colvector Estimate, K_h, u_vec, d_vec
    real colvector x_sorted, y_sorted, weight_sorted
    real colvector sort_idx
    real colvector x_sub, y_sub, w_sub, K_h_sub, u_sub, d_sub
    real rowvector index
    
    // Determinant threshold for numerical stability
    // When det < threshold, the system is degenerate (too few obs or collinear)
    real scalar det_threshold
    det_threshold = 1e-10
    
    n = rows(y)
    num_eval = rows(eval)
    Estimate = J(num_eval, 1, .)
    
    if (kernel == "epa") {
        // Pre-sort for binary search optimization (compact support)
        sort_idx = order(x, 1)
        x_sorted = x[sort_idx]
        y_sorted = y[sort_idx]
        weight_sorted = weight[sort_idx]
        
        for (j = 1; j <= num_eval; j++) {
            k1 = _dh_bsearch_ge(x_sorted, eval[j] - h_local[j])
            k2 = _dh_bsearch_le(x_sorted, eval[j] + h_local[j])
            
            if (k1 > k2) {
                Estimate[j] = .
            }
            else {
                x_sub = x_sorted[|k1 \ k2|]
                y_sub = y_sorted[|k1 \ k2|]
                w_sub = weight_sorted[|k1 \ k2|]
                
                // Kernel weights including user weights
                K_h_sub = w_sub :* didhetero_kernel_eval((x_sub :- eval[j]) / h_local[j], "epa")
                
                // Filter exact boundary zeros
                index = didhetero_selectindex(K_h_sub :!= 0)
                
                if (length(index) == 0) {
                    Estimate[j] = .
                }
                else if (length(index) == 1) {
                    Estimate[j] = y_sub[index]
                }
                else if (length(index) < 4) {
                    // n_eff < 2*(p+1)=4: insufficient observations for stable p=1 LPR
                    Estimate[j] = .
                }
                else {
                    K_h_sub = K_h_sub[index']
                    y_sub   = y_sub[index']
                    d_sub   = x_sub[index'] :- eval[j]
                    
                    // Sufficient statistics for 2x2 system
                    S0 = sum(K_h_sub)
                    S1 = sum(K_h_sub :* d_sub)
                    S2 = sum(K_h_sub :* d_sub :^ 2)
                    T0 = sum(K_h_sub :* y_sub)
                    T1 = sum(K_h_sub :* y_sub :* d_sub)
                    
                    det = S0 * S2 - S1^2
                    
                    if (abs(det) < det_threshold | S0 < 1e-10) {
                        Estimate[j] = .
                    }
                    else {
                        if (deriv == 0) {
                            Estimate[j] = (S2 * T0 - S1 * T1) / det
                        }
                        else {
                            // deriv == 1
                            Estimate[j] = (S0 * T1 - S1 * T0) / det
                        }
                    }
                }
            }
        }
    }
    else {
        // Gaussian kernel (infinite support): full O(n) per eval point
        for (j = 1; j <= num_eval; j++) {
            u_vec = (x :- eval[j]) / h_local[j]
            K_h = weight :* didhetero_kernel_eval(u_vec, kernel)
            
            index = didhetero_selectindex(K_h :!= 0)
            
            if (length(index) == 0) {
                Estimate[j] = .
            }
            else if (length(index) == 1) {
                Estimate[j] = y[index]
            }
            else if (length(index) < 4) {
                // n_eff < 2*(p+1)=4: insufficient observations for stable p=1 LPR
                Estimate[j] = .
            }
            else {
                K_h_sub = K_h[index']
                y_sub   = y[index']
                d_vec   = x[index'] :- eval[j]
                
                // Sufficient statistics for 2x2 system
                S0 = sum(K_h_sub)
                S1 = sum(K_h_sub :* d_vec)
                S2 = sum(K_h_sub :* d_vec :^ 2)
                T0 = sum(K_h_sub :* y_sub)
                T1 = sum(K_h_sub :* y_sub :* d_vec)
                
                det = S0 * S2 - S1^2
                
                if (abs(det) < det_threshold | S0 < 1e-10) {
                    Estimate[j] = .
                }
                else {
                    if (deriv == 0) {
                        Estimate[j] = (S2 * T0 - S1 * T1) / det
                    }
                    else {
                        // deriv == 1
                        Estimate[j] = (S0 * T1 - S1 * T0) / det
                    }
                }
            }
        }
    }
    
    return(Estimate)
}


// -----------------------------------------------------------------------------
// didhetero_lpr()
// Local polynomial regression estimator
//
// Solves weighted least squares at each evaluation point:
//   beta(z) = (Z' K_h Z)^{-1} Z' K_h y
// Returns the nu-th derivative estimate: nu! * e_nu' * beta(z)
//
// For p=1 and deriv<=1, dispatches to analytic 2x2 Cramer solver for speed.
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
    real scalar n, num_eval, j, k1, k2, cond_num, reg_lambda
    real colvector Estimate, K_h, y_sub, K_h_sub, e_vec, beta, h_local
    real colvector x_sorted, y_sorted, weight_sorted, x_sub, w_sub
    real matrix Z_sub, Gamma, Gamma_inv, Z_work
    real rowvector index
    real colvector sort_idx
    
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
    
    // =========================================================================
    // FAST PATH: Analytic 2x2 Cramer solver for p=1, deriv<=1
    // Mathematically equivalent to Cholesky but avoids matrix construction
    // and inversion overhead. ~30-50% faster for the common CATT case.
    // =========================================================================
    if (p == 1 & deriv <= 1) {
        return(_dh_lpr_analytic_p1(y, x, eval, deriv, kernel, h_local, weight))
    }
    
    // Unit vector to extract the (deriv+1)-th coefficient from beta
    e_vec = J(p + 1, 1, 0)
    e_vec[deriv + 1] = 1
    
    // Initialize output vector
    Estimate = J(num_eval, 1, .)
    
    // Pre-sort data for Epanechnikov kernel binary search optimization
    // For Epa kernel (compact support |u|<=1), K_h((x-z)/h) != 0 iff |x-z| < h
    // Binary search on sorted data: O(log n) per eval point vs O(n) selectindex
    if (kernel == "epa") {
        sort_idx = order(x, 1)
        x_sorted = x[sort_idx]
        y_sorted = y[sort_idx]
        weight_sorted = weight[sort_idx]
    }
    
    // Pre-allocate design matrix workspace to avoid per-iteration allocation
    // Mathematical equivalence: each iteration overwrites all elements via polynomial_inplace()
    Z_work = J(n, p + 1, .)
    
    // Main loop over evaluation points
    for (j = 1; j <= num_eval; j++) {
        
        if (kernel == "epa") {
            // =========================================================
            // FAST PATH: Epanechnikov kernel with binary search
            // Find contiguous range [k1, k2] where |x - eval[j]| <= h
            // =========================================================
            k1 = _dh_bsearch_ge(x_sorted, eval[j] - h_local[j])
            k2 = _dh_bsearch_le(x_sorted, eval[j] + h_local[j])
            
            if (k1 > k2) {
                // No observations in bandwidth window
                Estimate[j] = .
            }
            else {
                // Extract contiguous sub-arrays from sorted data
                x_sub   = x_sorted[|k1 \ k2|]
                y_sub   = y_sorted[|k1 \ k2|]
                w_sub   = weight_sorted[|k1 \ k2|]
                
                // Compute kernel weights (only for observations in window)
                K_h_sub = w_sub :* didhetero_kernel_eval((x_sub :- eval[j]) / h_local[j], "epa")
                
                // Filter exact boundary zeros (|x-eval|==h => K=0) for equivalence
                index = didhetero_selectindex(K_h_sub :!= 0)
                
                if (length(index) == 0) {
                    Estimate[j] = .
                }
                else if (length(index) == 1) {
                    Estimate[j] = y_sub[index]
                }
                else if (length(index) < 2 * (p + 1)) {
                    // n_eff < 2*(p+1): insufficient observations for stable LPR
                    Estimate[j] = .
                }
                else {
                    // Build design matrix and solve WLS on the non-zero subset
                    Z_sub   = didhetero_polynomial(x_sub[index'] :- eval[j], p)
                    y_sub   = y_sub[index']
                    K_h_sub = K_h_sub[index']
                    
                    Gamma = cross(Z_sub, K_h_sub, Z_sub)
                    
                    // Condition number check and Tikhonov regularization
                    cond_num = cond(Gamma)
                    if (cond_num != . & cond_num > 1e12) {
                        reg_lambda = 1e-8 * trace(Gamma) / cols(Gamma)
                        Gamma = Gamma + reg_lambda * I(cols(Gamma))
                    }
                    
                    Gamma_inv = cholinv(Gamma)
                    
                    if (hasmissing(Gamma_inv)) {
                        Estimate[j] = .
                    }
                    else {
                        beta = Gamma_inv * cross(Z_sub, K_h_sub, y_sub)
                        Estimate[j] = factorial(deriv) * (e_vec' * beta)
                    }
                }
            }
        }
        else {
            // =========================================================
            // STANDARD PATH: Gaussian kernel (infinite support)
            // Full O(n) computation per evaluation point
            // =========================================================
            
            // Construct polynomial design matrix in-place at evaluation point eval[j]
            didhetero_polynomial_inplace(x :- eval[j], p, Z_work)
            
            // Compute kernel weights multiplied by user-specified weights
            K_h = weight :* didhetero_kernel_eval((x :- eval[j]) / h_local[j], kernel)
            
            // Select observations with positive kernel weights
            index = didhetero_selectindex(K_h :!= 0)
            
            // Compute estimate based on number of effective observations
            if (length(index) == 0) {
                Estimate[j] = .
            }
            else if (length(index) == 1) {
                Estimate[j] = y[index]
            }
            else if (length(index) < 2 * (p + 1)) {
                // n_eff < 2*(p+1): insufficient observations for stable LPR
                Estimate[j] = .
            }
            else {
                // Multiple observations: solve weighted least squares
                Z_sub   = Z_work[index', .]
                y_sub   = y[index']
                K_h_sub = K_h[index']
                
                Gamma = cross(Z_sub, K_h_sub, Z_sub)
                
                // Condition number check and Tikhonov regularization
                cond_num = cond(Gamma)
                if (cond_num != . & cond_num > 1e12) {
                    reg_lambda = 1e-8 * trace(Gamma) / cols(Gamma)
                    Gamma = Gamma + reg_lambda * I(cols(Gamma))
                }
                
                Gamma_inv = cholinv(Gamma)
                
                if (hasmissing(Gamma_inv)) {
                    Estimate[j] = .
                }
                else {
                    beta = Gamma_inv * cross(Z_sub, K_h_sub, y_sub)
                    Estimate[j] = factorial(deriv) * (e_vec' * beta)
                }
            }
        }
    }
    
    return(Estimate)
}

end
