mata:

// =============================================================================
// Generalized Propensity Score (GPS) Estimation
//
// This module implements GPS estimation via logistic regression using
// iteratively reweighted least squares (IRLS). The GPS is used to adjust
// for covariate imbalance between treated and control units in heterogeneous
// treatment effect estimation.
//
// Functions:
//   - didhetero_invlogit()      : Numerically stable inverse logit function
//   - didhetero_gps_logit()     : Logit estimation via IRLS/Newton-Raphson
//   - didhetero_gps_estimate()  : GPS estimation dispatch for treatment groups
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_invlogit()
// Numerically stable inverse logit (sigmoid) function.
// For x >= 0: 1/(1+exp(-x))  — avoids overflow in exp(x)
// For x < 0: exp(x)/(1+exp(x)) — avoids overflow in exp(-x)
//
// Args:
//   x - n x 1 column vector of linear predictors
//
// Returns:
//   n x 1 column vector of probabilities in (0, 1)
// -----------------------------------------------------------------------------
real colvector didhetero_invlogit(real colvector x)
{
    real scalar n, i
    real colvector p
    real scalar ex
    
    n = rows(x)
    p = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        if (x[i] >= 0) {
            p[i] = 1 / (1 + exp(-x[i]))
        }
        else {
            ex = exp(x[i])
            p[i] = ex / (1 + ex)
        }
    }
    
    return(p)
}

// -----------------------------------------------------------------------------
// didhetero_gps_logit()
// Logit estimation via IRLS/Newton-Raphson on subset, prediction on full sample.
// pi_hat is populated as output (Mata passes by reference).
//
// Args:
//   y_sub  - n_sub x 1 binary dependent variable (on estimation subset)
//   X_sub  - n_sub x k covariate matrix (estimation subset, includes intercept)
//   X_full - n x k covariate matrix (full sample, for prediction)
//   pi_hat - k x 1 coefficient vector (OUTPUT: populated by this function)
//
// Returns:
//   n x 1 predicted probabilities truncated to [1e-12, 1-1e-12]
// -----------------------------------------------------------------------------
real colvector didhetero_gps_logit(real colvector y_sub, real matrix X_sub,
                                    real matrix X_full, real colvector pi_hat)
{
    real scalar n_sub, n_full, k, max_iter, iter, converged
    real scalar tol
    real colvector p_sub, w, score, delta, p_full
    real matrix H
    
    n_sub  = rows(X_sub)
    n_full = rows(X_full)
    k      = cols(X_sub)
    
    // Initialize coefficients at zero
    pi_hat = J(k, 1, 0)
    
    // Algorithm parameters
    max_iter  = 25
    tol       = 1e-8
    converged = 0
    
    // Verify outcome variation for model identification
    if (sum(y_sub) == 0 | sum(y_sub) == n_sub) {
        _error("GPS logit: y_sub has no variation (all 0 or all 1)")
    }
    
    // IRLS-Newton iteration
    for (iter = 1; iter <= max_iter; iter++) {
        
        // Compute predicted probabilities
        p_sub = didhetero_invlogit(X_sub * pi_hat)
        
        // Truncate to avoid numerical boundary issues
        p_sub = rowmax((J(n_sub, 1, 1e-12), rowmin((J(n_sub, 1, 1-1e-12), p_sub))))
        
        // Construct IRLS weights
        w = p_sub :* (1 :- p_sub)
        
        // Compute score vector
        score = cross(X_sub, (y_sub - p_sub))
        
        // Compute Hessian matrix
        H = cross(X_sub, w, X_sub)
        
        // Solve for parameter update
        delta = lusolve(H, score)
        
        // Check for singular Hessian
        if (hasmissing(delta)) {
            _error("GPS logit: singular Hessian matrix at iteration " + strofreal(iter))
        }
        
        // Check for perfect or quasi-separation
        if (min(p_sub) > 1 - 1e-6 | max(p_sub) < 1e-6) {
            _error("GPS logit: perfect or quasi-separation detected at iteration " + strofreal(iter))
        }
        
        // Update parameter estimates
        pi_hat = pi_hat + delta
        
        // Check convergence criterion
        if (max(abs(delta)) < tol) {
            converged = 1
            break
        }
    }
    
    // Warn if convergence not achieved
    if (converged == 0) {
        printf("{txt}Warning: GPS logit did not converge in %g iterations\n", max_iter)
    }
    
    // Full sample prediction
    p_full = didhetero_invlogit(X_full * pi_hat)
    
    // Truncate to [1e-12, 1-1e-12]
    p_full = rowmax((J(n_full, 1, 1e-12), rowmin((J(n_full, 1, 1-1e-12), p_full))))
    
    return(p_full)
}

// -----------------------------------------------------------------------------
// didhetero_gps_estimate()
// GPS estimation dispatch for treatment groups.
//
// For each treatment group (nevertreated control) or each (g,t) pair
// (not-yet-treated control), estimates a separate logistic regression
// for the GPS and predicts on the full sample.
//
// Args:
//   data          - DidHeteroData struct containing panel data
//   gteval        - K x 2 matrix of valid (g,t) pairs
//   geval         - K_g x 1 vector of valid treatment groups
//   control_group - "nevertreated" or "notyettreated"
//   anticipation  - Anticipation periods (integer >= 0)
//   gps_coef      - (OUTPUT) Logit coefficient matrix
//
// Returns:
//   GPS matrix in long format:
//     nevertreated:  (id, g, p_hat)      - 3 columns
//     notyettreated: (id, g, t, p_hat)   - 4 columns
// -----------------------------------------------------------------------------

real matrix didhetero_gps_estimate(struct DidHeteroData scalar data,
                                    real matrix gteval,
                                    real colvector geval,
                                    string scalar control_group,
                                    real scalar anticipation,
                                    real matrix gps_coef)
{
    real matrix gps_mat, X_sub
    real colvector G, id, pi_hat, p_hat, y_sub, subset_mask, G_ord
    real rowvector idx_sub
    real scalar n, k, n_g, n_gt, i, j, g1, t1, threshold_ord, g_ord
    
    G  = data.G
    id = data.id
    n  = data.n
    k  = cols(data.X)
    
    // Initialize output matrices
    gps_mat  = J(0, 0, .)
    gps_coef = J(0, 0, .)
    
    if (control_group == "nevertreated") {
        // =====================================================================
        // Never-treated control: estimate one GPS per treatment group
        // =====================================================================
        n_g = rows(geval)
        gps_mat  = J(0, 3, .)
        gps_coef = J(0, 1 + k, .)
        
        for (i = 1; i <= n_g; i++) {
            g1 = geval[i]
            
            // Subset: units in group g1 or never-treated (G==0)
            subset_mask = (G :== g1) :| (G :== 0)
            idx_sub = didhetero_selectindex(subset_mask)
            
            if (cols(idx_sub) == 0) {
                _error("GPS estimate: empty subset for g=" + strofreal(g1))
            }
            
            // Binary outcome: 1 if treated group g1, 0 if never-treated
            y_sub = (G[idx_sub'] :== g1)
            X_sub = data.X[idx_sub', .]
            
            // Estimate logit and predict on full sample
            pi_hat = J(k, 1, .)
            p_hat = didhetero_gps_logit(y_sub, X_sub, data.X, pi_hat)
            
            // Append to gps_mat: (id, g, p_hat)
            gps_mat = gps_mat \ (id, J(n, 1, g1), p_hat)
            
            // Append to gps_coef: (g, pi_hat')
            gps_coef = gps_coef \ (g1, pi_hat')
        }
    }
    else if (control_group == "notyettreated") {
        // =====================================================================
        // Not-yet-treated control: estimate one GPS per (g,t) pair
        // =====================================================================
        n_gt = rows(gteval)
        gps_mat  = J(0, 4, .)
        gps_coef = J(0, 2 + k, .)

        G_ord = J(n, 1, 0)
        for (i = 1; i <= n; i++) {
            if (G[i] != 0) {
                g_ord = didhetero_period_ord(G[i], data.t_vals)
                if (g_ord >= .) {
                    _error("GPS estimate: group " + strofreal(G[i]) + " not found in t_vals")
                }
                G_ord[i] = g_ord
            }
        }
        
        for (j = 1; j <= n_gt; j++) {
            g1 = gteval[j, 1]
            t1 = gteval[j, 2]
            
            // Compute threshold period for not-yet-treated subset construction
            threshold_ord = didhetero_period_ord(t1, data.t_vals) + anticipation
            
            // Subset: units in group g1, never-treated (G==0), or not-yet-treated
            // units whose treatment starts strictly after the threshold period.
            subset_mask = (G :== g1) :| (G :== 0) :| (G_ord :> threshold_ord)
            idx_sub = didhetero_selectindex(subset_mask)
            
            if (cols(idx_sub) == 0) {
                _error("GPS estimate: empty subset for g=" + strofreal(g1) + " t=" + strofreal(t1))
            }
            
            // Binary outcome: 1 if treated group g1, 0 otherwise
            y_sub = (G[idx_sub'] :== g1)
            X_sub = data.X[idx_sub', .]
            
            // Estimate logit and predict on full sample
            pi_hat = J(k, 1, .)
            p_hat = didhetero_gps_logit(y_sub, X_sub, data.X, pi_hat)
            
            // Append to gps_mat: (id, g, t, p_hat)
            gps_mat = gps_mat \ (id, J(n, 1, g1), J(n, 1, t1), p_hat)
            
            // Append to gps_coef: (g, t, pi_hat')
            gps_coef = gps_coef \ (g1, t1, pi_hat')
        }
    }
    else {
        _error("GPS estimate: invalid control_group '" + control_group + "'")
    }
    
    return(gps_mat)
}

end
