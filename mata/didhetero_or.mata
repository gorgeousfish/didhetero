mata:

// =============================================================================
// Outcome Regression (OR) Estimation
//
// This module implements outcome regression estimation for heterogeneous
// treatment effects in difference-in-differences designs.
//
// Functions:
//   didhetero_period_col()  - Locate period index in time vector
//   didhetero_or_ols()      - OLS estimation with full-sample prediction
//   didhetero_or_estimate() - Main estimation routine
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_period_col()
//
// Locate the column index corresponding to a specific time period.
//
// Parameters:
//   t      - scalar, target period value
//   t_vals - colvector, sorted unique period values
//
// Returns:
//   scalar, 1-based column index in the wide-format outcome matrix
// -----------------------------------------------------------------------------
real scalar didhetero_period_col(real scalar t, real colvector t_vals)
{
    real scalar i

    i = didhetero_period_ord(t, t_vals)
    if (i < .) return(i)

    _error("OR estimate: period " + strofreal(t) + " not found in t_vals")
    return(.)
}

// -----------------------------------------------------------------------------
// didhetero_or_ols()
//
// Estimate OLS coefficients on a control group subset and compute
// predicted values for the full sample.
//
// Parameters:
//   y_sub    - colvector, dependent variable for estimation subset
//   X_sub    - matrix, covariates for estimation subset (includes intercept)
//   X_full   - matrix, covariates for full sample
//   beta_hat - colvector, estimated coefficients (output)
//
// Returns:
//   colvector, predicted values for the full sample
// -----------------------------------------------------------------------------
real colvector didhetero_or_ols(real colvector y_sub, real matrix X_sub,
                                real matrix X_full, real colvector beta_hat)
{
    real matrix XtX
    real colvector Xty
    
    // Solve OLS normal equations: beta = (X'X)^{-1} X'y
    XtX = cross(X_sub, X_sub)
    Xty = cross(X_sub, y_sub)
    beta_hat = lusolve(XtX, Xty)
    
    // Check for singular design matrix
    if (hasmissing(beta_hat)) {
        _error("OR OLS: singular design matrix (X'X)")
    }
    
    // Compute predictions for full sample
    return(X_full * beta_hat)
}

// -----------------------------------------------------------------------------
// didhetero_or_estimate()
//
// Compute outcome regression estimates for all (group, time) pairs.
// Supports both never-treated and not-yet-treated control groups.
//
// Parameters:
//   data          - struct DidHeteroData, estimation data
//   gteval        - matrix, K x 2 matrix of (group, time) pairs
//   control_group - string, control group type ("nevertreated" or "notyettreated")
//   anticipation  - scalar, number of anticipation periods
//   or_coef       - matrix, OLS coefficients by (group, time) (output)
//
// Returns:
//   matrix, estimation results in long format (id, group, time, estimate)
// -----------------------------------------------------------------------------
real matrix didhetero_or_estimate(struct DidHeteroData scalar data,
                                   real matrix gteval,
                                   string scalar control_group,
                                   real scalar anticipation,
                                   real matrix or_coef)
{
    real matrix or_mat, X_sub
    real colvector G, id, beta_hat, m_hat, y_diff, y_sub, subset_mask, G_ord
    real rowvector idx_sub, idx_sub_nev
    real scalar n, k, K, j, i, g1, t1, col_t, col_base, threshold_ord, base_label
    
    G  = data.G
    id = data.id
    n  = data.n
    k  = cols(data.X)
    K  = rows(gteval)
    
    // Initialize output matrices
    or_mat  = J(0, 4, .)
    or_coef = J(0, 2 + k, .)
    
    // Never-treated: control group subset is constant across all (g,t) pairs
    if (control_group == "nevertreated") {
        // Identify never-treated units (G == 0)
        idx_sub_nev = didhetero_selectindex(G :== 0)
        
        if (cols(idx_sub_nev) == 0) {
            _error("OR estimate: no never-treated units (G==0)")
        }
        
        // Extract covariates for control group (constant across pairs)
        X_sub = data.X[idx_sub_nev', .]
        
        // Loop over all (g,t) pairs
        for (j = 1; j <= K; j++) {
            g1 = gteval[j, 1]
            t1 = gteval[j, 2]
            
            // Construct dependent variable: Y_t - Y_{g-anticipation-1}
            col_t    = didhetero_period_col(t1, data.t_vals)
            base_label = didhetero_period_at(
                didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
                data.t_vals,
                "OR estimate base period")
            col_base = didhetero_period_col(base_label, data.t_vals)
            y_diff   = data.Y_wide[., col_t] - data.Y_wide[., col_base]
            
            // Extract dependent variable for control group
            y_sub = y_diff[idx_sub_nev']
            
            // Estimate on control group, predict for full sample
            beta_hat = J(k, 1, .)
            m_hat = didhetero_or_ols(y_sub, X_sub, data.X, beta_hat)
            
            // Store coefficient estimates
            or_coef = or_coef \ (g1, t1, beta_hat')
            
            // Append results to output matrix
            or_mat = or_mat \ (id, J(n, 1, g1), J(n, 1, t1), m_hat)
        }
    }
    // Not-yet-treated: control group subset varies by (g,t) pair
    else if (control_group == "notyettreated") {
        G_ord = J(n, 1, 0)
        for (i = 1; i <= n; i++) {
            if (G[i] != 0) {
                G_ord[i] = didhetero_period_ord(G[i], data.t_vals)
                if (G_ord[i] >= .) {
                    _error("OR estimate: group " + strofreal(G[i]) + " not found in t_vals")
                }
            }
        }

        for (j = 1; j <= K; j++) {
            g1 = gteval[j, 1]
            t1 = gteval[j, 2]
            
            // Compute threshold for not-yet-treated condition
            threshold_ord = didhetero_period_ord(t1, data.t_vals) + anticipation
            
            // Identify control units: never-treated or not-yet-treated
            subset_mask = (G :== 0) :| (G_ord :> threshold_ord)
            idx_sub = didhetero_selectindex(subset_mask)
            
            if (cols(idx_sub) == 0) {
                _error("OR estimate: empty subset for (g,t)=(" + strofreal(g1) + "," + strofreal(t1) + ")")
            }
            
            // Construct dependent variable: Y_t - Y_{g-anticipation-1}
            col_t    = didhetero_period_col(t1, data.t_vals)
            base_label = didhetero_period_at(
                didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
                data.t_vals,
                "OR estimate base period")
            col_base = didhetero_period_col(base_label, data.t_vals)
            y_diff   = data.Y_wide[., col_t] - data.Y_wide[., col_base]
            
            // Extract data for control group
            y_sub = y_diff[idx_sub']
            X_sub = data.X[idx_sub', .]
            
            // Estimate on control group, predict for full sample
            beta_hat = J(k, 1, .)
            m_hat = didhetero_or_ols(y_sub, X_sub, data.X, beta_hat)
            
            // Store coefficient estimates
            or_coef = or_coef \ (g1, t1, beta_hat')
            
            // Append results to output matrix
            or_mat = or_mat \ (id, J(n, 1, g1), J(n, 1, t1), m_hat)
        }
    }
    else {
        _error("OR estimate: invalid control_group '" + control_group + "'")
    }
    
    return(or_mat)
}

end
