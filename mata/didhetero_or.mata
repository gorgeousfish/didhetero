mata:

// =============================================================================
// didhetero_or.mata
// OR (Outcome Regression) estimation via OLS
//
// Functions:
//   1. didhetero_period_col()    - Find column index of period in Y_wide
//   2. didhetero_or_ols()        - OLS core: estimate on subset, predict full
//   3. didhetero_or_estimate()   - OR dispatch (nevertreated/notyettreated)
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Section 4.2.1 (OR estimation via OLS)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_period_col()
// Find 1-based column index of period t in t_vals vector.
//
// Args:
//   t      - target period value
//   t_vals - T x 1 sorted vector of unique period values
//
// Returns:
//   1-based column index
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
// OLS estimation on subset, prediction on full sample.
// beta_hat is populated as output (Mata passes by reference).
//
// Args:
//   y_sub    - n_sub x 1 dependent variable (estimation subset)
//   X_sub    - n_sub x k covariate matrix (estimation subset, with intercept)
//   X_full   - n x k covariate matrix (full sample, for prediction)
//   beta_hat - k x 1 coefficient vector (OUTPUT)
//
// Returns:
//   n x 1 vector of predicted values (m_hat = X_full * beta_hat)
//
// Paper ref: Section 4.2.1, outcome regression via OLS
// -----------------------------------------------------------------------------
real colvector didhetero_or_ols(real colvector y_sub, real matrix X_sub,
                                real matrix X_full, real colvector beta_hat)
{
    real matrix XtX
    real colvector Xty
    
    // OLS normal equations: beta = (X'X)^{-1} X'y
    XtX = cross(X_sub, X_sub)
    Xty = cross(X_sub, y_sub)
    beta_hat = lusolve(XtX, Xty)
    
    // Check for singular design matrix
    if (hasmissing(beta_hat)) {
        _error("OR OLS: singular design matrix (X'X)")
    }
    
    // Full sample prediction
    return(X_full * beta_hat)
}

// -----------------------------------------------------------------------------
// didhetero_or_estimate()
// OR estimation dispatch for all (g,t) pairs.
// Always loops over (g,t) pairs regardless of control_group.
//
// Args:
//   data          - DidHeteroData struct
//   gteval        - K x 2 matrix of valid (g,t) pairs
//   control_group - "nevertreated" or "notyettreated"
//   anticipation  - anticipation periods (integer >= 0)
//   or_coef       - (OUTPUT) K x (2+k) matrix of OLS coefficients
//
// Returns:
//   OR result matrix (id, g, t, est), always 4 columns, long format
//
// Paper ref: Section 4.2.1, OR estimation dispatch
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
    
    // =====================================================================
    // nevertreated: subset defined OUTSIDE loop
    // =====================================================================
    if (control_group == "nevertreated") {
        // Subset: only never-treated units (G == 0)
        idx_sub_nev = didhetero_selectindex(G :== 0)
        
        if (cols(idx_sub_nev) == 0) {
            _error("OR estimate: no never-treated units (G==0)")
        }
        
        // Extract subset covariates once (fixed across all (g,t))
        X_sub = data.X[idx_sub_nev', .]
        
        // Loop over all (g,t) pairs
        for (j = 1; j <= K; j++) {
            g1 = gteval[j, 1]
            t1 = gteval[j, 2]
            
            // Dependent variable: Y_t - Y_{g-delta-1}
            col_t    = didhetero_period_col(t1, data.t_vals)
            base_label = didhetero_period_at(
                didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
                data.t_vals,
                "OR estimate base period")
            col_base = didhetero_period_col(base_label, data.t_vals)
            y_diff   = data.Y_wide[., col_t] - data.Y_wide[., col_base]
            
            // Subset dependent variable
            y_sub = y_diff[idx_sub_nev']
            
            // OLS: estimate on subset, predict on full sample
            beta_hat = J(k, 1, .)
            m_hat = didhetero_or_ols(y_sub, X_sub, data.X, beta_hat)
            
            // Collect OLS coefficients
            or_coef = or_coef \ (g1, t1, beta_hat')
            
            // Append (id, g, t, est) block
            or_mat = or_mat \ (id, J(n, 1, g1), J(n, 1, t1), m_hat)
        }
    }
    // =====================================================================
    // notyettreated: subset defined INSIDE loop
    // =====================================================================
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
            
            // Appendix D uses D_{t+delta} = 0, which depends on the ordinal
            // position of t within the observed time support.
            threshold_ord = didhetero_period_ord(t1, data.t_vals) + anticipation
            
            // Subset: never-treated OR not-yet-treated
            subset_mask = (G :== 0) :| (G_ord :> threshold_ord)
            idx_sub = didhetero_selectindex(subset_mask)
            
            if (cols(idx_sub) == 0) {
                _error("OR estimate: empty subset for (g,t)=(" + strofreal(g1) + "," + strofreal(t1) + ")")
            }
            
            // Dependent variable: Y_t - Y_{g-delta-1}
            col_t    = didhetero_period_col(t1, data.t_vals)
            base_label = didhetero_period_at(
                didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
                data.t_vals,
                "OR estimate base period")
            col_base = didhetero_period_col(base_label, data.t_vals)
            y_diff   = data.Y_wide[., col_t] - data.Y_wide[., col_base]
            
            // Subset
            y_sub = y_diff[idx_sub']
            X_sub = data.X[idx_sub', .]
            
            // OLS: estimate on subset, predict on full sample
            beta_hat = J(k, 1, .)
            m_hat = didhetero_or_ols(y_sub, X_sub, data.X, beta_hat)
            
            // Collect OLS coefficients
            or_coef = or_coef \ (g1, t1, beta_hat')
            
            // Append (id, g, t, est) block
            or_mat = or_mat \ (id, J(n, 1, g1), J(n, 1, t1), m_hat)
        }
    }
    else {
        _error("OR estimate: invalid control_group '" + control_group + "'")
    }
    
    return(or_mat)
}

end
