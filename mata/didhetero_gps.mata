mata:

// =============================================================================
// didhetero_gps.mata
// Generalized Propensity Score (GPS) estimation via logistic regression
//
// Provides:
//   - didhetero_invlogit()      // Numerically stable inverse logit function
//   - didhetero_gps_logit()     // Logit estimation via IRLS/Newton-Raphson
//   - didhetero_gps_estimate()  // GPS estimation dispatch for treatment groups
//
// Requires:
//   - didhetero_types.mata          (DidHeteroData struct)
//   - didhetero_utils_formula.mata  (didhetero_selectindex)
//   - didhetero_utils_domain.mata   (didhetero_period_ord)
//
// Paper reference: Section 2.3, generalized propensity score
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
//   strict - (optional) 0=warn on non-convergence (default), 1=error on non-convergence
//
// Returns:
//   n x 1 predicted probabilities truncated to [1e-12, 1-1e-12]
// -----------------------------------------------------------------------------
real colvector didhetero_gps_logit(real colvector y_sub, real matrix X_sub,
                                    real matrix X_full, real colvector pi_hat,
                                    | real scalar strict, real rowvector diag_row,
                                      real scalar verbose)
{
    real scalar n_sub, n_full, k, max_iter, iter, converged
    real scalar tol, tol_grad, tol_rel_ll
    real scalar ll_old, ll_new, rel_improve
    real colvector p_sub, w, score, delta, p_full
    real matrix H, iter_history
    
    n_sub  = rows(X_sub)
    n_full = rows(X_full)
    k      = cols(X_sub)
    
    // Default: soft failure (backward compatible)
    if (args() < 5) strict = 0
    if (args() < 7) verbose = 0
    
    // Initialize coefficients at zero
    pi_hat = J(k, 1, 0)
    
    // Algorithm parameters
    max_iter   = 25
    tol        = 1e-8          // Parameter step convergence threshold
    tol_grad   = 1e-8          // Gradient norm convergence threshold
    tol_rel_ll = 1e-10         // Relative log-likelihood improvement threshold
    converged  = 0
    ll_old     = .             // Initialize log-likelihood tracker
    
    // Iteration history: columns = (iter, max_abs_delta, max_abs_score, log_likelihood)
    iter_history = J(max_iter, 4, .)
    
    // Verify outcome variation for model identification
    if (sum(y_sub) == 0 | sum(y_sub) == n_sub) {
        _error("GPS logit: y_sub has no variation (all 0 or all 1)")
    }
    
    // ========== IRLS-Newton iteration loop ==========
    for (iter = 1; iter <= max_iter; iter++) {
        
        // Step 1: Compute predicted probabilities
        p_sub = didhetero_invlogit(X_sub * pi_hat)
        
        // Step 2: Numerical stability truncation
        p_sub = rowmax((J(n_sub, 1, 1e-12), rowmin((J(n_sub, 1, 1-1e-12), p_sub))))
        
        // Step 3: Compute log-likelihood for convergence monitoring
        ll_new = sum(y_sub :* ln(p_sub) + (1 :- y_sub) :* ln(1 :- p_sub))
        
        // Step 4: Construct IRLS weights
        w = p_sub :* (1 :- p_sub)
        
        // Step 5: Compute Score vector (gradient)
        score = cross(X_sub, (y_sub - p_sub))
        
        // Step 6: Compute Hessian matrix
        H = cross(X_sub, w, X_sub)
        
        // Step 7: Solve linear system H * delta = score
        delta = lusolve(H, score)
        
        // Step 8: Singular matrix check
        if (hasmissing(delta)) {
            _error("GPS logit: singular Hessian matrix at iteration " + strofreal(iter))
        }
        
        // Step 9: Perfect/quasi-separation check
        if (min(p_sub) > 1 - 1e-6 | max(p_sub) < 1e-6) {
            _error("GPS logit: perfect or quasi-separation detected at iteration " + strofreal(iter))
        }
        
        // Step 10: Parameter update
        pi_hat = pi_hat + delta
        
        // Step 11: Record iteration history
        iter_history[iter, 1] = iter
        iter_history[iter, 2] = max(abs(delta))
        iter_history[iter, 3] = max(abs(score))
        iter_history[iter, 4] = ll_new
        
        // ========== Enhanced convergence criteria ==========
        // Criterion 1: Parameter step size AND gradient norm (MLE first-order condition)
        // Criterion 2: Relative log-likelihood improvement (practical convergence)
        rel_improve = (ll_old != . ? abs(ll_new - ll_old) / max((1, abs(ll_old))) : .)
        
        if (max(abs(delta)) < tol & max(abs(score)) < tol_grad) {
            converged = 1
            break
        }
        if (rel_improve != . & rel_improve < tol_rel_ll & iter > 3) {
            converged = 1
            break
        }
        
        ll_old = ll_new
    }
    
    // ========== Post-iteration handling ==========
    if (converged == 0) {
        real scalar _final_grad_norm, _final_step_norm, _final_rel_ll
        _final_step_norm = iter_history[max_iter, 2]
        _final_grad_norm = iter_history[max_iter, 3]
        _final_rel_ll = (max_iter >= 2 ? 
            abs(iter_history[max_iter, 4] - iter_history[max_iter-1, 4]) / 
            max((1, abs(iter_history[max_iter-1, 4]))) : .)
        
        if (strict) {
            errprintf("\n{err}Error: GPS logit did not converge in %g iterations\n", max_iter)
            errprintf("{err}  Diagnostics at final iteration:\n")
            errprintf("{err}    Gradient norm (max|score|):      %12.6e\n", _final_grad_norm)
            errprintf("{err}    Parameter step (max|delta|):     %12.6e\n", _final_step_norm)
            errprintf("{err}    Relative log-likelihood change:  %12.6e\n", _final_rel_ll)
            errprintf("{err}    Log-likelihood:                  %12.4f\n", iter_history[max_iter, 4])
            errprintf("{err}  Possible causes:\n")
            errprintf("{err}    - Complete or quasi-separation in covariates\n")
            errprintf("{err}    - Multicollinearity in the covariate matrix\n")
            errprintf("{err}    - Model misspecification (overlap violation)\n")
            errprintf("{err}  Suggestion: check data for perfect predictors of treatment assignment\n")
            _error(430, "GPS logit failed to converge (strict mode)")
        }
        
        printf("{txt}Warning: GPS logit did not converge in %g iterations\n", max_iter)
        if (verbose) {
            printf("{txt}  Diagnostics at final iteration:\n")
            printf("{txt}    Gradient norm (max|score|):      %12.6e\n", _final_grad_norm)
            printf("{txt}    Parameter step (max|delta|):     %12.6e\n", _final_step_norm)
            printf("{txt}    Relative log-likelihood change:  %12.6e\n", _final_rel_ll)
            printf("{txt}    Log-likelihood:                  %12.4f\n", iter_history[max_iter, 4])
            if (max_iter >= 3) {
                printf("{txt}    Last 3 iterations ll: %12.4f -> %12.4f -> %12.4f\n",
                    iter_history[max_iter-2, 4], iter_history[max_iter-1, 4], iter_history[max_iter, 4])
            }
            printf("{txt}  Note: Results may be unreliable. Check for complete/quasi-separation in data.\n")
            printf("{txt}  Hint: Use gpsstrict option to enforce convergence as an error.\n")
        }
    }
    
    // Full-sample prediction
    p_full = didhetero_invlogit(X_full * pi_hat)
    
    // Full-sample probability truncation
    p_full = rowmax((J(n_full, 1, 1e-12), rowmin((J(n_full, 1, 1-1e-12), p_full))))
    
    // ========== Post-convergence separation diagnostic ==========
    // Check for quasi-separation in full-sample predictions
    // Theoretical basis: Assumption 2.5 requires p_{g,t}(X) bounded away from 0 and 1.
    // Quasi-separated fits violate this assumption and produce unreliable R_{i,g,t} weights.
    {
        real scalar _n_extreme, _sep_threshold, _cond_number
        _sep_threshold = 0.01
        _n_extreme = sum(p_full :< _sep_threshold) + sum(p_full :> (1 - _sep_threshold))
        if (_n_extreme > 0 & verbose) {
            printf("{txt}  GPS quasi-separation detail: %g/%g obs extreme\n",
                _n_extreme, n_full)
            printf("{txt}         %5.1f pct have predicted P outside [%5.3f, %5.3f]\n",
                100 * _n_extreme / n_full, _sep_threshold, 1 - _sep_threshold)
            if (_n_extreme > 0.25 * n_full) {
                printf("{res}         CAUTION: >25 pct extreme predictions suggest model misspecification or overlap violation\n")
            }
        }

        // Compute Hessian condition number for diagnostics
        _cond_number = cond(H)

        // Assemble diagnostics row: [converged, iterations, max_gradient, ll_final, n_extreme, cond_number]
        diag_row = (converged, min((iter, max_iter)), max(abs(score)), ll_new, _n_extreme, _cond_number)
    }
    
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
    real rowvector idx_sub, diag_row
    real scalar n, k, n_g, n_gt, i, j, g1, t1, threshold_ord, g_ord
    real scalar _gps_strict
    real scalar _trim_lo, _trim_hi, _n_below, _n_above, _total_trimmed
    
    G  = data.G
    id = data.id
    n  = data.n
    k  = cols(data.X)
    
    // Resolve strict mode from struct (default 0 = soft failure)
    _gps_strict = (data.gps_strict >= . ? 0 : data.gps_strict)
    
    // Resolve GPS trim bounds (default: no trimming)
    _trim_lo = (data.gps_trim_lo >= . ? 0 : data.gps_trim_lo)
    _trim_hi = (data.gps_trim_hi >= . ? 1 : data.gps_trim_hi)
    _total_trimmed = 0
    
    // Initialize output matrices
    gps_mat  = J(0, 0, .)
    gps_coef = J(0, 0, .)
    
    // Quasi-separation aggregation tracking
    real scalar _qs_count, _qs_total, _qs_max_extreme
    _qs_count = 0
    _qs_total = 0
    _qs_max_extreme = 0
    
    if (control_group == "nevertreated") {
        // =====================================================================
        // Never-treated control: estimate one GPS per treatment group
        // =====================================================================
        n_g = rows(geval)
        gps_mat  = J(0, 3, .)
        gps_coef = J(0, 1 + k, .)
        data.gps_diagnostics = J(n_g, 6, .)
        
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
            diag_row = J(1, 6, .)
            p_hat = didhetero_gps_logit(y_sub, X_sub, data.X, pi_hat, _gps_strict, diag_row, data.verbose)
            
            // Store diagnostics for this group
            data.gps_diagnostics[i, .] = diag_row
            
            // Track quasi-separation for aggregated warning
            _qs_total = _qs_total + 1
            if (diag_row[5] > 0) {
                _qs_count = _qs_count + 1
                if (diag_row[5] > _qs_max_extreme) _qs_max_extreme = diag_row[5]
            }
            
            // Verbose output: GPS convergence diagnostics per group
            if (data.verbose) {
                printf("{txt}  GPS (g=%g): converged=%s, iter=%g, ||grad||=%9.2e, LL=%12.6f\n",
                    g1, (diag_row[1] ? "yes" : "no"), diag_row[2], diag_row[3], diag_row[4])
            }
            
            // Apply GPS trimming if specified
            if (_trim_lo > 0 | _trim_hi < 1) {
                _n_below = sum(p_hat :< _trim_lo)
                _n_above = sum(p_hat :> _trim_hi)
                p_hat = rowmax((p_hat, J(rows(p_hat), 1, _trim_lo)))
                p_hat = rowmin((p_hat, J(rows(p_hat), 1, _trim_hi)))
                _total_trimmed = _total_trimmed + _n_below + _n_above
            }
            
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
        data.gps_diagnostics = J(n_gt, 6, .)

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
            diag_row = J(1, 6, .)
            p_hat = didhetero_gps_logit(y_sub, X_sub, data.X, pi_hat, _gps_strict, diag_row, data.verbose)
            
            // Store diagnostics for this (g,t) pair
            data.gps_diagnostics[j, .] = diag_row
            
            // Track quasi-separation for aggregated warning
            _qs_total = _qs_total + 1
            if (diag_row[5] > 0) {
                _qs_count = _qs_count + 1
                if (diag_row[5] > _qs_max_extreme) _qs_max_extreme = diag_row[5]
            }
            
            // Verbose output: GPS convergence diagnostics per (g,t) pair
            if (data.verbose) {
                printf("{txt}  GPS (g=%g, t=%g): converged=%s, iter=%g, ||grad||=%9.2e, LL=%12.6f\n",
                    g1, t1, (diag_row[1] ? "yes" : "no"), diag_row[2], diag_row[3], diag_row[4])
            }
            
            // Apply GPS trimming if specified
            if (_trim_lo > 0 | _trim_hi < 1) {
                _n_below = sum(p_hat :< _trim_lo)
                _n_above = sum(p_hat :> _trim_hi)
                p_hat = rowmax((p_hat, J(rows(p_hat), 1, _trim_lo)))
                p_hat = rowmin((p_hat, J(rows(p_hat), 1, _trim_hi)))
                _total_trimmed = _total_trimmed + _n_below + _n_above
            }
            
            // Append to gps_mat: (id, g, t, p_hat)
            gps_mat = gps_mat \ (id, J(n, 1, g1), J(n, 1, t1), p_hat)
            
            // Append to gps_coef: (g, t, pi_hat')
            gps_coef = gps_coef \ (g1, t1, pi_hat')
        }
    }
    else {
        _error("GPS estimate: invalid control_group '" + control_group + "'")
    }
    
    // Aggregated quasi-separation warning (single summary line)
    if (_qs_count > 0) {
        printf("{txt}Warning: GPS quasi-separation in %g of %g (g,t) pairs (max: %g/%g obs extreme)\n",
            _qs_count, _qs_total, _qs_max_extreme, n)
    }
    
    // Store total trimmed count in data struct
    data.gps_n_trimmed = _total_trimmed
    
    // Display trimming diagnostics
    if (_total_trimmed > 0) {
        real scalar _total_obs, _trim_pct
        if (control_group == "nevertreated") {
            _total_obs = n * n_g
        }
        else {
            _total_obs = n * n_gt
        }
        _trim_pct = 100 * _total_trimmed / _total_obs
        printf("{txt}Note: %g observations had GPS values trimmed to [%g, %g]\n",
            _total_trimmed, _trim_lo, _trim_hi)
        if (_trim_pct > 10) {
            printf("{err}Warning: %5.1f pct of GPS values were trimmed (>10 pct).\n", _trim_pct)
            printf("{err}         Consider checking the overlap assumption (Assumption 3).\n")
        }
    }
    
    return(gps_mat)
}

end
