// =============================================================================
// Mata aggregation module for heterogeneous treatment effect estimation
//
// This module implements the construction and filtering of (g, t, eval) triples
// for four aggregation types: dynamic, group, calendar, and simple. The aggregation
// combines conditional average treatment effects (CATT) across groups (g),
// time periods (t), and evaluation points (eval).
// =============================================================================

mata:

// -----------------------------------------------------------------------------
// _gteeval_lt_gbar()
//
// Comparison helper treating missing gbar as infinity. In Mata, the expression
// val < . evaluates to missing rather than 1; this function ensures a binary
// return value.
//
// Arguments:
//   val  - scalar value to compare
//   gbar - upper bound (missing denotes infinity)
//
// Returns:
//   1 if val < gbar or gbar is missing; 0 otherwise
// -----------------------------------------------------------------------------
real scalar _gteeval_lt_gbar(real scalar val, real scalar gbar)
{
    if (gbar >= .) return(1)
    return(val < gbar)
}

// -----------------------------------------------------------------------------
// _aggte_time_support_from_gteval()
//
// Construct the ordered time support from the (g, t) pairs in gteval. Returns
// the sorted union of all unique group and time period values.
// -----------------------------------------------------------------------------
real colvector _aggte_time_support_from_gteval(real matrix gteval)
{
    real colvector vals
    real scalar i

    vals = J(rows(gteval) * 2, 1, .)
    for (i = 1; i <= rows(gteval); i++) {
        vals[2 * i - 1] = gteval[i, 1]
        vals[2 * i]     = gteval[i, 2]
    }

    return(sort(uniqrows(vals), 1))
}

// -----------------------------------------------------------------------------
// _aggte_dynamic_event_time()
//
// Compute the event time as the ordinal distance between group adoption time
// g and evaluation period t on the ordered time support.
// -----------------------------------------------------------------------------
real scalar _aggte_dynamic_event_time(
    real scalar g0,
    real scalar t0,
    real colvector t_support,
    string scalar context)
{
    real scalar g_ord, t_ord

    g_ord = didhetero_period_ord(g0, t_support)
    t_ord = didhetero_period_ord(t0, t_support)

    if (g_ord >= . | t_ord >= .) {
        _error(3498, context + ": found a (g,t) pair outside the ordered time support")
    }

    return(t_ord - g_ord)
}

// -----------------------------------------------------------------------------
// didhetero_build_gteeval()
//
// Construct the (g, t, eval) triple matrix for aggregation by filtering (g, t)
// pairs according to aggregation type and evaluation points. For non-dynamic
// aggregation types, pre-filters to retain only post-treatment periods (t - g >= 0).
// Then for each (g, t) pair and evaluation point, applies type-specific inclusion
// criteria.
//
// Arguments:
//   gteval   - K x 2 matrix of (g, t) pairs
//   gbar     - upper bound for control group (missing denotes infinity)
//   type     - aggregation type: "dynamic", "group", "calendar", or "simple"
//   eval_pts - column vector of evaluation points (missing scalar for "simple")
//
// Returns:
//   M x 3 matrix of (g, t, eval) triples satisfying inclusion conditions,
//   or J(0, 3, .) if no triples match
// -----------------------------------------------------------------------------
real matrix didhetero_build_gteeval(
    real matrix gteval,
    real scalar gbar,
    string scalar type,
    real colvector eval_pts)
{
    real matrix     gteval0, all_gteeval
    real colvector  keep, dyn_support
    real scalar     num_gt, num_eval, id_gt, id_eval
    real scalar     g0, t0, e0, i, dyn_e, t_ord, gbar_ord

    // Input validation
    if (type != "dynamic" & type != "group" &
        type != "calendar" & type != "simple") {
        _error(3498, "invalid aggregation type: " + type)
    }

    // Pre-filter: exclude pre-treatment periods (t - g < 0) for group/calendar/simple
    if (type != "dynamic") {
        keep = J(rows(gteval), 1, 0)
        for (i = 1; i <= rows(gteval); i++) {
            if (gteval[i, 2] - gteval[i, 1] >= 0) {
                keep[i] = 1
            }
        }
        gteval0 = select(gteval, keep)
    }
    else {
        gteval0 = gteval
    }

    if (rows(gteval0) == 0) return(J(0, 3, .))

    if (type == "dynamic") {
        dyn_support = _aggte_time_support_from_gteval(gteval0)
        if (gbar >= .) gbar_ord = .
        else {
            gbar_ord = didhetero_period_ord(gbar, dyn_support)
            if (gbar_ord >= .) {
                _error(3498, "aggte_gt: failed to locate gbar on the ordered time support")
            }
        }
    }

    // Initialize
    num_gt   = rows(gteval0)
    num_eval = rows(eval_pts)
    all_gteeval = J(0, 3, .)

    // Build triples by nested loop over (g,t) pairs and evaluation points
    for (id_gt = 1; id_gt <= num_gt; id_gt++) {

        g0 = gteval0[id_gt, 1]
        t0 = gteval0[id_gt, 2]

        for (id_eval = 1; id_eval <= num_eval; id_eval++) {

            e0 = eval_pts[id_eval]

            // Dynamic aggregation: select pairs matching event time
            if (type == "dynamic") {
                dyn_e = _aggte_dynamic_event_time(
                    g0, t0, dyn_support, "didhetero_build_gteeval()")
                t_ord = didhetero_period_ord(t0, dyn_support)
                if (e0 == dyn_e & _gteeval_lt_gbar(t_ord, gbar_ord)) {
                    all_gteeval = all_gteeval \ (g0, t0, e0)
                }
            }
            // Group aggregation: select pairs matching group identifier
            else if (type == "group") {
                if (g0 == e0 & t0 >= e0 & _gteeval_lt_gbar(t0, gbar)) {
                    all_gteeval = all_gteeval \ (g0, t0, e0)
                }
            }
            // Calendar aggregation: select pairs matching calendar time
            else if (type == "calendar") {
                if (g0 <= e0 & t0 == e0) {
                    all_gteeval = all_gteeval \ (g0, t0, e0)
                }
            }
            // Simple aggregation: select all post-treatment pairs
            else if (type == "simple") {
                if (g0 <= t0 & _gteeval_lt_gbar(t0, gbar)) {
                    all_gteeval = all_gteeval \ (g0, t0, e0)
                }
            }
        }
    }

    return(all_gteeval)
}

// -----------------------------------------------------------------------------
// didhetero_aggte_weights()
//
// Compute aggregation weights for a given evaluation point. For dynamic,
// calendar, and simple aggregation, weights are proportional to the conditional
// group density: w_{g,t}(z_r) = mu_G_g(z_r) / sum(mu_G_g(z_r)). For group
// aggregation, weights are uniform: w_{g,t} = 1 / T_post(g').
//
// Arguments:
//   mu_G_g_sub   - R x num_gte conditional group density submatrix
//   gteeval      - M x 3 matrix of (g, t, eval) triples
//   type         - aggregation type
//   gbar         - upper bound for control group (missing denotes infinity)
//   num_zeval    - number of evaluation points R
//   num_gte      - number of (g,t) pairs
//
// Outputs (passed by reference):
//   aggte_weight - R x num_gte weight matrix
//   aggte_kappa  - R x 1 normalization constant
// -----------------------------------------------------------------------------
void didhetero_aggte_weights(
    real matrix mu_G_g_sub,
    real matrix gteeval,
    string scalar type,
    real scalar gbar,
    real scalar num_zeval,
    real scalar num_gte,
    real matrix aggte_weight,
    real colvector aggte_kappa)
{
    real colvector mu_rowsum
    real scalar k, id_gte, e1, T_post, j

    // Input validation
    if (type != "dynamic" & type != "group" &
        type != "calendar" & type != "simple") {
        _error(3498, "invalid aggregation type: " + type)
    }

    // Initialize outputs
    aggte_weight = J(num_zeval, num_gte, .)
    aggte_kappa = J(num_zeval, 1, .)

    // Precompute row sums of mu_G_g_sub
    mu_rowsum = J(num_zeval, 1, 0)
    for (k = 1; k <= num_gte; k++) {
        mu_rowsum = mu_rowsum + mu_G_g_sub[., k]
    }

    // Weight calculation main loop
    for (id_gte = 1; id_gte <= num_gte; id_gte++) {

        if (type == "dynamic" | type == "calendar") {
            // Conditional group probability weights
            aggte_weight[., id_gte] = mu_G_g_sub[., id_gte] :/ mu_rowsum
        }
        else if (type == "group") {
            // Known uniform weights 1/T_post
            e1 = gteeval[id_gte, 3]
            T_post = 0
            for (j = 1; j <= rows(gteeval); j++) {
                if (gteeval[j, 1] == e1 & e1 <= gteeval[j, 2] &
                    _gteeval_lt_gbar(gteeval[j, 2], gbar)) {
                    T_post = T_post + 1
                }
            }
            if (T_post > 0) {
                aggte_weight[., id_gte] = J(num_zeval, 1, 1/T_post)
            }
            else {
                aggte_weight[., id_gte] = J(num_zeval, 1, 0)
            }
        }
        else if (type == "simple") {
            // Same formula as dynamic/calendar
            aggte_weight[., id_gte] = mu_G_g_sub[., id_gte] :/ mu_rowsum
        }
    }

    // Compute kappa = rowSums(aggte_weight)
    aggte_kappa = J(num_zeval, 1, 0)
    for (k = 1; k <= num_gte; k++) {
        aggte_kappa = aggte_kappa + aggte_weight[., k]
    }
}

// -----------------------------------------------------------------------------
// didhetero_aggte_xi()
//
// Compute influence functions xi_{i,g,t} for weight estimation uncertainty.
// The functional form depends on the aggregation type: two-term for dynamic
// and calendar, zero for group, and three-term for simple aggregation.
//
// Arguments:
//   G_g_sub      - n x num_gte group indicator submatrix
//   mu_G_g_sub   - R x num_gte conditional group density submatrix
//   aggte_weight - R x num_gte weight matrix
//   aggte_kappa  - R x 1 normalization constant
//   type         - aggregation type
//   n            - sample size
//   num_zeval    - number of evaluation points R
//   num_gte      - number of (g,t) pairs
//
// Outputs (passed by reference):
//   xi_g_t - pointer vector of length num_gte, each pointing to an n x R matrix
// -----------------------------------------------------------------------------
void didhetero_aggte_xi(
    real matrix G_g_sub,
    real matrix mu_G_g_sub,
    real matrix aggte_weight,
    real colvector aggte_kappa,
    string scalar type,
    real scalar n,
    real scalar num_zeval,
    real scalar num_gte,
    pointer(real matrix) colvector xi_g_t)
{
    real colvector mu_rowsum, term1, term2, term3, term3_inner
    real scalar i, k, id_gte, sum_G_i

    // Initialize output: pointer vector of length num_gte
    xi_g_t = J(num_gte, 1, NULL)
    for (k = 1; k <= num_gte; k++) {
        xi_g_t[k] = &(J(n, num_zeval, 0))
    }

    // Precompute mu_rowsum = sum(mu_G_g_sub, cols)
    mu_rowsum = J(num_zeval, 1, 0)
    for (k = 1; k <= num_gte; k++) {
        mu_rowsum = mu_rowsum + mu_G_g_sub[., k]
    }

    // Group type: xi = 0 (weights are deterministic)
    if (type == "group") {
        return
    }

    // Input validation
    if (type != "dynamic" & type != "calendar" & type != "simple") {
        _error(3498, "invalid aggregation type for xi: " + type)
    }

    // Dynamic / calendar type: two-term formula
    if (type == "dynamic" | type == "calendar") {

        for (i = 1; i <= n; i++) {

            // Precompute sum_G_i = sum(G_g[i, ])
            sum_G_i = 0
            for (k = 1; k <= num_gte; k++) {
                sum_G_i = sum_G_i + G_g_sub[i, k]
            }

            for (id_gte = 1; id_gte <= num_gte; id_gte++) {

                // Term 1: G_{i,g} / Sigma_mu
                term1 = J(num_zeval, 1, G_g_sub[i, id_gte]) :/ mu_rowsum

                // Term 2: mu_{G_g} / Sigma_mu^2 * sum(G_{i,g'})
                term2 = mu_G_g_sub[., id_gte] :/ (mu_rowsum :^ 2) :* sum_G_i

                // xi_{i,g,t}(z_r) = Term1 - Term2
                (*xi_g_t[id_gte])[i, .] = (term1 - term2)'
            }
        }
    }

    // Simple type: three-term formula
    else if (type == "simple") {

        for (i = 1; i <= n; i++) {

            // Precompute sum_G_i = sum(G_g[i, ])
            sum_G_i = 0
            for (k = 1; k <= num_gte; k++) {
                sum_G_i = sum_G_i + G_g_sub[i, k]
            }

            // Precompute term3_inner:
            // sum_k { G_{i,g_k} - mu_{G_{g_k}} / Sigma_mu * sum_G_i }
            // This is an R x 1 vector
            term3_inner = J(num_zeval, 1, 0)
            for (k = 1; k <= num_gte; k++) {
                term3_inner = term3_inner + (J(num_zeval, 1, G_g_sub[i, k]) - mu_G_g_sub[., k] :/ mu_rowsum :* sum_G_i)
            }

            for (id_gte = 1; id_gte <= num_gte; id_gte++) {

                // Term 1: G_{i,g} / (kappa * Sigma_mu)
                term1 = J(num_zeval, 1, G_g_sub[i, id_gte]) :/ (aggte_kappa :* mu_rowsum)

                // Term 2: mu_{G_g} / (kappa * Sigma_mu^2) * sum_G_i
                term2 = mu_G_g_sub[., id_gte] :/ (aggte_kappa :* (mu_rowsum :^ 2)) :* sum_G_i

                // Term 3: mu_{G_g} / (kappa * Sigma_mu)^2 * term3_inner
                term3 = mu_G_g_sub[., id_gte] :/ ((aggte_kappa :* mu_rowsum) :^ 2) :* term3_inner

                // xi^{OW}_{i,g,t}(z_r) = Term1 - Term2 - Term3
                (*xi_g_t[id_gte])[i, .] = (term1 - term2 - term3)'
            }
        }
    }
}

// -----------------------------------------------------------------------------
// didhetero_aggte_J()
//
// Compute the aggregated influence function J_i and point estimate by combining
// CATT influence functions B_{i,g,t} with weight influence functions xi_{i,g,t}.
//
// Arguments:
//   B_g_t_sub    - pointer vector of length num_gte, each pointing to n x R
//   xi_g_t       - pointer vector of length num_gte, each pointing to n x R
//   aggte_weight - R x num_gte weight matrix
//   catt         - R x num_gte CATT point estimate matrix
//   n            - sample size
//   num_zeval    - number of evaluation points R
//   num_gte      - number of (g,t) pairs
//
// Outputs (passed by reference):
//   J          - n x R aggregated influence function matrix
//   aggte_est  - R x 1 aggregated point estimate vector
// -----------------------------------------------------------------------------
void didhetero_aggte_J(
    pointer(real matrix) colvector B_g_t_sub,
    pointer(real matrix) colvector xi_g_t,
    real matrix aggte_weight,
    real matrix catt,
    real scalar n,
    real scalar num_zeval,
    real scalar num_gte,
    real matrix J,
    real colvector aggte_est)
{
    real colvector J_temp_row
    real scalar i, k, id_gte

    // Initialize outputs
    J = J(n, num_zeval, 0)
    aggte_est = J(num_zeval, 1, 0)

    // Compute aggregated point estimate
    for (k = 1; k <= num_gte; k++) {
        aggte_est = aggte_est + catt[., k] :* aggte_weight[., k]
    }

    // Compute J_i by summing over all (g,t) pairs
    for (id_gte = 1; id_gte <= num_gte; id_gte++) {

        for (i = 1; i <= n; i++) {

            // Contribution for (g,t) pair id_gte
            J_temp_row = aggte_weight[., id_gte] :* (*B_g_t_sub[id_gte])[i, .]' + catt[., id_gte] :* (*xi_g_t[id_gte])[i, .]'

            // Accumulate into J
            J[i, .] = J[i, .] + J_temp_row'
        }
    }
}

// -----------------------------------------------------------------------------
// _aggte_quantile_type7()
//
// Compute the type-7 sample quantile using linear interpolation.
// -----------------------------------------------------------------------------
real scalar _aggte_quantile_type7(real colvector x_sorted, real scalar prob)
{
    real scalar n, idx_r, j_r, h_r

    n = rows(x_sorted)
    if (n <= 0) return(.)
    if (n == 1) return(x_sorted[1])

    idx_r = 1 + (n - 1) * prob
    j_r = floor(idx_r)
    h_r = idx_r - j_r

    if (j_r >= n) return(x_sorted[n])
    return((1 - h_r) * x_sorted[j_r] + h_r * x_sorted[j_r + 1])
}

// -----------------------------------------------------------------------------
// _aggte_lprobust_bw()
//
// Bandwidth selection helper for aggregation standard error computation.
// Implements the lprobust bandwidth algorithm for local polynomial regression.
// -----------------------------------------------------------------------------
real scalar _aggte_lprobust_bw(
    real colvector Y,
    real colvector X,
    real scalar c_eval,
    real scalar o,
    real scalar nu,
    real scalar o_B,
    real scalar h_V,
    real scalar h_B1,
    real scalar h_B2,
    real scalar scale,
    real scalar nnmatch,
    string scalar kernel,
    real colvector dups,
    real colvector dupsid,
    real scalar V_out,
    real scalar B1_out,
    real scalar B2_out,
    real scalar R_out)
{
    real scalar N, n_V, n_B, j
    real scalar V_V, BConst1, BConst2, BWreg, V_B
    real scalar V_final, r_exp, rB, rV, bw
    real colvector w, ind_V, eY, eX, eW
    real colvector dups_V, dupsid_V, res_V
    real colvector dups_B, dupsid_B, res_B
    real matrix R_V, invG_V, RW
    real colvector Hp, v1, v2
    real colvector w_B, ind_B, eY_B, eX_B, eW_B
    real matrix R_B1, invG_B1
    real colvector beta_B1
    real matrix R_B2, invG_B2
    real colvector beta_B2
    real matrix vce_mat
    real colvector u_V

    N = rows(X)

    u_V = (X :- c_eval) / h_V
    w = didhetero_kernel_eval(u_V, kernel) / h_V
    ind_V = (w :> 0)
    n_V = sum(ind_V)
    if (n_V < o + 2) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }

    eY = select(Y, ind_V)
    eX = select(X, ind_V)
    eW = select(w, ind_V)

    R_V = J(n_V, o + 1, .)
    for (j = 1; j <= o + 1; j++) {
        R_V[., j] = (eX :- c_eval) :^ (j - 1)
    }

    invG_V = cholinv(cross(R_V :* sqrt(eW), R_V :* sqrt(eW)))
    if (hasmissing(invG_V)) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }

    dups_V = select(dups, ind_V)
    dupsid_V = select(dupsid, ind_V)
    res_V = _didhetero_lprobust_res(eX, eY, nnmatch, dups_V, dupsid_V)

    RW = R_V :* eW
    vce_mat = _didhetero_lprobust_vce(RW, res_V)
    V_V = (invG_V * vce_mat * invG_V)[nu + 1, nu + 1]

    Hp = J(o + 1, 1, .)
    for (j = 1; j <= o + 1; j++) {
        Hp[j] = h_V^(j - 1)
    }

    v1 = cross(R_V, eW :* (((eX :- c_eval) / h_V) :^ (o + 1)))
    v2 = cross(R_V, eW :* (((eX :- c_eval) / h_V) :^ (o + 2)))

    BConst1 = (Hp :* (invG_V * v1))[nu + 1]
    BConst2 = (Hp :* (invG_V * v2))[nu + 1]

    w_B = didhetero_kernel_eval((X :- c_eval) / h_B1, kernel)
    ind_B = (w_B :> 0)
    n_B = sum(ind_B)
    if (n_B < o_B + 2) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }

    eY_B = select(Y, ind_B)
    eX_B = select(X, ind_B)
    eW_B = select(w_B, ind_B)

    R_B1 = J(n_B, o_B + 1, .)
    for (j = 1; j <= o_B + 1; j++) {
        R_B1[., j] = (eX_B :- c_eval) :^ (j - 1)
    }
    invG_B1 = cholinv(cross(R_B1 :* sqrt(eW_B), R_B1 :* sqrt(eW_B)))
    if (hasmissing(invG_B1)) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }
    beta_B1 = invG_B1 * cross(R_B1, eW_B :* eY_B)

    BWreg = 0
    if (scale > 0) {
        dups_B = select(dups, ind_B)
        dupsid_B = select(dupsid, ind_B)
        res_B = _didhetero_lprobust_res(eX_B, eY_B, nnmatch, dups_B, dupsid_B)
        vce_mat = _didhetero_lprobust_vce(R_B1 :* eW_B, res_B)
        V_B = (invG_B1 * vce_mat * invG_B1)[o + 2, o + 2]
        BWreg = 3 * BConst1^2 * V_B
    }

    w_B = didhetero_kernel_eval((X :- c_eval) / h_B2, kernel)
    ind_B = (w_B :> 0)
    n_B = sum(ind_B)
    if (n_B < o_B + 3) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }

    eY_B = select(Y, ind_B)
    eX_B = select(X, ind_B)
    eW_B = select(w_B, ind_B)

    R_B2 = J(n_B, o_B + 2, .)
    for (j = 1; j <= o_B + 2; j++) {
        R_B2[., j] = (eX_B :- c_eval) :^ (j - 1)
    }
    invG_B2 = cholinv(cross(R_B2 :* sqrt(eW_B), R_B2 :* sqrt(eW_B)))
    if (hasmissing(invG_B2)) {
        V_out = .
        B1_out = .
        B2_out = .
        R_out = 0
        return(.)
    }
    beta_B2 = invG_B2 * cross(R_B2, eW_B :* eY_B)

    B1_out = BConst1 * beta_B1[o + 2]
    B2_out = BConst2 * beta_B2[o + 3]
    V_final = N * h_V^(2 * nu + 1) * V_V
    V_out = V_final
    R_out = BWreg

    r_exp = 1 / (2 * o + 3)
    rB = 2 * (o + 1 - nu)
    rV = 2 * nu + 1

    if (abs(B1_out^2 + scale * R_out) < 1e-30) return(.)

    bw = ((rV * V_final) / (N * rB * (B1_out^2 + scale * R_out)))^r_exp
    return(bw)
}

// -----------------------------------------------------------------------------
// _aggte_lpbw_mse_details_odd()
//
// Compute pointwise MSE-DPI bandwidth details for the odd polynomial order case.
// Used in aggregation standard error estimation with p=1 and deriv=0.
// -----------------------------------------------------------------------------
void _aggte_lpbw_mse_details_odd(
    real colvector Y,
    real colvector X,
    real scalar eval_pt,
    real scalar p,
    real scalar deriv,
    string scalar kernel,
    real scalar h_out,
    real scalar Vh_out,
    real scalar Bh_out)
{
    real scalar N, q, C_c, x_sd, x_iq, x_min, x_max, range_X
    real scalar c_bw, bw_max_r, bw_min
    real scalar bw_mp2, bw_mp3, b_mse_dpi, h_val
    real scalar V_tmp, B1_tmp, B2_tmp, R_tmp
    real scalar j, rV_h, rB_h
    real colvector sort_idx, X_s, Y_s, dups, dupsid, abs_dist, X_sorted

    h_out = .
    Vh_out = .
    Bh_out = .

    N = rows(X)
    if (N <= 1) return

    q = p + 1

    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else return

    x_sd = sqrt(variance(X))
    X_sorted = sort(X, 1)
    x_iq = _aggte_quantile_type7(X_sorted, 0.75) - _aggte_quantile_type7(X_sorted, 0.25)

    x_min = min(X)
    x_max = max(X)
    range_X = x_max - x_min
    if (range_X <= 0) return

    sort_idx = order(X, 1)
    X_s = X[sort_idx]
    Y_s = Y[sort_idx]

    dups = J(N, 1, 0)
    dupsid = J(N, 1, 0)
    for (j = 1; j <= N; j++) {
        dups[j] = sum(X_s :== X_s[j])
    }
    j = 1
    while (j <= N) {
        real scalar k_dup
        for (k_dup = 0; k_dup < dups[j]; k_dup++) {
            dupsid[j + k_dup] = k_dup + 1
        }
        j = j + dups[j]
    }

    bw_max_r = max((abs(eval_pt - x_min), abs(eval_pt - x_max)))
    c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
    c_bw = min((c_bw, bw_max_r))

    abs_dist = sort(abs(X_s :- eval_pt), 1)
    if (21 <= N) bw_min = abs_dist[21]
    else bw_min = abs_dist[N]
    c_bw = max((c_bw, bw_min))

    bw_mp2 = _aggte_lprobust_bw(Y_s, X_s, eval_pt,
                q + 1, q + 1, q + 2,
                c_bw, range_X, range_X, 0,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (bw_mp2 == .) bw_mp2 = range_X
    bw_mp2 = min((bw_mp2, bw_max_r))
    bw_mp2 = max((bw_mp2, bw_min))

    bw_mp3 = _aggte_lprobust_bw(Y_s, X_s, eval_pt,
                q + 2, q + 2, q + 3,
                c_bw, range_X, range_X, 0,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (bw_mp3 == .) bw_mp3 = range_X
    bw_mp3 = min((bw_mp3, bw_max_r))
    bw_mp3 = max((bw_mp3, bw_min))

    b_mse_dpi = _aggte_lprobust_bw(Y_s, X_s, eval_pt,
                q, p + 1, q + 1,
                c_bw, bw_mp2, bw_mp3, 1,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (b_mse_dpi == .) b_mse_dpi = range_X
    b_mse_dpi = min((b_mse_dpi, bw_max_r))
    b_mse_dpi = max((b_mse_dpi, bw_min))

    h_val = _aggte_lprobust_bw(Y_s, X_s, eval_pt,
                p, deriv, q,
                c_bw, b_mse_dpi, bw_mp2, 1,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (h_val == .) h_val = range_X
    h_val = min((h_val, bw_max_r))
    h_val = max((h_val, bw_min))

    if (V_tmp >= . | B1_tmp >= .) return

    rV_h = 2 * deriv + 1
    rB_h = 2 * (p + 1 - deriv)

    h_out = h_val
    Vh_out = rV_h * V_tmp
    Bh_out = rB_h * B1_tmp^2
}

// -----------------------------------------------------------------------------
// _aggte_lpbwselect_mse_odd()
//
// Compute pointwise MSE-DPI bandwidths for the odd polynomial order case.
// -----------------------------------------------------------------------------
real colvector _aggte_lpbwselect_mse_odd(
    real colvector Y,
    real colvector X,
    real colvector eval,
    real scalar p,
    real scalar deriv,
    string scalar kernel)
{
    real scalar R_eval, r_idx, h_tmp, Vh_tmp, Bh_tmp
    real colvector h_mse

    R_eval = rows(eval)
    h_mse = J(R_eval, 1, .)

    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        _aggte_lpbw_mse_details_odd(Y, X, eval[r_idx], p, deriv,
                                    kernel, h_tmp, Vh_tmp, Bh_tmp)
        h_mse[r_idx] = h_tmp
    }

    return(h_mse)
}

// -----------------------------------------------------------------------------
// _aggte_lpbwselect_imse_odd()
//
// Compute the common IMSE-DPI bandwidth for the odd polynomial order case.
// -----------------------------------------------------------------------------
real scalar _aggte_lpbwselect_imse_odd(
    real colvector Y,
    real colvector X,
    real colvector eval,
    real scalar p,
    real scalar deriv,
    string scalar kernel)
{
    real scalar N, R_eval, r_idx, h_tmp, Vh_tmp, Bh_tmp
    real scalar mean_Vh, mean_Bh, n_valid_V, n_valid_B, range_X, r_exp
    real scalar x_min, x_max
    real colvector imse_grid

    N = rows(X)
    x_min = min(X)
    x_max = max(X)
    range_X = x_max - x_min
    if (N <= 1 | range_X <= 0) return(.)

    // IMSE-DPI averages pilot quantities over an internal 30-point equally
    // spaced grid over the support of X, ignoring the user-supplied eval.
    imse_grid = x_min :+ ((0::29) :/ 29) :* range_X
    R_eval = rows(imse_grid)

    mean_Vh = 0
    mean_Bh = 0
    n_valid_V = 0
    n_valid_B = 0
    r_exp = 1 / (2 * p + 3)

    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        _aggte_lpbw_mse_details_odd(Y, X, imse_grid[r_idx], p, deriv,
                                    kernel, h_tmp, Vh_tmp, Bh_tmp)
        if (Vh_tmp < .) {
            mean_Vh = mean_Vh + Vh_tmp
            n_valid_V = n_valid_V + 1
        }
        if (Bh_tmp < .) {
            mean_Bh = mean_Bh + Bh_tmp
            n_valid_B = n_valid_B + 1
        }
    }

    if (n_valid_V == 0 | n_valid_B == 0) return(.)

    mean_Vh = mean_Vh / n_valid_V
    mean_Bh = mean_Bh / n_valid_B

    if (mean_Bh < 1e-20) return(range_X / 2)
    if (mean_Vh <= 0) return(.)

    return((mean_Vh / (N * mean_Bh))^r_exp)
}

// -----------------------------------------------------------------------------
// didhetero_aggte_bw_pass1()
//
// Compute the IMSE-optimal bandwidth for the aggregated parameter (Pass 1).
// Applies the IMSE-DPI pipeline to the aggregated influence function J_i,
// supporting IMSE1, IMSE2, and US1 bandwidth selection methods.
//
// Arguments:
//   J          - n x R aggregated influence function matrix
//   Z          - n x 1 covariate vector
//   zeval      - R x 1 sorted evaluation points
//   kd0_Z      - R x 1 kernel density estimates at zeval
//   kd1_Z      - R x 1 kernel density derivative estimates at zeval
//   n          - sample size
//   bwselect   - bandwidth selection method: "IMSE1", "IMSE2", or "US1"
//   kernel     - kernel type ("epa" or "gau")
//   uniformall - 1 for common bandwidth
//
// Returns:
//   scalar bandwidth, or missing on failure
// -----------------------------------------------------------------------------
real scalar didhetero_aggte_bw_pass1(
    real matrix J,
    real colvector Z,
    real colvector zeval,
    real colvector kd0_Z,
    real colvector kd1_Z,
    real scalar n,
    string scalar bwselect,
    string scalar kernel,
    real scalar uniformall)
{
    real scalar R_eval, r
    real scalar I_2_K1, I_4_K1, I_6_K1, I_0_K2, I_2_K2, I_4_K2, I_6_K2
    real scalar const_V1, const_V2, const_V
    real scalar C_B_LQ, cb_num, cb_den, cv_num, cv_den
    real colvector mathcal_B, mathcal_V
    real colvector y_r
    real colvector mu_J_2_bw_vec, mu_J_3_bw_vec, mu_J_4_bw_vec
    real scalar mu_J_2_bw, mu_J_2
    real scalar mu_J_3_bw, mu_J_3, mu_J_4_bw, mu_J_4
    real scalar int_bias, int_var, h_opt, h_max

    // Input validation
    if (bwselect != "IMSE1" & bwselect != "IMSE2" & bwselect != "US1") {
        _error(3498, "invalid bwselect for aggte Pass 1: " + bwselect)
    }

    R_eval = rows(zeval)

    // Need at least 2 evaluation points for trapezoidal integration
    if (R_eval < 2) return(.)

    // Get kernel constants
    _didhetero_kernel_constants(kernel, I_2_K1, I_4_K1, I_6_K1,
                                I_0_K2, I_2_K2, I_4_K2, I_6_K2)

    // Compute derived constants
    const_V1 = I_0_K2

    if (bwselect == "IMSE2") {
        // C_{B,LQ} = (I_4^2 - I_2*I_6) / (I_4 - I_2^2)
        cb_num = I_4_K1^2 - I_2_K1 * I_6_K1
        cb_den = I_4_K1 - I_2_K1^2
        if (abs(cb_den) < 1e-30) return(.)
        C_B_LQ = cb_num / cb_den

        // C_{V2} = (I_4^2*I_0_K2 - 2*I_2*I_4*I_2_K2 + I_2^2*I_4_K2) / (I_4 - I_2^2)^2
        cv_num = I_4_K1^2 * I_0_K2 - 2 * I_2_K1 * I_4_K1 * I_2_K2 + I_2_K1^2 * I_4_K2
        cv_den = (I_4_K1 - I_2_K1^2)^2
        if (abs(cv_den) < 1e-30) return(.)
        const_V2 = cv_num / cv_den
    }

    // Select const_V by bwselect
    if (bwselect == "IMSE1" | bwselect == "US1") {
        const_V = const_V1
    }
    else {
        const_V = const_V2
    }

    // Initialize bias vector
    mathcal_B = J(R_eval, 1, .)

    // Bias estimation loop
    for (r = 1; r <= R_eval; r++) {

        y_r = J[., r]

        // Skip if density is non-positive
        if (kd0_Z[r] >= . | kd0_Z[r] <= 0) continue

        if (bwselect == "IMSE1" | bwselect == "US1") {

            // IMSE1/US1 bias: 2nd derivative (p=3, deriv=2)

            // MSE-DPI bandwidth for mu_J^{(2)}
            mu_J_2_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                                zeval[r..r], 3, 2, kernel)
            mu_J_2_bw = mu_J_2_bw_vec[1]

            if (mu_J_2_bw >= . | mu_J_2_bw <= 0) continue

            // Estimate 2nd derivative via LPR (p=3, deriv=2)
            {
                real colvector mu_J_2_vec
                mu_J_2_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                                 3, 2, kernel, mu_J_2_bw)
                mu_J_2 = mu_J_2_vec[1]
            }

            if (mu_J_2 >= .) continue

            // IMSE1 bias for LLR: mathcal_B[r] = mu^{(2)} * I_{2,K} / 2
            mathcal_B[r] = mu_J_2 * I_2_K1 / 2
        }
        else {

            // IMSE2 bias: 3rd and 4th derivatives

            // 3rd derivative: MSE-DPI bw (p=4, deriv=3) -> LPR
            mu_J_3_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                                zeval[r..r], 4, 3, kernel)
            mu_J_3_bw = mu_J_3_bw_vec[1]

            if (mu_J_3_bw >= . | mu_J_3_bw <= 0) continue

            {
                real colvector mu_J_3_vec
                mu_J_3_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                                 4, 3, kernel, mu_J_3_bw)
                mu_J_3 = mu_J_3_vec[1]
            }

            if (mu_J_3 >= .) continue

            // 4th derivative: MSE-DPI bw (p=5, deriv=4) -> LPR
            mu_J_4_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                                zeval[r..r], 5, 4, kernel)
            mu_J_4_bw = mu_J_4_bw_vec[1]

            if (mu_J_4_bw >= . | mu_J_4_bw <= 0) continue

            {
                real colvector mu_J_4_vec
                mu_J_4_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                                 5, 4, kernel, mu_J_4_bw)
                mu_J_4 = mu_J_4_vec[1]
            }

            if (mu_J_4 >= .) continue

            // IMSE2 bias formula
            mathcal_B[r] = (1 / (24 * kd0_Z[r])) * (2 * mu_J_3 * kd1_Z[r] + mu_J_4 * kd0_Z[r]) * C_B_LQ
        }
    }

    // Variance estimation
    _didhetero_bwselect_var_est(J, Z, zeval, kd0_Z, const_V,
                                           kernel, mathcal_V)

    // Trapezoidal integration
    int_bias = _didhetero_trapz(zeval, mathcal_B :^ 2)
    int_var  = _didhetero_trapz(zeval, mathcal_V)

    // Edge cases
    h_max = (max(Z) - min(Z)) / 2

    if (int_bias >= . | int_bias < 1e-20) return(h_max)
    if (int_var >= . | int_var <= 0) return(.)

    // Bandwidth formula
    if (bwselect == "IMSE1") {
        h_opt = (int_var / (4 * int_bias))^(1/5) * n^(-1/5)
    }
    else if (bwselect == "IMSE2") {
        h_opt = (int_var / (8 * int_bias))^(1/9) * n^(-1/9)
    }
    else if (bwselect == "US1") {
        h_opt = (int_var / (4 * int_bias))^(1/5) * n^(-2/7)
    }

    if (h_opt >= . | h_opt <= 0) return(.)

    return(h_opt)
}

// -----------------------------------------------------------------------------
// didhetero_aggte_se()
//
// Compute standard errors and analytical uniform confidence bands for
// aggregated treatment effect estimates. The variance estimation pipeline
// consists of: (1) bandwidth selection for the residual mean, (2) local
// polynomial regression of the influence function on the covariate,
// (3) residual computation and variance estimation, and (4) standard
// error calculation. Analytical uniform confidence bands are constructed
// using the critical value formula for suprema of Gaussian processes.
//
// Arguments:
//   J          - n x num_zeval aggregated influence function matrix
//   aggte_est  - num_zeval x 1 aggregated point estimate vector
//   Z          - n x 1 covariate vector
//   zeval      - num_zeval x 1 evaluation points
//   Z_supp     - support grid for IMSE-DPI bandwidth selection
//   kd0_Z      - num_zeval x 1 kernel density estimates
//   bw         - scalar bandwidth
//   n          - sample size
//   porder     - polynomial order (1 or 2)
//   kernel     - kernel type ("epa" or "gau")
//   alp        - significance level
//   num_zeval  - number of evaluation points
//
// Outputs (passed by reference):
//   se         - num_zeval x 1 standard errors
//   mathcal_V  - num_zeval x 1 variance estimates
//   U_hat      - n x num_zeval residual matrix
//   ci1_lower  - num_zeval x 1 analytical CI lower bounds
//   ci1_upper  - num_zeval x 1 analytical CI upper bounds
// -----------------------------------------------------------------------------
void didhetero_aggte_se(
    real matrix J,
    real colvector aggte_est,
    real colvector Z,
    real colvector zeval,
    real colvector Z_supp,
    real colvector kd0_Z,
    real scalar bw,
    real scalar n,
    real scalar porder,
    string scalar kernel,
    real scalar alp,
    real scalar num_zeval,
    real colvector se,
    real colvector mathcal_V,
    real matrix U_hat,
    real colvector ci1_lower,
    real colvector ci1_upper)
{
    real scalar const_V, lambda
    real scalar r

    // Initialize outputs with missing values
    se        = J(num_zeval, 1, .)
    mathcal_V = J(num_zeval, 1, .)
    U_hat     = J(n, num_zeval, .)
    ci1_lower = J(num_zeval, 1, .)
    ci1_upper = J(num_zeval, 1, .)

    // const_V selection by polynomial order
    if (porder == 1) {
        if (kernel == "epa") {
            const_V = 3/5                          // 0.6
        }
        else {
            const_V = 1 / (2 * sqrt(c("pi")))      // ~0.2820948
        }
    }
    else if (porder == 2) {
        if (kernel == "epa") {
            const_V = 5/4                           // 1.25
        }
        else {
            const_V = 27 / (32 * sqrt(c("pi")))    // ~0.4760350
        }
    }
    else {
        _error(3498, "invalid porder for aggte SE: must be 1 or 2")
    }

    // Lambda selection for analytical uniform confidence bands
    if (kernel == "epa") {
        lambda = 2.5
    }
    else if (kernel == "gau") {
        lambda = 0.5
    }
    else {
        _error(3498, "invalid kernel for aggte SE: " + kernel)
    }

    // Guard: invalid inputs
    if (n <= 0 | bw <= 0) {
        return
    }

    // Main loop over evaluation points
    {
        real scalar mu_J_bw, sigma2_bw_scalar, sigma2, V_hat_r
        real colvector J_r, mu_J_hat, sigma2_bw_vec

        for (r = 1; r <= num_zeval; r++) {

            // Extract J_i(z_r) column
            J_r = J[., r]

            // Step 1: IMSE-DPI bandwidth for residual mean (p=1)
            mu_J_bw = _aggte_lpbwselect_imse_odd(J_r, Z, Z_supp,
                          1, 0, kernel)

            if (mu_J_bw >= . | mu_J_bw <= 0) {
                se[r] = .
                continue
            }

            // Step 2: LPR estimate at all observations (p=1)
            mu_J_hat = didhetero_lpr(J_r, Z, Z, 1, 0, kernel, mu_J_bw)

            // Step 3a: Residuals
            U_hat[., r] = J_r - mu_J_hat

            // Step 3b: MSE-DPI bandwidth for squared residuals (p=1)
            sigma2_bw_vec = _aggte_lpbwselect_mse_odd(U_hat[., r] :^ 2, Z,
                                zeval[r], 1, 0, kernel)
            sigma2_bw_scalar = sigma2_bw_vec[1]

            if (sigma2_bw_scalar >= . | sigma2_bw_scalar <= 0) {
                se[r] = .
                continue
            }

            // Step 3c: LPR estimate of sigma^2 (p=1)
            real scalar sigma2_raw_aggte
            sigma2_raw_aggte = didhetero_lpr(U_hat[., r] :^ 2, Z, zeval[r],
                                   1, 0, kernel, sigma2_bw_scalar)

            // Truncate negative variance to zero
            sigma2 = sigma2_raw_aggte
            if (sigma2[1] < 0) sigma2 = 0

            // Step 4: Variance and standard error
            if (kd0_Z[r] > 0 & kd0_Z[r] < .) {
                V_hat_r = const_V * sigma2[1] / kd0_Z[r]
            }
            else {
                V_hat_r = .
            }

            mathcal_V[r] = V_hat_r

            if (sigma2_raw_aggte[1] < 0) {
                // Negative variance indicates insufficient local data
                se[r] = .
            }
            else if (V_hat_r != . & V_hat_r >= 0) {
                se[r] = sqrt(V_hat_r / (n * bw))
            }
            else {
                se[r] = .
            }
        }
    }

    // Analytical critical value and uniform confidence bands
    {
        real scalar c_hat, n_missing

        c_hat = didhetero_analytical_crit(zeval, bw, lambda, alp)

        if (c_hat < .) {
            ci1_lower = aggte_est - c_hat :* se
            ci1_upper = aggte_est + c_hat :* se
        }
        else {
            ci1_lower = J(num_zeval, 1, .)
            ci1_upper = J(num_zeval, 1, .)
            printf("{txt}Warning: analytical UCB critical value cannot be computed; leaving aggte analytical CI missing\n")
        }

        // Diagnostic: count missing SEs
        n_missing = sum(se :>= .)
        if (n_missing > 0) {
            printf("Warning: aggte SE=. at %g of %g evaluation point(s): insufficient local data\n",
                n_missing, num_zeval)
            printf("  (boundary z-values far from treated/control units; CI reported as missing)\n")
        }
    }
}

// =============================================================================
// didhetero_aggte_bootstrap()
//
// Aggregation-level multiplier bootstrap for uniform confidence bands (Phase A).
// For each evaluation point, computes kernel-weighted bootstrap t-statistics
// across all z evaluation points and stores the supremum. Phase B (critical
// value computation and confidence interval construction) is handled separately
// after all evaluation points have been processed.
// =============================================================================
void didhetero_aggte_bootstrap(
    real colvector aggte_est,
    real colvector se,
    real matrix U_hat,
    real colvector Z,
    real colvector zeval,
    real colvector kd0_Z,
    real scalar bw,
    real scalar n,
    real scalar porder,
    string scalar kernel,
    real scalar biters,
    real scalar alp,
    real scalar uniformall,
    real scalar num_zeval,
    real scalar num_eval,
    real scalar id_eval,
    real matrix mb_weight,
    real matrix mb_sup_t)
{
    // --- Local declarations ---
    real scalar r, e, b
    real scalar I_2_K1, I_4_K1, I_6_K1, I_0_K2, I_2_K2, I_4_K2, I_6_K2
    real scalar scale, c_global, c_per_eval
    real colvector u_r, kv_r, Psi_r, perturb
    real colvector mb_est, mb_t_r
    real matrix mb_t

    // --- Input validation ---
    if (biters < 1) {
        _error(3498, "biters must be >= 1")
    }

    // --- Get kernel constants ---
    _didhetero_kernel_constants(kernel, I_2_K1, I_4_K1, I_6_K1,
                                I_0_K2, I_2_K2, I_4_K2, I_6_K2)

    // Accumulate bootstrap t-statistics across zeval points
    mb_t = J(biters, num_zeval, .)

    for (r = 1; r <= num_zeval; r++) {

        // Scaled distance u_{ih} = (Z_i - zeval[r]) / bw
        u_r = (Z :- zeval[r]) / bw

        // Kernel values K(u_{ih})
        kv_r = didhetero_kernel_eval(u_r, kernel)

        // Equivalent kernel Psi
        if (porder == 1) {
            Psi_r = J(n, 1, 1)
        }
        else if (porder == 2) {
            Psi_r = (I_4_K1 :- u_r:^2 :* I_2_K1) :/ (I_4_K1 - I_2_K1^2)
        }
        else {
            _error(3498, "porder must be 1 or 2")
        }

        // Perturbation vector = Psi_r .* U_hat[.,r] .* kv_r
        perturb = Psi_r :* U_hat[., r] :* kv_r

        // Scale factor = 1 / (kd0_Z[r] * n * bw)
        if (kd0_Z[r] <= 0 | kd0_Z[r] >= .) {
            mb_t[., r] = J(biters, 1, .)
            continue
        }
        scale = 1 / (kd0_Z[r] * n * bw)

        // Vectorized bootstrap estimates (all B iterations at once)
        mb_est = J(biters, 1, aggte_est[r]) :+ scale :* ((mb_weight :- 1) * perturb)

        // Bootstrap t-statistics
        if (se[r] > 0 & se[r] < .) {
            mb_t_r = abs(mb_est :- aggte_est[r]) :/ se[r]
        }
        else {
            mb_t_r = J(biters, 1, .)
        }

        // Store t-statistics for this zeval point
        mb_t[., r] = mb_t_r
    }

    // Sup-t for this eval: row-wise max over all zeval points
    for (b = 1; b <= biters; b++) {
        mb_sup_t[b, id_eval] = didhetero_max_nonmissing(mb_t[b, .]')
    }
}


// =============================================================================
// _didhetero_unflatten_B_g_t()
//
// Convert a flattened B_g_t matrix (n x R*K) into a pointer column vector of
// length K, where each pointer references an n x R matrix. The flattened format
// stores columns as [B_1, B_2, ..., B_K] where each B_k contains R columns
// corresponding to the k-th (g,t) pair.
//
// Arguments:
//   B_flat      - n x (R*K) flattened influence function matrix
//   n           - sample size
//   num_zeval   - R, number of z evaluation points
//   num_gteval  - K, number of (g,t) pairs
//
// Returns:
//   pointer column vector of length K
// =============================================================================
pointer(real matrix) colvector _didhetero_unflatten_B_g_t(
    real matrix B_flat,
    real scalar n,
    real scalar num_zeval,
    real scalar num_gteval)
{
    pointer(real matrix) colvector B_g_t
    real scalar k, col_lo, col_hi

    B_g_t = J(num_gteval, 1, NULL)

    for (k = 1; k <= num_gteval; k++) {
        col_lo = (k - 1) * num_zeval + 1
        col_hi = k * num_zeval
        // Force independent copy to avoid pointer aliasing
        B_g_t[k] = &(1 * B_flat[., col_lo..col_hi])
    }

    return(B_g_t)
}


// =============================================================================
// didhetero_aggte_pass2()
//
// Pass 2 re-estimation for a single aggregation evaluation point. Re-estimates
// CATT parameters using the aggregation-optimal bandwidth, then recomputes
// weights, influence functions, and aggregated estimates. When num_gte equals 1,
// the CATT estimate is returned directly as the aggregated result.
//
// Arguments:
//   data        - DidHeteroData struct
//   gps_mat     - GPS result matrix from Stage 1
//   or_mat      - OR result matrix from Stage 1
//   gteeval_sub - M x 3 matrix of (g, t, eval) triples
//   h_agg       - aggregation-optimal bandwidth
//   type        - aggregation type
//   gbar        - upper bound for control group
//
// Outputs (passed by reference):
//   aggte_est      - R x 1 aggregated point estimate
//   J              - n x R aggregated influence function
//   aggte_weight   - R x num_gte weight matrix
//   aggte_kappa    - R x 1 normalization constant
//   catt_est_new   - R x num_gte CATT point estimates
//   B_g_t_cv       - pointer column vector of length num_gte
//   G_g_new        - n x num_gte group indicator matrix
//   mu_G_g_new     - R x num_gte conditional group density
// =============================================================================
void didhetero_aggte_pass2(
    struct DidHeteroData scalar data,
    real matrix gps_mat,
    real matrix or_mat,
    real matrix gteeval_sub,
    real scalar h_agg,
    string scalar type,
    real scalar gbar,
    real scalar porder,
    string scalar kernel,
    real colvector aggte_est,
    real matrix J,
    real matrix aggte_weight,
    real colvector aggte_kappa,
    real matrix catt_est_new,
    pointer(real matrix) colvector B_g_t_cv,
    real matrix G_g_new,
    real matrix mu_G_g_new)
{
    // --- Local declarations ---
    struct DidHeteroCattResult scalar catt_result
    real scalar num_gte, n, num_zeval, k
    real matrix gteval_gt
    real colvector bw_scalar
    pointer(real matrix) colvector xi_g_t

    // Input validation
    if (h_agg <= 0 | h_agg >= .) {
        _error(3498, "didhetero_aggte_pass2: h_agg must be > 0 and non-missing")
    }
    if (rows(gteeval_sub) == 0) {
        _error(3498, "didhetero_aggte_pass2: gteeval_sub is empty")
    }
    if (cols(gteeval_sub) != 3) {
        _error(3498, "didhetero_aggte_pass2: gteeval_sub must have 3 columns (g, t, eval)")
    }

    // Variable initialization
    num_gte   = rows(gteeval_sub)
    n         = data.n
    num_zeval = data.num_zeval

    // Short-circuit: when only one (g,t) pair maps to this evaluation point,
    // the aggregated parameter equals the CATT parameter
    if (num_gte == 1) {

        gteval_gt = gteeval_sub[., 1::2]
        bw_scalar = J(1, 1, h_agg)

        catt_result = didhetero_catt_core(data, gps_mat, or_mat,
                          gteval_gt, bw_scalar, "manual",
                          porder, kernel)

        catt_est_new = catt_result.catt_est
        B_g_t_cv     = J(1, 1, NULL)
        B_g_t_cv[1]  = catt_result.B_g_t[1]
        G_g_new      = catt_result.G_g
        mu_G_g_new   = catt_result.mu_G_g

        aggte_est    = catt_est_new[., 1]
        J            = J(n, num_zeval, .)
        aggte_weight = J(num_zeval, 1, 1)
        aggte_kappa  = J(num_zeval, 1, 1)

        return
    }

    // Internal CATT call with aggregation-optimal bandwidth
    gteval_gt = gteeval_sub[., 1::2]
    bw_scalar = J(1, 1, h_agg)

    // Use aggte-level porder/kernel for the re-estimation pass
    catt_result = didhetero_catt_core(data, gps_mat, or_mat,
                      gteval_gt, bw_scalar, "manual",
                      porder, kernel)

    // Extract Pass 2 results
    catt_est_new = catt_result.catt_est
    G_g_new      = catt_result.G_g
    mu_G_g_new   = catt_result.mu_G_g

    // Convert pointer rowvector to column vector
    B_g_t_cv = J(num_gte, 1, NULL)
    for (k = 1; k <= num_gte; k++) {
        B_g_t_cv[k] = catt_result.B_g_t[k]
    }

    // Dimension validation
    if (cols(catt_est_new) != num_gte) {
        _error(3498, sprintf("aggte_pass2: catt_est cols=%g, expected num_gte=%g",
                             cols(catt_est_new), num_gte))
    }
    if (cols(G_g_new) != num_gte) {
        _error(3498, sprintf("aggte_pass2: G_g cols=%g, expected num_gte=%g",
                             cols(G_g_new), num_gte))
    }
    if (cols(mu_G_g_new) != num_gte) {
        _error(3498, sprintf("aggte_pass2: mu_G_g cols=%g, expected num_gte=%g",
                             cols(mu_G_g_new), num_gte))
    }

    // Pass 2 weight recalculation
    didhetero_aggte_weights(mu_G_g_new, gteeval_sub, type, gbar,
                            num_zeval, num_gte,
                            aggte_weight, aggte_kappa)

    // Pass 2 influence function recalculation
    didhetero_aggte_xi(G_g_new, mu_G_g_new, aggte_weight, aggte_kappa,
                       type, n, num_zeval, num_gte,
                       xi_g_t)

    // Pass 2 J and aggregated point estimate recalculation
    didhetero_aggte_J(B_g_t_cv, xi_g_t, aggte_weight, catt_est_new,
                      n, num_zeval, num_gte,
                      J, aggte_est)
}


// =============================================================================
// _didhetero_aggte_build_eval()
//
// Build the default evaluation vector for aggregation based on type, or validate
// and use user-specified evaluation values. Default construction by type:
//   dynamic:  unique sorted values of (t - g)
//   group:    unique sorted values of g (post-treatment only)
//   calendar: unique sorted values of t (post-treatment only)
//   simple:   single missing value
//
// User-specified values override defaults after validation against the default set.
//
// Arguments:
//   gteval    - K x 2 matrix of (g, t) pairs
//   type      - aggregation type
//   user_eval - column vector of user-specified values (empty if not specified)
//
// Returns:
//   column vector of evaluation points
// =============================================================================
real colvector _didhetero_aggte_build_eval(
    real matrix gteval,
    string scalar type,
    real colvector user_eval)
{
    real colvector raw_vals, default_eval, sorted_vals, dyn_support
    real scalar K, i, n_raw, n_unique, j, found, is_dup
    real colvector unique_vals

    K = rows(gteval)

    // Build default eval based on type

    if (type == "simple") {
        default_eval = J(1, 1, .)
    }
    else if (type == "dynamic") {
        raw_vals = J(K, 1, .)
        dyn_support = _aggte_time_support_from_gteval(gteval)
        for (i = 1; i <= K; i++) {
            raw_vals[i] = _aggte_dynamic_event_time(
                gteval[i, 1], gteval[i, 2], dyn_support,
                "_didhetero_aggte_build_eval()")
        }
        default_eval = _didhetero_unique_sorted(raw_vals)
    }
    else if (type == "group") {
        // Group: unique sorted values of g for post-treatment pairs only
        n_raw = 0
        for (i = 1; i <= K; i++) {
            if (gteval[i, 2] >= gteval[i, 1]) {
                n_raw++
            }
        }
        if (n_raw == 0) {
            _error(3498, "aggte_gt: no post-treatment pairs found for group type")
        }
        raw_vals = J(n_raw, 1, .)
        j = 0
        for (i = 1; i <= K; i++) {
            if (gteval[i, 2] >= gteval[i, 1]) {
                j++
                raw_vals[j] = gteval[i, 1]
            }
        }
        default_eval = _didhetero_unique_sorted(raw_vals)
    }
    else if (type == "calendar") {
        // Calendar: unique sorted values of t for post-treatment periods only
        n_raw = 0
        for (i = 1; i <= K; i++) {
            if (gteval[i, 2] >= gteval[i, 1]) {
                n_raw++
            }
        }
        if (n_raw == 0) {
            _error(3498, "aggte_gt: no post-treatment periods found for calendar type")
        }
        raw_vals = J(n_raw, 1, .)
        j = 0
        for (i = 1; i <= K; i++) {
            if (gteval[i, 2] >= gteval[i, 1]) {
                j++
                raw_vals[j] = gteval[i, 2]
            }
        }
        default_eval = _didhetero_unique_sorted(raw_vals)
    }
    else {
        _error(3498, "invalid aggregation type: " + type)
    }

    // User override or return default

    if (rows(user_eval) == 0) {
        return(default_eval)
    }

    if (type == "simple") {
        _error(198, "aggte_gt: eval() is not allowed when type(simple)")
    }

    // =====================================================================
    // Step 3: Validate user-specified eval values
    // Check that each user value exists in the default set
    // =====================================================================
    for (i = 1; i <= rows(user_eval); i++) {
        found = 0
        for (j = 1; j <= rows(default_eval); j++) {
            if (user_eval[i] == default_eval[j]) {
                found = 1
                break
            }
        }
        if (!found) {
            printf("{err}aggte_gt: eval value %g is not in the valid set for type '%s'\n",
                   user_eval[i], type)
            printf("{err}  Valid values: ")
            for (j = 1; j <= rows(default_eval); j++) {
                printf("%g ", default_eval[j])
            }
            printf("\n")
            _error(198, "aggte_gt: invalid eval value specified")
        }
    }

    return(user_eval)
}


// =============================================================================
// _aggte_requires_pass1_bw()
//
// Return 1 when any requested eval point maps to multiple (g,t) pairs and
// therefore requires Pass 1 automatic aggregation bandwidth selection.
// =============================================================================
real scalar _aggte_requires_pass1_bw(
    real matrix gteval,
    real scalar gbar,
    string scalar type,
    real colvector eval_points)
{
    real scalar i

    for (i = 1; i <= rows(eval_points); i++) {
        if (rows(didhetero_build_gteeval(gteval, gbar, type, J(1, 1, eval_points[i]))) > 1) {
            return(1)
        }
    }

    return(0)
}


// =============================================================================
// _didhetero_unique_sorted()
//
// Return unique sorted values from a column vector.
// Simple O(n^2) implementation suitable for small vectors (eval points).
//
// Parameters:
//   v - column vector of real values
//
// Returns:
//   sorted column vector of unique values
// =============================================================================
real colvector _didhetero_unique_sorted(real colvector v)
{
    real colvector sorted_v, result
    real scalar n, i, n_unique

    n = rows(v)
    if (n == 0) return(J(0, 1, .))
    if (n == 1) return(v)

    // Sort the input
    sorted_v = sort(v, 1)

    // Count unique values
    n_unique = 1
    for (i = 2; i <= n; i++) {
        if (sorted_v[i] != sorted_v[i - 1]) {
            n_unique++
        }
    }

    // Extract unique values
    result = J(n_unique, 1, .)
    result[1] = sorted_v[1]
    n_unique = 1
    for (i = 2; i <= n; i++) {
        if (sorted_v[i] != sorted_v[i - 1]) {
            n_unique++
            result[n_unique] = sorted_v[i]
        }
    }

    return(result)
}


// =============================================================================
// didhetero_aggte_main()
//
// Main orchestration function for the aggte_gt command. Implements the
// three-phase aggregation pipeline: (1) Pass 1 estimation for all evaluation
// points including triplet filtering, weight computation, influence functions,
// and bandwidth selection; (2) Pass 2 re-estimation with optimal bandwidths
// for non-shortcircuit evaluation points, including standard error computation
// and bootstrap Phase A; (3) Bootstrap Phase B for critical value computation.
// =============================================================================
struct AggtResult scalar didhetero_aggte_main(
    string scalar type,
    real colvector eval_points,
    real scalar bstrap,
    real scalar biters,
    real scalar seed,
    real scalar porder,
    string scalar kernel,
    string scalar bwselect,
    real colvector bw_manual,
    real scalar uniformall,
    real scalar alp,
    pointer(real matrix) colvector B_g_t,
    real matrix G_g,
    real colvector Z,
    real matrix mu_G_g,
    real matrix gteval,
    real matrix catt_est,
    real matrix catt_se,
    real colvector zeval,
    real colvector bw_catt,
    real colvector kd0_Z,
    real colvector kd1_Z,
    real colvector Z_supp,
    real scalar gbar,
    real scalar n,
    real colvector c_hat_catt,
    real colvector c_check_catt,
    struct DidHeteroData scalar pass2_data,
    real matrix pass2_gps_mat,
    real matrix pass2_or_mat)
{
    // --- Local declarations ---
    struct AggtResult scalar result
    real scalar num_eval, num_zeval, num_gteval
    real scalar id_eval, id_gt, k
    real scalar e1, h_agg, h_common, num_gte

    real matrix gteeval
    real colvector gt_indices
    real colvector h_agg_vec, shortcircuit

    // Pass 1 per-eval temporaries
    real matrix aggte_weight, J_mat, mu_G_g_sub, G_g_sub, catt_est_sub
    real colvector aggte_kappa, aggte_est_e
    pointer(real matrix) colvector B_g_t_sub, xi_g_t

    // Pass 2 per-eval temporaries
    real matrix J_pass2, aggte_weight_p2, catt_est_new
    real matrix G_g_new, mu_G_g_new
    real colvector aggte_est_p2, aggte_kappa_p2
    pointer(real matrix) colvector B_g_t_cv

    // SE temporaries
    real colvector se_e, mathcal_V_e, ci1_lower_e, ci1_upper_e
    real matrix U_hat_e

    // Bootstrap temporaries
    real matrix mb_weight, mb_sup_t
    real colvector ci2_lower_e, ci2_upper_e

    // Extract dimensions and initialize result
    num_eval   = rows(eval_points)
    num_zeval  = rows(zeval)
    num_gteval = rows(gteval)

    result.aggte_est = J(num_eval, num_zeval, .)
    result.aggte_se  = J(num_eval, num_zeval, .)
    result.ci1_lower = J(num_eval, num_zeval, .)
    result.ci1_upper = J(num_eval, num_zeval, .)
    result.ci2_lower = J(num_eval, num_zeval, .)
    result.ci2_upper = J(num_eval, num_zeval, .)
    result.aggte_bw  = J(num_eval, 1, .)
    result.eval_info = eval_points
    result.zeval     = zeval
    result.type      = type

    // Phase control vectors
    h_agg_vec    = J(num_eval, 1, .)
    shortcircuit = J(num_eval, 1, 0)

    printf("{txt}aggte_gt: starting aggregation (%s, %g eval points)\n",
           type, num_eval)

    // Pre-generate bootstrap weights (shared across all eval points)
    if (bstrap) {
        if (seed >= 0 & seed < .) {
            rseed(seed)
        }
        mb_weight = _didhetero_mammen_weights_batch(biters, n)
        mb_sup_t  = J(biters, num_eval, .)
        printf("{txt}  Bootstrap: %g iterations, weights pre-generated\n",
               biters)
    }

    // Phase 1: Pass 1 estimation — triplet filtering, weights, xi/J, BW selection
    printf("{txt}  Phase 1: Pass 1 estimation...\n")

    for (id_eval = 1; id_eval <= num_eval; id_eval++) {

        e1 = eval_points[id_eval]

        // Triplet filtering
        gteeval = didhetero_build_gteeval(gteval, gbar, type,
                      J(1, 1, e1))
        num_gte = rows(gteeval)

        // Skip if no matching triples
        if (num_gte == 0) {
            printf("{txt}    eval[%g]=%g: no matching (g,t), skipping\n",
                   id_eval, e1)
            shortcircuit[id_eval] = 1
            continue
        }

        // Short-circuit: num_gte == 1, directly use CATT-level results
        if (num_gte == 1) {
            id_gt = _aggte_find_gt_index(gteeval[1,1], gteeval[1,2],
                        gteval)
            if (id_gt > 0) {
                result.aggte_est[id_eval, .] = catt_est[., id_gt]'
                result.aggte_se[id_eval, .]  = catt_se[., id_gt]'
                if (id_gt <= rows(c_hat_catt) & c_hat_catt[id_gt] < .) {
                    result.ci1_lower[id_eval, .] = (catt_est[., id_gt] -
                        c_hat_catt[id_gt] :* catt_se[., id_gt])'
                    result.ci1_upper[id_eval, .] = (catt_est[., id_gt] +
                        c_hat_catt[id_gt] :* catt_se[., id_gt])'
                }
                else {
                    result.ci1_lower[id_eval, .] = J(1, num_zeval, .)
                    result.ci1_upper[id_eval, .] = J(1, num_zeval, .)
                    printf("{txt}Warning: upstream analytical c_hat is unavailable for short-circuit aggte eval[%g]=%g; leaving analytical CI missing\n",
                        id_eval, e1)
                }
                // Bandwidth for short-circuit evals:
                // - If bwselect=="manual": honor user-specified scalar/vector
                // - Else: inherit CATT-level bandwidth
                if (bwselect == "manual") {
                    if (rows(bw_manual) == 1) {
                        result.aggte_bw[id_eval] = bw_manual[1]
                        h_agg_vec[id_eval] = bw_manual[1]
                    }
                    else if (rows(bw_manual) == num_eval) {
                        result.aggte_bw[id_eval] = bw_manual[id_eval]
                        h_agg_vec[id_eval] = bw_manual[id_eval]
                    }
                    else {
                        _error(3001, "aggte manual bw: length(bw) must be 1 or num_eval (short-circuit)")
                    }
                }
                else {
                    result.aggte_bw[id_eval] = bw_catt[id_gt]
                    h_agg_vec[id_eval] = bw_catt[id_gt]
                }
            }
            // Copy CATT-level bootstrap CI2 for short-circuit evals
            if (bstrap & id_gt > 0 & id_gt <= rows(c_check_catt)) {
                if (c_check_catt[id_gt] < .) {
                    result.ci2_lower[id_eval, .] = (catt_est[., id_gt] -
                        c_check_catt[id_gt] :* catt_se[., id_gt])'
                    result.ci2_upper[id_eval, .] = (catt_est[., id_gt] +
                        c_check_catt[id_gt] :* catt_se[., id_gt])'
                }
            }
            shortcircuit[id_eval] = 1
            printf("{txt}    eval[%g]=%g: num_gte=1, short-circuit\n",
                   id_eval, e1)
            continue
        }

        // Extract submatrices for this eval's (g,t) pairs
        gt_indices   = _aggte_find_gt_indices(gteeval, gteval)
        mu_G_g_sub   = mu_G_g[., gt_indices]
        G_g_sub      = G_g[., gt_indices]
        catt_est_sub = catt_est[., gt_indices]
        B_g_t_sub    = J(num_gte, 1, NULL)
        for (k = 1; k <= num_gte; k++) {
            B_g_t_sub[k] = B_g_t[gt_indices[k]]
        }

        // Weight computation
        didhetero_aggte_weights(mu_G_g_sub, gteeval, type, gbar,
                                num_zeval, num_gte,
                                aggte_weight, aggte_kappa)

        // Influence function computation
        didhetero_aggte_xi(G_g_sub, mu_G_g_sub, aggte_weight,
                           aggte_kappa, type, n, num_zeval, num_gte,
                           xi_g_t)

        didhetero_aggte_J(B_g_t_sub, xi_g_t, aggte_weight,
                          catt_est_sub, n, num_zeval, num_gte,
                          J_mat, aggte_est_e)

        // Pass 1 bandwidth selection
        // Manual bandwidth handling:
        // - If bwselect == "manual" and length(bw) == 1 → use scalar for all eval
        // - If bwselect == "manual" and length(bw) == num_eval → per-eval values
        // - Otherwise → compute IMSE/US1 bandwidth
        if (bwselect == "manual") {
            if (rows(bw_manual) == 1) {
                h_agg_vec[id_eval] = bw_manual[1]
            }
            else if (rows(bw_manual) == num_eval) {
                h_agg_vec[id_eval] = bw_manual[id_eval]
            }
            else {
                _error(3001, "aggte manual bw: length(bw) must be 1 or num_eval")
            }
        }
        else {
            h_agg_vec[id_eval] = didhetero_aggte_bw_pass1(
                J_mat, Z, zeval, kd0_Z, kd1_Z,
                n, bwselect, kernel, uniformall)
        }

        printf("{txt}    eval[%g]=%g: num_gte=%g, h_pass1=%9.6f\n",
               id_eval, e1, num_gte, h_agg_vec[id_eval])
    }

    // Phase 1.5: uniformall common bandwidth
    if (uniformall == 1 | (bwselect == "manual" & rows(bw_manual) == 1)) {
        // If manual scalar provided, use the scalar for all evals,
        // regardless of any short-circuit CATT bandwidths.
        if (bwselect == "manual" & rows(bw_manual) == 1) {
            h_common = bw_manual[1]
        }
        else {
            h_common = .
            for (id_eval = 1; id_eval <= num_eval; id_eval++) {
                if (h_agg_vec[id_eval] < .) {
                    if (h_agg_vec[id_eval] < h_common) {
                        h_common = h_agg_vec[id_eval]
                    }
                }
            }
        }
        if (h_common < .) {
            for (id_eval = 1; id_eval <= num_eval; id_eval++) {
                h_agg_vec[id_eval] = h_common
            }
            printf("{txt}  uniformall: common h = %9.6f\n", h_common)
        }
    }

    result.aggte_bw = h_agg_vec

    // Phase 2: Pass 2 re-estimation with optimal bandwidth
    printf("{txt}  Phase 2: Pass 2 re-estimation + SE...\n")

    for (id_eval = 1; id_eval <= num_eval; id_eval++) {

        if (shortcircuit[id_eval]) continue

        e1    = eval_points[id_eval]
        h_agg = h_agg_vec[id_eval]

        // Skip if bandwidth is missing
        if (h_agg >= . | h_agg <= 0) {
            printf("{txt}    eval[%g]=%g: bw missing, skipping\n",
                   id_eval, e1)
            continue
        }

        // Rebuild gteeval for this eval point
        gteeval = didhetero_build_gteeval(gteval, gbar, type,
                      J(1, 1, e1))
        num_gte = rows(gteeval)

        // Pass 2 re-estimation
        didhetero_aggte_pass2(
            pass2_data, pass2_gps_mat, pass2_or_mat,
            gteeval, h_agg, type, gbar,
            porder, kernel,
            aggte_est_p2, J_pass2, aggte_weight_p2, aggte_kappa_p2,
            catt_est_new, B_g_t_cv, G_g_new, mu_G_g_new)

        // Standard errors
        didhetero_aggte_se(
            J_pass2, aggte_est_p2, Z, zeval, Z_supp, kd0_Z,
            h_agg, n, porder, kernel, alp, num_zeval,
            se_e, mathcal_V_e, U_hat_e, ci1_lower_e, ci1_upper_e)

        // Store results
        result.aggte_est[id_eval, .] = aggte_est_p2'
        result.aggte_se[id_eval, .]  = se_e'
        result.ci1_lower[id_eval, .] = ci1_lower_e'
        result.ci1_upper[id_eval, .] = ci1_upper_e'

        // Bootstrap Phase A
        if (bstrap) {
            didhetero_aggte_bootstrap(
                aggte_est_p2, se_e, U_hat_e,
                Z, zeval, kd0_Z, h_agg,
                n, porder, kernel, biters, alp, uniformall,
                num_zeval, num_eval, id_eval,
                mb_weight,
                mb_sup_t)
        }

        printf("{txt}    eval[%g]=%g: h=%9.6f, est[1]=%9.6f, se[1]=%9.6f\n",
               id_eval, e1, h_agg, aggte_est_p2[1], se_e[1])
    }

    // =====================================================================
    // PHASE 3: Bootstrap CI2 for all eval points
    // =====================================================================
    if (bstrap) {
        _aggte_fill_ci2(result, mb_sup_t, shortcircuit,
                        num_eval, num_zeval, biters, alp, uniformall)
    }

    printf("{txt}aggte_gt: aggregation complete\n")

    return(result)
}


// =============================================================================
// _aggte_find_gt_index()
// Find the row index of a (g,t) pair in the full gteval matrix.
// Returns 0 if not found.
// =============================================================================
real scalar _aggte_find_gt_index(
    real scalar g, real scalar t, real matrix gteval)
{
    real scalar i
    for (i = 1; i <= rows(gteval); i++) {
        if (gteval[i, 1] == g & gteval[i, 2] == t) return(i)
    }
    return(0)
}


// =============================================================================
// _aggte_find_gt_indices()
// Find row indices in full gteval for all (g,t) pairs in gteeval.
// gteeval is M x 3 (g, t, eval), gteval is K x 2 (g, t).
// Returns M x 1 vector of indices.
// =============================================================================
real colvector _aggte_find_gt_indices(
    real matrix gteeval, real matrix gteval)
{
    real scalar m, i
    real colvector indices
    m = rows(gteeval)
    indices = J(m, 1, 0)
    for (i = 1; i <= m; i++) {
        indices[i] = _aggte_find_gt_index(gteeval[i, 1],
                         gteeval[i, 2], gteval)
    }
    return(indices)
}


// =============================================================================
// _aggte_fill_ci2()
//
// Fill bootstrap confidence intervals (CI2) for all evaluation points using
// accumulated supremum t-statistics. For uniformall=1, a single global critical
// value is applied to all evaluation points including short-circuit cases.
// For uniformall=0, per-evaluation critical values are used.
// =============================================================================
void _aggte_fill_ci2(
    struct AggtResult scalar result,
    real matrix mb_sup_t,
    real colvector shortcircuit,
    real scalar num_eval,
    real scalar num_zeval,
    real scalar biters,
    real scalar alp,
    real scalar uniformall)
{
    real scalar id_eval, b
    real scalar c_global, c_per_eval
    real colvector mb_sup_t_global, mb_sup_t_col
    real colvector est_row, se_row

    if (uniformall == 1) {
        // Global UCB: sup-t over all eval points with bootstrap data
        // (short-circuit evals have missing mb_sup_t, handled by
        //  didhetero_max_nonmissing which ignores missing values)
        mb_sup_t_global = J(biters, 1, .)
        for (b = 1; b <= biters; b++) {
            mb_sup_t_global[b] = didhetero_max_nonmissing(
                mb_sup_t[b, .]')
        }
        c_global = _didhetero_bs_quantile(mb_sup_t_global, 1 - alp)

        // Apply global critical value to ALL evals including short-circuit
        for (id_eval = 1; id_eval <= num_eval; id_eval++) {
            est_row = result.aggte_est[id_eval, .]'
            se_row  = result.aggte_se[id_eval, .]'
            // Skip evals with missing est/se (num_gte==0 case)
            if (est_row[1] >= . | se_row[1] >= .) continue
            result.ci2_lower[id_eval, .] =
                (est_row - c_global :* se_row)'
            result.ci2_upper[id_eval, .] =
                (est_row + c_global :* se_row)'
        }
    }
    else {
        // Per-eval UCB: short-circuit evals have no bootstrap data
        // (mb_sup_t column is all missing -> quantile returns missing -> ci2 stays missing)
        for (id_eval = 1; id_eval <= num_eval; id_eval++) {
            if (shortcircuit[id_eval]) continue
            mb_sup_t_col = mb_sup_t[., id_eval]
            c_per_eval = _didhetero_bs_quantile(mb_sup_t_col,
                             1 - alp)
            est_row = result.aggte_est[id_eval, .]'
            se_row  = result.aggte_se[id_eval, .]'
            result.ci2_lower[id_eval, .] =
                (est_row - c_per_eval :* se_row)'
            result.ci2_upper[id_eval, .] =
                (est_row + c_per_eval :* se_row)'
        }
    }
}


// =============================================================================
// _didhetero_aggte_ado_entry()
//
// Entry point called from aggte_gt.ado. Reads CATT-level results from Stata
// e() matrices, executes the aggregation pipeline, and stores results back
// into Stata matrices for the ADO layer to process.
//
// Arguments (passed from ADO):
//   type         - aggregation type
//   eval_matname - name of Stata matrix holding evaluation points
//   bstrap       - 1 to enable bootstrap, 0 otherwise
//   biters       - bootstrap iterations
//   seed         - random seed for bootstrap
//   porder       - polynomial order (1 or 2)
//   kernel       - kernel type
//   bwselect     - bandwidth selection method
//   bw_matname   - name of Stata matrix for manual bandwidth (empty for auto)
//   uniformall   - 1 for global UCB, 0 for per-evaluation
//   alp          - significance level
//   control_group- control group specification
//   anticipation - anticipation period
//
// Side effects:
//   Stores results in Stata matrices: __aggte_Estimate, __aggte_est,
//   __aggte_se, __aggte_ci1_lower, __aggte_ci1_upper, __aggte_ci2_lower,
//   __aggte_ci2_upper, __aggte_bw, __aggte_eval, __aggte_zeval
// =============================================================================
void _didhetero_aggte_ado_entry(
    string scalar type,
    string scalar eval_matname,
    real scalar bstrap,
    real scalar biters,
    real scalar seed,
    real scalar porder,
    string scalar kernel,
    string scalar bwselect,
    string scalar bw_matname,
    real scalar uniformall,
    real scalar alp,
    string scalar control_group,
    real scalar anticipation)
{
    // --- Local declarations ---
    struct AggtResult scalar R
    struct DidHeteroData scalar pass2_data
    real matrix B_g_t_flat, G_g, mu_G_g, gteval_mat, catt_est_mat, catt_se_mat
    real matrix Z_mat, zeval_mat, bw_mat, kd0_mat, kd1_mat, Z_supp_mat
    real matrix Y_wide_mat, G_unit_mat, t_vals_mat, gps_mat, or_mat
    real colvector Z, zeval, bw_catt, kd0_Z, kd1_Z, Z_supp, eval_points
    real colvector t_vals_vec
    pointer(real matrix) colvector B_g_t
    real scalar n, num_zeval, num_gteval, gbar
    real scalar num_eval, id_eval, r, row_idx
    real matrix Estimate
    real colvector bw_manual

    // Read all e() matrices into Mata
    B_g_t_flat = st_matrix("e(B_g_t)")
    G_g        = st_matrix("e(G_g)")
    Z_mat      = st_matrix("e(Z)")
    mu_G_g     = st_matrix("e(mu_G_g)")
    gteval_mat = st_matrix("e(gteval)")
    catt_est_mat = st_matrix("e(catt_est)")
    catt_se_mat  = st_matrix("e(catt_se)")
    zeval_mat  = st_matrix("e(zeval)")
    bw_mat     = st_matrix("e(bw)")
    kd0_mat    = st_matrix("e(kd0_Z)")
    kd1_mat    = st_matrix("e(kd1_Z)")
    Z_supp_mat = st_matrix("e(Z_supp)")
    Y_wide_mat = st_matrix("e(dh_Y_wide)")
    G_unit_mat = st_matrix("e(dh_G_unit)")
    t_vals_mat = st_matrix("e(dh_t_vals)")
    gps_mat    = st_matrix("e(dh_gps_mat)")
    or_mat     = st_matrix("e(dh_or_mat)")
    gbar       = st_numscalar("e(gbar)")

    // Placeholder declarations for analytical UCB critical values (CATT level)
    real matrix c_hat_mat
    real colvector c_hat_catt

    // Normalize vector orientations to column vectors
    Z      = (cols(Z_mat) == 1)      ? Z_mat      : Z_mat'
    zeval  = (cols(zeval_mat) == 1)   ? zeval_mat  : zeval_mat'
    bw_catt = (cols(bw_mat) == 1)    ? bw_mat     : bw_mat'
    kd0_Z  = (cols(kd0_mat) == 1)    ? kd0_mat    : kd0_mat'
    kd1_Z  = (cols(kd1_mat) == 1)    ? kd1_mat    : kd1_mat'
    Z_supp = (cols(Z_supp_mat) == 1) ? Z_supp_mat : Z_supp_mat'

    // Read eval points from Stata temp matrix
    eval_points = st_matrix(eval_matname)
    if (cols(eval_points) > 1 & rows(eval_points) == 1) {
        eval_points = eval_points'
    }
    // Read manual bandwidth vector (optional)
    if (bw_matname != "") {
        real matrix _bw_tmp
        _bw_tmp = st_matrix(bw_matname)
        bw_manual = (cols(_bw_tmp) == 1) ? _bw_tmp : _bw_tmp'
    }
    else {
        bw_manual = J(0, 1, .) // empty indicates auto-selection when not manual
    }

    // Extract dimensions
    n           = rows(Z)
    num_zeval   = rows(zeval)
    num_gteval  = rows(gteval_mat)
    num_eval    = rows(eval_points)

    // Robustness: read e(c_hat) if available, else create missing vector
    // Expected shapes: 1 x K (row) or K x 1 (col). Fallback to J(K,1,.) when absent/mismatch.
    c_hat_mat = st_matrix("e(c_hat)")
    if (rows(c_hat_mat) == 1 & cols(c_hat_mat) == num_gteval) {
        c_hat_catt = c_hat_mat'
    }
    else if (rows(c_hat_mat) == num_gteval & cols(c_hat_mat) == 1) {
        c_hat_catt = c_hat_mat
    }
    else {
        c_hat_catt = J(num_gteval, 1, .)
    }

    // Read CATT-level bootstrap critical values e(c_check) for shortcircuit ci2
    real matrix c_check_mat
    real colvector c_check_catt
    c_check_mat = st_matrix("e(c_check)")
    if (rows(c_check_mat) == 1 & cols(c_check_mat) == num_gteval) {
        c_check_catt = c_check_mat'
    }
    else if (rows(c_check_mat) == num_gteval & cols(c_check_mat) == 1) {
        c_check_catt = c_check_mat
    }
    else {
        c_check_catt = J(num_gteval, 1, .)
    }

    // Dimension validation
    real scalar __err
    __err = 0
    // gteval: K x 2
    if (cols(gteval_mat) != 2) {
        printf("{err}aggte_gt Mata: e(gteval) must have 2 columns, found %g\n", cols(gteval_mat))
        __err++
    }
    // B_g_t_flat: n x (R*K)
    if (rows(B_g_t_flat) != n) {
        printf("{err}aggte_gt Mata: e(B_g_t) rows=%g, expected n=%g\n", rows(B_g_t_flat), n)
        __err++
    }
    if (cols(B_g_t_flat) != num_zeval * num_gteval) {
        printf("{err}aggte_gt Mata: e(B_g_t) cols=%g, expected R*K=%g\n", cols(B_g_t_flat), num_zeval * num_gteval)
        __err++
    }
    // G_g: n x K
    if (rows(G_g) != n) {
        printf("{err}aggte_gt Mata: e(G_g) rows=%g, expected n=%g\n", rows(G_g), n)
        __err++
    }
    if (cols(G_g) != num_gteval) {
        printf("{err}aggte_gt Mata: e(G_g) cols=%g, expected K=%g\n", cols(G_g), num_gteval)
        __err++
    }
    // mu_G_g: R x K
    if (rows(mu_G_g) != num_zeval) {
        printf("{err}aggte_gt Mata: e(mu_G_g) rows=%g, expected R=%g\n", rows(mu_G_g), num_zeval)
        __err++
    }
    if (cols(mu_G_g) != num_gteval) {
        printf("{err}aggte_gt Mata: e(mu_G_g) cols=%g, expected K=%g\n", cols(mu_G_g), num_gteval)
        __err++
    }
    // catt_est / catt_se: R x K
    if (rows(catt_est_mat) != num_zeval | cols(catt_est_mat) != num_gteval) {
        printf("{err}aggte_gt Mata: e(catt_est) %gx%g, expected %gx%g\n", rows(catt_est_mat), cols(catt_est_mat), num_zeval, num_gteval)
        __err++
    }
    if (rows(catt_se_mat) != num_zeval | cols(catt_se_mat) != num_gteval) {
        printf("{err}aggte_gt Mata: e(catt_se) %gx%g, expected %gx%g\n", rows(catt_se_mat), cols(catt_se_mat), num_zeval, num_gteval)
        __err++
    }
    // bw_catt: K x 1
    if (rows(bw_catt) != num_gteval) {
        printf("{err}aggte_gt Mata: e(bw) length=%g, expected K=%g\n", rows(bw_catt), num_gteval)
        __err++
    }
    // kd0_Z / kd1_Z: R x 1
    if (rows(kd0_Z) != num_zeval) {
        printf("{err}aggte_gt Mata: e(kd0_Z) length=%g, expected R=%g\n", rows(kd0_Z), num_zeval)
        __err++
    }
    if (rows(kd1_Z) != num_zeval) {
        printf("{err}aggte_gt Mata: e(kd1_Z) length=%g, expected R=%g\n", rows(kd1_Z), num_zeval)
        __err++
    }
    if (__err > 0) {
        _error(198, sprintf("aggte_gt: %g dimension mismatch(es) detected", __err))
    }

    // Persisted Pass 2 inputs must be available for re-estimation
    if (rows(Y_wide_mat) != n) {
        _error(198, sprintf("aggte_gt: e(dh_Y_wide) rows=%g, expected n=%g", rows(Y_wide_mat), n))
    }
    if (rows(G_unit_mat) != n & cols(G_unit_mat) != n) {
        _error(198, sprintf("aggte_gt: e(dh_G_unit) length mismatch, expected n=%g", n))
    }
    t_vals_vec = (cols(t_vals_mat) == 1) ? t_vals_mat : t_vals_mat'
    if (cols(Y_wide_mat) != rows(t_vals_vec)) {
        _error(198, "aggte_gt: e(dh_t_vals) length must equal number of columns in e(dh_Y_wide)")
    }
    if (rows(gps_mat) == 0 | rows(or_mat) == 0) {
        _error(198, "aggte_gt: persisted Pass 2 inputs are incomplete; please re-run catt_gt or didhetero")
    }

    // Convert B_g_t from flat to pointer array
    B_g_t = _didhetero_unflatten_B_g_t(B_g_t_flat, n, num_zeval, num_gteval)

    // Validate unflattened B_g_t pointer array
    real scalar __k
    for (__k = 1; __k <= num_gteval; __k++) {
        if (B_g_t[__k] == NULL) {
            _error(198, sprintf("aggte_gt: B_g_t[%g] is NULL after unflatten", __k))
        }
        if (rows(*B_g_t[__k]) != n | cols(*B_g_t[__k]) != num_zeval) {
            printf("{err}aggte_gt: B_g_t[%g] is %gx%g, expected %gx%g\n", __k, rows(*B_g_t[__k]), cols(*B_g_t[__k]), n, num_zeval)
            _error(198, "aggte_gt: B_g_t unflatten dimension error")
        }
    }

    // Rehydrate the minimum DidHeteroData fields required by aggte Pass 2.
    pass2_data.Y_wide = Y_wide_mat
    pass2_data.G = (cols(G_unit_mat) == 1) ? G_unit_mat : G_unit_mat'
    pass2_data.Z = Z
    pass2_data.zeval = zeval
    pass2_data.t_vals = t_vals_vec
    pass2_data.n = n
    pass2_data.num_zeval = num_zeval
    pass2_data.T_num = cols(Y_wide_mat)
    pass2_data.period1 = (rows(t_vals_vec) > 0) ? t_vals_vec[1] : .
    pass2_data.control_group = control_group
    pass2_data.anticipation = anticipation

    // Call the orchestrator
    R = didhetero_aggte_main(
        type, eval_points, bstrap, biters, seed, porder,
        kernel, bwselect, bw_manual, uniformall, alp,
        B_g_t, G_g, Z, mu_G_g, gteval_mat,
        catt_est_mat, catt_se_mat, zeval, bw_catt,
        kd0_Z, kd1_Z, Z_supp, gbar, n,
        c_hat_catt, c_check_catt, pass2_data, gps_mat, or_mat)

    // Build combined results matrix.
    //
    // Schema follows Imai, Qin, and Yanagi (2025, Section 5):
    //   - type "simple":   Theta^O(z) has no aggregation-dimension index,
    //                      so Estimate is num_zeval x 8 with columns
    //                      (z, est, se, ci1_lower, ci1_upper,
    //                       ci2_lower, ci2_upper, bw).
    //   - other types: Estimate is (num_eval * num_zeval) x 9 with
    //                      leading column (eval) holding the
    //                      aggregation index (event time, cohort, or
    //                      calendar period).
    if (type == "simple") {
        Estimate = J(num_zeval, 8, .)
        for (r = 1; r <= num_zeval; r++) {
            Estimate[r, 1] = R.zeval[r]
            Estimate[r, 2] = R.aggte_est[1, r]
            Estimate[r, 3] = R.aggte_se[1, r]
            Estimate[r, 4] = R.ci1_lower[1, r]
            Estimate[r, 5] = R.ci1_upper[1, r]
            Estimate[r, 6] = R.ci2_lower[1, r]
            Estimate[r, 7] = R.ci2_upper[1, r]
            Estimate[r, 8] = R.aggte_bw[1]
        }
    }
    else {
        Estimate = J(num_eval * num_zeval, 9, .)
        for (id_eval = 1; id_eval <= num_eval; id_eval++) {
            for (r = 1; r <= num_zeval; r++) {
                row_idx = (id_eval - 1) * num_zeval + r
                Estimate[row_idx, 1] = R.eval_info[id_eval]
                Estimate[row_idx, 2] = R.zeval[r]
                Estimate[row_idx, 3] = R.aggte_est[id_eval, r]
                Estimate[row_idx, 4] = R.aggte_se[id_eval, r]
                Estimate[row_idx, 5] = R.ci1_lower[id_eval, r]
                Estimate[row_idx, 6] = R.ci1_upper[id_eval, r]
                Estimate[row_idx, 7] = R.ci2_lower[id_eval, r]
                Estimate[row_idx, 8] = R.ci2_upper[id_eval, r]
                Estimate[row_idx, 9] = R.aggte_bw[id_eval]
            }
        }
    }

    // Store results as Stata matrices
    st_matrix("__aggte_Estimate",  Estimate)
    st_matrix("__aggte_est",       R.aggte_est)
    st_matrix("__aggte_se",        R.aggte_se)
    st_matrix("__aggte_ci1_lower", R.ci1_lower)
    st_matrix("__aggte_ci1_upper", R.ci1_upper)
    st_matrix("__aggte_ci2_lower", R.ci2_lower)
    st_matrix("__aggte_ci2_upper", R.ci2_upper)
    st_matrix("__aggte_bw",        R.aggte_bw)
    st_matrix("__aggte_eval",      R.eval_info)
    st_matrix("__aggte_zeval",     R.zeval)
}


// =============================================================================
// _didhetero_aggte_display_table()
//
// Display formatted aggregation results table. Reads results from Stata matrices
// __aggte_Estimate and prints a formatted table with evaluation points,
// point estimates, standard errors, and confidence intervals.
// =============================================================================
void _didhetero_aggte_display_table(
    string scalar type,
    real scalar porder,
    string scalar kernel,
    string scalar bwselect,
    real scalar alp,
    real scalar bstrap,
    real scalar biters,
    real scalar uniformall)
{
    real matrix Est
    real scalar nrows, r, is_simple
    real scalar eval_val, z_val, est_val, se_val
    real scalar ci1l, ci1u, ci2l, ci2u, bw_val
    string scalar bstrap_str, uniform_str

    Est = st_matrix("__aggte_Estimate")
    nrows = rows(Est)
    is_simple = (type == "simple")

    bstrap_str  = (bstrap)     ? sprintf("%g iterations", biters) : "disabled"
    uniform_str = (uniformall) ? "global" : "per-eval"

    // Header
    printf("\n")
    printf("{txt}{hline 78}\n")
    printf("{txt}Aggregated Treatment Effect Estimates (aggte_gt)\n")
    printf("{txt}{hline 78}\n")
    printf("{txt}  Type: {res}%s{txt} | Polynomial order: {res}%g{txt} | Kernel: {res}%s\n",
           type, porder, kernel)
    printf("{txt}  Bandwidth selection: {res}%s{txt} | Significance level: {res}%g\n",
           bwselect, alp)
    printf("{txt}  Bootstrap: {res}%s{txt} | Uniform band: {res}%s\n",
           bstrap_str, uniform_str)
    printf("{txt}{hline 78}\n")

    // Column headers (simple has no eval column)
    if (is_simple) {
        printf("{txt}  %8s | %9s | %9s | %9s | %9s |",
               "zeval", "est", "se", "ci1_lb", "ci1_ub")
    }
    else {
        printf("{txt}  %8s | %8s | %9s | %9s | %9s | %9s |",
               "eval", "zeval", "est", "se", "ci1_lb", "ci1_ub")
    }
    if (bstrap) {
        printf(" %9s | %9s |", "ci2_lb", "ci2_ub")
    }
    printf(" %9s\n", "bw")

    // Separator
    if (is_simple) {
        printf("{txt}  %s+%s+%s+%s+%s+",
               "{hline 9}", "{hline 10}", "{hline 10}",
               "{hline 10}", "{hline 10}")
    }
    else {
        printf("{txt}  %s+%s+%s+%s+%s+%s+",
               "{hline 9}", "{hline 9}", "{hline 10}", "{hline 10}",
               "{hline 10}", "{hline 10}")
    }
    if (bstrap) {
        printf("%s+%s+", "{hline 10}", "{hline 10}")
    }
    printf("%s\n", "{hline 10}")

    // Data rows
    for (r = 1; r <= nrows; r++) {
        if (is_simple) {
            z_val    = Est[r, 1]
            est_val  = Est[r, 2]
            se_val   = Est[r, 3]
            ci1l     = Est[r, 4]
            ci1u     = Est[r, 5]
            ci2l     = Est[r, 6]
            ci2u     = Est[r, 7]
            bw_val   = Est[r, 8]

            printf("{res}  %8.3f | %9.4f | %9.4f | %9.4f | %9.4f |",
                   z_val, est_val, se_val, ci1l, ci1u)
        }
        else {
            eval_val = Est[r, 1]
            z_val    = Est[r, 2]
            est_val  = Est[r, 3]
            se_val   = Est[r, 4]
            ci1l     = Est[r, 5]
            ci1u     = Est[r, 6]
            ci2l     = Est[r, 7]
            ci2u     = Est[r, 8]
            bw_val   = Est[r, 9]

            printf("{res}  %8.3f | %8.3f | %9.4f | %9.4f | %9.4f | %9.4f |",
                   eval_val, z_val, est_val, se_val, ci1l, ci1u)
        }
        if (bstrap) {
            printf(" %9.4f | %9.4f |", ci2l, ci2u)
        }
        printf(" %9.4f\n", bw_val)
    }

    // Footer
    printf("{txt}{hline 78}\n")
    if (is_simple) {
        printf("{txt}  %g z-eval points = %g rows (type=simple has no eval dimension)\n",
               rows(st_matrix("__aggte_zeval")), nrows)
    }
    else {
        printf("{txt}  %g eval points x %g z-eval points = %g rows\n",
               rows(st_matrix("__aggte_eval")),
               rows(st_matrix("__aggte_zeval")),
               nrows)
    }
    printf("{txt}{hline 78}\n")
}


end
