mata:

// =============================================================================
// didhetero_bwselect_lp.mata
// MSE-optimal bandwidth selection for local polynomial regression
//
// Provides:
//   - _didhetero_bw_even_optimize()     // Brent optimizer for even-order BW
//   - _didhetero_bw_imse_optimize()     // IMSE bandwidth optimizer
//   - _dh_lpbw_mse_details()            // Bias/variance pilot estimates
//   - _didhetero_lpbwselect_mse()       // MSE-DPI pointwise bandwidth
//   - _didhetero_lpbwselect_imse()      // IMSE-DPI global bandwidth
//   - _didhetero_nn_variance()          // Nearest-neighbor variance estimate
//   - _didhetero_nn_variance_fast()     // Fast NN variance (sorted input)
//   - _didhetero_bwcheck()              // Minimum bandwidth enforcement
//   - _didhetero_bwregul()              // Bandwidth regularization
//
// Requires:
//   - didhetero_kernel.mata         (didhetero_kernel_eval)
//   - didhetero_lpr.mata            (didhetero_lpr)
//   - didhetero_utils_formula.mata  (didhetero_selectindex)
//
// Paper reference: Section 3.2, IMSE-optimal bandwidth selection
// =============================================================================


real scalar _didhetero_bw_even_optimize(real scalar V_in,
                                        real scalar B1_in,
                                        real scalar B2_in,
                                        real scalar R_in,
                                        real scalar N,
                                        real scalar range_X,
                                        real scalar poly_order,
                                        real scalar deriv_order,
                                        real scalar bwregul_scale)
{
    real scalar exp1, exp2
    real scalar B2_eff, R_eff, a_eff, c_eff
    real scalar lower_bd, upper_bd
    real scalar a, b, x, w, v, fx, fw, fv
    real scalar d, e, xm, tol1, tol2, u, fu
    real scalar p_num, q_num, r_num, etemp
    real scalar i, cgold, zeps

    if (V_in == . | B1_in == . | N <= 0 | range_X <= 0) return(.)

    B2_eff = (B2_in < . ? B2_in : 0)
    R_eff = (R_in < . ? R_in : 0)
    a_eff = B1_in + bwregul_scale * R_eff
    c_eff = V_in / N

    exp1 = 2 * poly_order + 2 - 2 * deriv_order
    exp2 = 1 + 2 * deriv_order

    lower_bd = epsilon(1)
    upper_bd = range_X
    a = lower_bd
    b = upper_bd
    cgold = 0.3819660112501051
    zeps = 1e-15

    x = a + cgold * (b - a)
    w = x
    v = x
    fx = abs(x^exp1 * (a_eff + x * B2_eff)^2 + c_eff / x^exp2)
    fw = fx
    fv = fx
    d = 0
    e = 0

    for (i = 1; i <= 200; i++) {
        xm = 0.5 * (a + b)
        tol1 = 1e-12 * abs(x) + zeps
        tol2 = 2 * tol1

        if (abs(x - xm) <= (tol2 - 0.5 * (b - a))) break

        if (abs(e) > tol1) {
            r_num = (x - w) * (fx - fv)
            q_num = (x - v) * (fx - fw)
            p_num = (x - v) * q_num - (x - w) * r_num
            q_num = 2 * (q_num - r_num)

            if (q_num > 0) p_num = -p_num
            q_num = abs(q_num)
            etemp = e
            e = d

            if (abs(p_num) >= abs(0.5 * q_num * etemp) | p_num <= q_num * (a - x) | p_num >= q_num * (b - x)) {
                if (x >= xm) e = a - x
                else         e = b - x
                d = cgold * e
            }
            else {
                d = p_num / q_num
                u = x + d
                if (u - a < tol2 | b - u < tol2) {
                    if (xm - x >= 0) d = tol1
                    else             d = -tol1
                }
            }
        }
        else {
            if (x >= xm) e = a - x
            else         e = b - x
            d = cgold * e
        }

        if (abs(d) >= tol1) u = x + d
        else if (d >= 0)    u = x + tol1
        else                u = x - tol1

        fu = abs(u^exp1 * (a_eff + u * B2_eff)^2 + c_eff / u^exp2)

        if (fu <= fx) {
            if (u >= x) a = x
            else        b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
        }
        else {
            if (u < x) a = u
            else       b = u
            if (fu <= fw | w == x) {
                v = w
                fv = fw
                w = u
                fw = fu
            }
            else if (fu <= fv | v == x | v == w) {
                v = u
                fv = fu
            }
        }
    }

    return(x)
}


real scalar _didhetero_bw_imse_optimize(real scalar mean_Vh,
                                        real scalar mean_Bh,
                                        real scalar N,
                                        real scalar range_X,
                                        real scalar poly_order,
                                        real scalar deriv_order)
{
    real scalar exp1, exp2
    real scalar a_gs, b_gs, gr, c_gs, d_gs
    real scalar fc, fd, tol_gs
    real scalar i

    if (mean_Vh == . | mean_Bh == . | N <= 0 | range_X <= 0) return(.)
    if (mean_Vh <= 0 | mean_Bh <= 0) return(.)

    exp1 = 2 * poly_order + 2 - 2 * deriv_order
    exp2 = 1 + 2 * deriv_order
    gr = (sqrt(5) - 1) / 2
    a_gs = epsilon(1)
    b_gs = range_X
    tol_gs = 1e-12 * range_X

    c_gs = b_gs - gr * (b_gs - a_gs)
    d_gs = a_gs + gr * (b_gs - a_gs)

    for (i = 1; i <= 400; i++) {
        fc = abs(c_gs^exp1 * mean_Bh + mean_Vh / (N * c_gs^exp2))
        fd = abs(d_gs^exp1 * mean_Bh + mean_Vh / (N * d_gs^exp2))

        if (fc < fd) b_gs = d_gs
        else         a_gs = c_gs

        if (abs(b_gs - a_gs) <= tol_gs) break

        c_gs = b_gs - gr * (b_gs - a_gs)
        d_gs = a_gs + gr * (b_gs - a_gs)
    }

    return((a_gs + b_gs) / 2)
}


// --- 2.4b Pointwise MSE-DPI details ---
// Computes pointwise MSE-DPI bandwidth and variance/bias components
// for a single evaluation point.
void _dh_lpbw_mse_details(real colvector Y, real colvector X,
                          real scalar eval_pt, real scalar p,
                          real scalar deriv,
                          string scalar kernel,
                          real scalar h_out,
                          real scalar V_h,
                          real scalar B_h)
{
    real scalar N, q, even, C_c, x_iq, x_sd
    real scalar c_bw, bw_max_r, bw_min
    real scalar bw_mp2, bw_mp3, b_mse_dpi, h_val
    real scalar V_tmp, B1_tmp, B2_tmp, R_tmp
    real scalar j
    real colvector X_s, Y_s, dups, dupsid
    real colvector sort_idx
    real scalar range_X, x_min, x_max
    real scalar rV_h, rB_h

    N = rows(X)
    q = p + 1
    even = (mod(p - deriv, 2) == 0)

    // Pilot bandwidth constant
    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else C_c = 2.34

    x_sd = sqrt(variance(X))
    // Interquartile range
    {
        real colvector X_sorted
        real scalar q25, q75, idx_r, j_r, h_r
        X_sorted = sort(X, 1)
        idx_r = 1 + (N - 1) * 0.25
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q25 = X_sorted[N]
        else q25 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        idx_r = 1 + (N - 1) * 0.75
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q75 = X_sorted[N]
        else q75 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        x_iq = q75 - q25
    }

    x_min = min(X)
    x_max = max(X)
    range_X = x_max - x_min

    // Sort data
    sort_idx = order(X, 1)
    X_s = X[sort_idx]
    Y_s = Y[sort_idx]

    // Compute duplicate counts and IDs
    dups = J(N, 1, 0)
    dupsid = J(N, 1, 0)
    for (j = 1; j <= N; j++) {
        dups[j] = sum(X_s :== X_s[j])
    }
    j = 1
    while (j <= N) {
        {
            real scalar k_dup
            for (k_dup = 0; k_dup < dups[j]; k_dup++) {
                dupsid[j + k_dup] = k_dup + 1
            }
            j = j + dups[j]
        }
    }

    h_out = .
    V_h = .
    B_h = .

    // Maximum bandwidth at evaluation point
    bw_max_r = max((abs(eval_pt - x_min), abs(eval_pt - x_max)))

    // Pilot bandwidth constant
    c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
    c_bw = min((c_bw, bw_max_r))

    // Minimum bandwidth (bwcheck)
    {
        real colvector abs_dist
        abs_dist = sort(abs(X_s :- eval_pt), 1)
        if (21 <= N) bw_min = abs_dist[21]
        else bw_min = abs_dist[N]
    }
    c_bw = max((c_bw, bw_min))

    // Stage 1
    bw_mp2 = _didhetero_lprobust_bw(Y_s, X_s, eval_pt,
                q + 1, q + 1, q + 2,
                c_bw, range_X, range_X, 0,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (even) {
        bw_mp2 = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, 0,
                    N, range_X, q + 1, q + 1, 0)
    }
    if (bw_mp2 == .) bw_mp2 = range_X
    bw_mp2 = min((bw_mp2, bw_max_r))
    bw_mp2 = max((bw_mp2, bw_min))

    // Stage 2
    bw_mp3 = _didhetero_lprobust_bw(Y_s, X_s, eval_pt,
                q + 2, q + 2, q + 3,
                c_bw, range_X, range_X, 0,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (even) {
        bw_mp3 = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, 0,
                    N, range_X, q + 2, q + 2, 0)
    }
    if (bw_mp3 == .) bw_mp3 = range_X
    bw_mp3 = min((bw_mp3, bw_max_r))
    bw_mp3 = max((bw_mp3, bw_min))

    // Stage 3
    b_mse_dpi = _didhetero_lprobust_bw(Y_s, X_s, eval_pt,
                q, p + 1, q + 1,
                c_bw, bw_mp2, bw_mp3, 1,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (even) {
        b_mse_dpi = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, R_tmp,
                    N, range_X, q, p + 1, 1)
    }
    if (b_mse_dpi == .) b_mse_dpi = range_X
    b_mse_dpi = min((b_mse_dpi, bw_max_r))
    b_mse_dpi = max((b_mse_dpi, bw_min))

    // Stage 4
    h_val = _didhetero_lprobust_bw(Y_s, X_s, eval_pt,
                p, deriv, q,
                c_bw, b_mse_dpi, bw_mp2, 1,
                3, kernel, dups, dupsid,
                V_tmp, B1_tmp, B2_tmp, R_tmp)
    if (even) {
        h_val = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, R_tmp,
                    N, range_X, p, deriv, 1)
    }
    if (h_val == .) h_val = range_X
    h_val = min((h_val, bw_max_r))
    h_val = max((h_val, bw_min))
    h_out = h_val

    if (even == 0) {
        rV_h = 2 * deriv + 1
        rB_h = 2 * (p + 1 - deriv)
        V_h = rV_h * V_tmp
        B_h = rB_h * B1_tmp^2
    }
    else {
        V_h = V_tmp
        if (B2_tmp >= .) B2_tmp = 0
        B_h = (B1_tmp + h_val * B2_tmp)^2
    }
}


// --- 2.4c Pointwise MSE-DPI bandwidth selector ---
real colvector _didhetero_lpbwselect_mse(real colvector Y, real colvector X,
                                         real colvector eval, real scalar p,
                                         real scalar deriv,
                                         string scalar kernel)
{
    real scalar R_eval, r_idx, h_tmp, Vh_tmp, Bh_tmp
    real colvector h_mse

    R_eval = rows(eval)
    h_mse = J(R_eval, 1, .)

    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        _dh_lpbw_mse_details(Y, X, eval[r_idx], p, deriv,
            kernel, h_tmp, Vh_tmp, Bh_tmp)
        h_mse[r_idx] = h_tmp
    }

    return(h_mse)
}


// --- 2.5 IMSE-DPI common bandwidth selector ---
// Returns scalar bandwidth minimizing integrated MSE over evaluation grid.
//
// For odd (p - deriv): V_h = rV * V, B_h = rB * B1^2, closed-form solution.
// For even (p - deriv): numerical optimization of IMSE objective.
real scalar _didhetero_lpbwselect_imse(real colvector Y, real colvector X,
                                       real colvector eval, real scalar p,
                                       real scalar deriv,
                                       string scalar kernel)
{
    real scalar N, R_eval, r_idx
    real scalar r_exp, rB, rV
    real scalar mean_Vh, mean_Bh, h_imse
    real scalar range_X
    real colvector B_hat_sq, V_hat
    real colvector Vh, Bh
    real scalar n_valid_V, n_valid_B
    real scalar is_even
    real colvector eval_grid

    N = rows(X)
    range_X = max(X) - min(X)

    if (rows(eval) < 1) {
        errprintf("_didhetero_lpbwselect_imse: no evaluation points\n")
        return(.)
    }

    // IMSE-DPI uses 30-point equispaced grid
    eval_grid = rangen(min(X), max(X), 30)
    R_eval = rows(eval_grid)

    is_even = (mod(p - deriv, 2) == 0)

    Vh = J(R_eval, 1, .)
    Bh = J(R_eval, 1, .)
    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        real scalar h_tmp, Vh_tmp, Bh_tmp
        _dh_lpbw_mse_details(Y, X, eval_grid[r_idx], p, deriv,
            kernel, h_tmp, Vh_tmp, Bh_tmp)
        Vh[r_idx] = Vh_tmp
        Bh[r_idx] = Bh_tmp
    }

    if (is_even == 0) {
        rV = 2 * deriv + 1
        rB = 2 * (p + 1 - deriv)
        r_exp = 1 / (2 * p + 3)
    }

    // Compute means, skipping missing values
    mean_Vh = 0
    mean_Bh = 0
    n_valid_V = 0
    n_valid_B = 0
    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        if (Vh[r_idx] < .) {
            mean_Vh = mean_Vh + Vh[r_idx]
            n_valid_V = n_valid_V + 1
        }
        if (Bh[r_idx] < .) {
            mean_Bh = mean_Bh + Bh[r_idx]
            n_valid_B = n_valid_B + 1
        }
    }

    if (n_valid_V == 0 | n_valid_B == 0) return(.)
    mean_Vh = mean_Vh / n_valid_V
    mean_Bh = mean_Bh / n_valid_B

    if (mean_Bh < 1e-20) return(range_X / 2)
    if (mean_Vh <= 0) return(.)

    if (is_even == 0) {
        // ODD: closed-form
        h_imse = (mean_Vh / (N * mean_Bh))^r_exp
    }
    else {
        // EVEN: numerically optimize the IMSE objective
        // |H^(2p+2-2v)*mean(B_h) + mean(V_h)/(N*H^(1+2v))|.
        h_imse = _didhetero_bw_imse_optimize(mean_Vh, mean_Bh, N,
                    range_X, p, deriv)
    }

    if (h_imse >= . | h_imse <= 0) return(.)

    return(h_imse)
}


end
mata:

// =============================================================================
// didhetero_bwselect_lp.mata
// LP regression bandwidth selection routines (MSE-DPI and IMSE-DPI)
//
// Dependencies: didhetero_kernel.mata (must be compiled first)
//
// Contains:
//   - Shared utility functions (kernel integrals, NN variance, bwcheck, etc.)
//   - lpbwselect MSE-DPI bandwidth for local polynomial regression
//   - lpbwselect IMSE-DPI bandwidth for local polynomial regression
// =============================================================================

// =====================================================================
// Section 1: Shared utility functions
// =====================================================================

real scalar _didhetero_kernel_int(real scalar l, real scalar m,
                                  string scalar kernel)
{
    real scalar k, dbl_fact, i

    if (mod(l, 2) == 1) return(0)
    k = l / 2

    if (kernel == "epa") {
        if (m == 1) return(3 / ((2*k + 1) * (2*k + 3)))
        else if (m == 2) return((9/8) * (1/(2*k + 1) - 2/(2*k + 3) + 1/(2*k + 5)))
    }
    else if (kernel == "gau") {
        dbl_fact = 1
        for (i = 1; i <= k; i++) dbl_fact = dbl_fact * (2*i - 1)
        if (m == 1) return(dbl_fact)
        else if (m == 2) return(dbl_fact / (2^k * 2 * sqrt(c("pi"))))
    }
    _dh_error_input("bwselect_lp", "invalid kernel or m")
    return(.)
}

void _didhetero_bw_kernel_matrices(real scalar p, string scalar kernel,
                                   real matrix S_p, real matrix C_p,
                                   real colvector s_p1)
{
    real scalar i, j
    S_p  = J(p + 1, p + 1, .)
    C_p  = J(p + 1, p + 1, .)
    s_p1 = J(p + 1, 1, .)
    for (i = 0; i <= p; i++) {
        for (j = 0; j <= p; j++) {
            S_p[i + 1, j + 1] = _didhetero_kernel_int(i + j, 1, kernel)
            C_p[i + 1, j + 1] = _didhetero_kernel_int(i + j, 2, kernel)
        }
        s_p1[i + 1] = _didhetero_kernel_int(p + 1 + i, 1, kernel)
    }
}

real colvector _didhetero_nn_variance(real colvector Y, real colvector X,
                                      | real scalar nnmatch)
{
    real scalar n, i, j, J_nn, count, ss
    real colvector sigma2, dist_i, Y_nn
    real matrix sort_order
    real scalar Ybar_nn

    if (args() < 3) nnmatch = 3
    n = rows(X)
    J_nn = nnmatch
    if (n <= J_nn) J_nn = n - 1
    if (J_nn < 1) return(J(n, 1, 0))

    // Dispatch to O(n log n) algorithm for larger samples
    if (n > 100) return(_didhetero_nn_variance_fast(Y, X, J_nn))

    sigma2 = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        dist_i = abs(X :- X[i])
        dist_i[i] = .
        sort_order = order(dist_i, 1)
        Y_nn = J(J_nn, 1, .)
        count = 0
        for (j = 1; j <= n; j++) {
            if (dist_i[sort_order[j]] < .) {
                count++
                Y_nn[count] = Y[sort_order[j]]
                if (count >= J_nn) break
            }
        }
        if (count >= 1) {
            Ybar_nn = mean(Y_nn[1..count])
            ss = 0
            for (j = 1; j <= count; j++) ss = ss + (Y_nn[j] - Ybar_nn)^2
            sigma2[i] = ss / count
        }
        else sigma2[i] = 0
    }
    return(sigma2)
}

// O(n log n) nearest-neighbor variance via sort + window scan
// Equivalent to _didhetero_nn_variance() but avoids O(n^2 log n) brute force.
real colvector _didhetero_nn_variance_fast(real colvector Y, real colvector X,
                                            real scalar J_nn)
{
    real scalar n, i, lo, hi, count, j
    real scalar Ybar_nn, ss
    real colvector sort_idx, X_sorted, Y_sorted
    real colvector sigma2_sorted, sigma2
    real colvector Y_nn

    n = rows(X)

    // Step 1: Global sort by X
    sort_idx = order(X, 1)
    X_sorted = X[sort_idx]
    Y_sorted = Y[sort_idx]

    // Step 2: For each sorted position, find J_nn nearest neighbors via window
    sigma2_sorted = J(n, 1, .)

    for (i = 1; i <= n; i++) {
        // Expand left (lo) and right (hi) pointers to collect J_nn neighbors
        // excluding self (position i)
        lo = i - 1
        hi = i + 1
        count = 0
        Y_nn = J(J_nn, 1, .)

        while (count < J_nn) {
            if (lo >= 1 & hi <= n) {
                // Both sides available: pick the closer one
                if ((X_sorted[i] - X_sorted[lo]) <= (X_sorted[hi] - X_sorted[i])) {
                    count++
                    Y_nn[count] = Y_sorted[lo]
                    lo--
                } else {
                    count++
                    Y_nn[count] = Y_sorted[hi]
                    hi++
                }
            } else if (lo >= 1) {
                // Only left available
                count++
                Y_nn[count] = Y_sorted[lo]
                lo--
            } else if (hi <= n) {
                // Only right available
                count++
                Y_nn[count] = Y_sorted[hi]
                hi++
            } else {
                // No more neighbors (shouldn't happen if J_nn < n)
                break
            }
        }

        // Compute NN variance
        if (count >= 1) {
            Ybar_nn = 0
            for (j = 1; j <= count; j++) Ybar_nn = Ybar_nn + Y_nn[j]
            Ybar_nn = Ybar_nn / count
            ss = 0
            for (j = 1; j <= count; j++) ss = ss + (Y_nn[j] - Ybar_nn)^2
            sigma2_sorted[i] = ss / count
        } else {
            sigma2_sorted[i] = 0
        }
    }

    // Step 3: Restore original order
    sigma2 = J(n, 1, .)
    sigma2[sort_idx] = sigma2_sorted
    return(sigma2)
}

real colvector _didhetero_bwcheck(real colvector h, real colvector X,
                                  real colvector eval, | real scalar bwcheck)
{
    real scalar R_eval, r, n_eff, range_X
    if (args() < 4) bwcheck = 21
    R_eval = rows(eval)
    range_X = max(X) - min(X)
    for (r = 1; r <= R_eval; r++) {
        if (h[r] == . | h[r] <= 0) continue
        n_eff = sum(abs(X :- eval[r]) :<= h[r])
        while (n_eff < bwcheck & h[r] < range_X) {
            h[r] = h[r] * 1.1
            n_eff = sum(abs(X :- eval[r]) :<= h[r])
        }
    }
    return(h)
}

real colvector _didhetero_bwregul(real colvector h, real colvector X,
                                  real colvector eval, real scalar p,
                                  | real scalar bwregul)
{
    real scalar R_eval, r, n_eff, min_obs, range_X
    if (args() < 5) bwregul = 1
    R_eval = rows(eval)
    range_X = max(X) - min(X)
    min_obs = p + 2 + bwregul
    for (r = 1; r <= R_eval; r++) {
        if (h[r] == . | h[r] <= 0) continue
        n_eff = sum(abs(X :- eval[r]) :<= h[r])
        while (n_eff < min_obs & h[r] < range_X) {
            h[r] = h[r] * 1.1
            n_eff = sum(abs(X :- eval[r]) :<= h[r])
        }
    }
    return(h)
}

real colvector _didhetero_ecdf(real colvector X)
{
    real scalar n, i
    real colvector F_n
    n = rows(X)
    F_n = J(n, 1, .)
    for (i = 1; i <= n; i++) F_n[i] = sum(X :<= X[i]) / n
    return(F_n)
}


// =====================================================================
// Section 2: lpbwselect MSE-DPI bandwidth for local polynomial regression
//
// Internal helper functions:
//   _didhetero_lpbwselect_mse   - pointwise MSE-DPI bandwidth
//   _didhetero_lprobust_bw      - core bandwidth computation
//   _didhetero_lprobust_res     - nearest-neighbor residuals
//   _didhetero_lprobust_vce     - sandwich variance estimator
// =====================================================================

// --- 2.1 NN residuals on sorted data ---
// Computes nearest-neighbor residuals. X must be sorted with
// duplicate counts (dups) and duplicate IDs (dupsid) precomputed.
real colvector _didhetero_lprobust_res(real colvector X, real colvector Y,
                                       real scalar nnmatch,
                                       real colvector dups,
                                       real colvector dupsid)
{
    real scalar n, pos, rpos, lpos, Ji
    real scalar y_J
    real colvector res
    real scalar n_ind, lo, hi

    n = rows(X)
    res = J(n, 1, .)

    for (pos = 1; pos <= n; pos++) {
        rpos = dups[pos] - dupsid[pos]
        lpos = dupsid[pos] - 1

        while (lpos + rpos < min((nnmatch, n - 1))) {
            if (pos - lpos - 1 <= 0) {
                rpos = rpos + dups[pos + rpos + 1]
            }
            else if (pos + rpos + 1 > n) {
                lpos = lpos + dups[pos - lpos - 1]
            }
            else if ((X[pos] - X[pos - lpos - 1]) >
                     (X[pos + rpos + 1] - X[pos])) {
                rpos = rpos + dups[pos + rpos + 1]
            }
            else if ((X[pos] - X[pos - lpos - 1]) <
                     (X[pos + rpos + 1] - X[pos])) {
                lpos = lpos + dups[pos - lpos - 1]
            }
            else {
                rpos = rpos + dups[pos + rpos + 1]
                lpos = lpos + dups[pos - lpos - 1]
            }
        }

        lo = pos - lpos
        hi = min((n, pos + rpos))
            y_J = colsum(Y[lo..hi]) - Y[pos]
        Ji = (hi - lo + 1) - 1
        if (Ji > 0) {
            res[pos] = sqrt(Ji / (Ji + 1)) * (Y[pos] - y_J / Ji)
        }
        else {
            res[pos] = 0
        }
    }
    return(res)
}

// --- 2.2 Sandwich variance estimator (no clustering) ---
real matrix _didhetero_lprobust_vce(real matrix RX, real colvector res)
{
    return(cross(res :* RX, res :* RX))
}

// --- Trapezoidal integration ---
// Computes approximate integral via trapezoidal rule on sorted grid.
// Skips intervals containing missing values.
real scalar _didhetero_trapz(real colvector grid, real colvector values)
{
    real scalar R, r, result, valid_intervals, dz

    R = rows(grid)
    if (R < 2) return(.)

    result = 0
    valid_intervals = 0

    for (r = 1; r <= R - 1; r++) {
        // Skip intervals with any missing values
        if (grid[r] >= . | grid[r + 1] >= . | values[r] >= . | values[r + 1] >= .) {
            continue
        }
        dz = grid[r + 1] - grid[r]
        result = result + dz * (values[r] + values[r + 1]) / 2
        valid_intervals = valid_intervals + 1
    }

    if (valid_intervals == 0) return(.)
    return(result)
}

// --- 2.3 Core bandwidth computation ---
// Implements pre-asymptotic MSE-DPI bandwidth formula for nearest-neighbor
// variance estimation without clustering.
real scalar _didhetero_lprobust_bw(real colvector Y, real colvector X,
                                    real scalar c_eval,
                                    real scalar o, real scalar nu,
                                    real scalar o_B,
                                    real scalar h_V, real scalar h_B1,
                                    real scalar h_B2, real scalar scale,
                                    real scalar nnmatch,
                                    string scalar kernel,
                                    real colvector dups,
                                    real colvector dupsid,
                                    real scalar V_out,
                                    real scalar B1_out,
                                    real scalar B2_out,
                                    real scalar R_out)
{
    real scalar N, n_V, n_B, j, i
    real scalar V_V, BConst1, BConst2, BWreg, V_B
    real scalar V_final, r_exp, rB, rV, bw
    real colvector w, ind_V, eY, eX, eW
    real colvector dups_V, dupsid_V, res_V
    real colvector dups_B, dupsid_B, res_B
    real matrix R_V, invG_V, RW
    real colvector beta_V, Hp, v1, v2
    real colvector w_B, ind_B, eY_B, eX_B, eW_B
    real matrix R_B1, invG_B1
    real colvector beta_B1
    real matrix R_B2, invG_B2
    real colvector beta_B2
    real scalar cond_num, reg_lambda
    real matrix Gram_tmp
    real matrix vce_mat
    real colvector u_V

    N = rows(X)

    // ===== VARIANCE PART (uses h_V) =====
    u_V = (X :- c_eval) / h_V
    w = didhetero_kernel_eval(u_V, kernel) / h_V
    ind_V = (w :> 0)
    n_V = sum(ind_V)
    if (n_V < o + 2) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }

    eY = select(Y, ind_V)
    eX = select(X, ind_V)
    eW = select(w, ind_V)

    // Polynomial basis matrix
    R_V = J(n_V, o + 1, .)
    for (j = 1; j <= o + 1; j++) {
        R_V[., j] = (eX :- c_eval) :^ (j - 1)
    }

    // Inverse Gram matrix via Cholesky decomposition
    Gram_tmp = cross(R_V :* sqrt(eW), R_V :* sqrt(eW))
    // Condition number check and Tikhonov regularization
    cond_num = cond(Gram_tmp)
    if (cond_num != . & cond_num > 1e12) {
        reg_lambda = 1e-8 * trace(Gram_tmp) / cols(Gram_tmp)
        Gram_tmp = Gram_tmp + reg_lambda * I(cols(Gram_tmp))
    }
    invG_V = cholinv(Gram_tmp)
    if (hasmissing(invG_V)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }

    // NN residuals on the variance subsample
    dups_V = select(dups, ind_V)
    dupsid_V = select(dupsid, ind_V)
    res_V = _didhetero_lprobust_res(eX, eY, nnmatch, dups_V, dupsid_V)

    // Variance component
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

    // Bias component B1
    w_B = didhetero_kernel_eval((X :- c_eval) / h_B1, kernel)
    ind_B = (w_B :> 0)
    n_B = sum(ind_B)
    if (n_B < o_B + 2) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }

    eY_B = select(Y, ind_B)
    eX_B = select(X, ind_B)
    eW_B = select(w_B, ind_B)

    R_B1 = J(n_B, o_B + 1, .)
    for (j = 1; j <= o_B + 1; j++) {
        R_B1[., j] = (eX_B :- c_eval) :^ (j - 1)
    }
    Gram_tmp = cross(R_B1 :* sqrt(eW_B), R_B1 :* sqrt(eW_B))
    // Condition number check and Tikhonov regularization
    cond_num = cond(Gram_tmp)
    if (cond_num != . & cond_num > 1e12) {
        reg_lambda = 1e-8 * trace(Gram_tmp) / cols(Gram_tmp)
        Gram_tmp = Gram_tmp + reg_lambda * I(cols(Gram_tmp))
    }
    invG_B1 = cholinv(Gram_tmp)
    if (hasmissing(invG_B1)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }
    beta_B1 = invG_B1 * cross(R_B1, eW_B :* eY_B)

    // Regularization (if scale > 0)
    BWreg = 0
    if (scale > 0) {
        dups_B = select(dups, ind_B)
        dupsid_B = select(dupsid, ind_B)
        res_B = _didhetero_lprobust_res(eX_B, eY_B, nnmatch,
                                         dups_B, dupsid_B)
        vce_mat = _didhetero_lprobust_vce(R_B1 :* eW_B, res_B)
        V_B = (invG_B1 * vce_mat * invG_B1)[o + 2, o + 2]
        BWreg = 3 * BConst1^2 * V_B
    }

    // Bias component B2
    w_B = didhetero_kernel_eval((X :- c_eval) / h_B2, kernel)
    ind_B = (w_B :> 0)
    n_B = sum(ind_B)
    if (n_B < o_B + 3) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }

    eY_B = select(Y, ind_B)
    eX_B = select(X, ind_B)
    eW_B = select(w_B, ind_B)

    R_B2 = J(n_B, o_B + 2, .)
    for (j = 1; j <= o_B + 2; j++) {
        R_B2[., j] = (eX_B :- c_eval) :^ (j - 1)
    }
    Gram_tmp = cross(R_B2 :* sqrt(eW_B), R_B2 :* sqrt(eW_B))
    // Condition number check and Tikhonov regularization
    cond_num = cond(Gram_tmp)
    if (cond_num != . & cond_num > 1e12) {
        reg_lambda = 1e-8 * trace(Gram_tmp) / cols(Gram_tmp)
        Gram_tmp = Gram_tmp + reg_lambda * I(cols(Gram_tmp))
    }
    invG_B2 = cholinv(Gram_tmp)
    if (hasmissing(invG_B2)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
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

    if (abs(B1_out^2 + scale * R_out) < 1e-30) {
        return(.)
    }

    bw = ((rV * V_final) / (N * rB * (B1_out^2 + scale * R_out)))^r_exp

    return(bw)
}


// --- 2.4a Pre-asymptotic pilot estimates ---
// Computes pre-asymptotic bias and variance terms at each evaluation point
// using a 4-stage pilot bandwidth cascade. Used by both MSE-DPI and IMSE-DPI.
void _didhetero_bw_pilot_estimates(real colvector Y, real colvector X,
                                    real colvector eval, real scalar p,
                                    real scalar deriv,
                                    string scalar kernel,
                                    real colvector B_hat_sq,
                                    real colvector V_hat)
{
    real scalar N, R_eval, q, C_c, x_iq, x_sd
    real scalar c_bw, bw_max_r, bw_min
    real scalar bw_mp2, bw_mp3, b_mse_dpi, h_val
    real scalar V_tmp, B1_tmp, B2_tmp, R_tmp
    real scalar r_idx, j
    real colvector X_s, Y_s, dups, dupsid
    real colvector sort_idx
    real scalar range_X, x_min, x_max

    N = rows(X)
    R_eval = rows(eval)
    q = p + 1

    // Initialize outputs
    B_hat_sq = J(R_eval, 1, .)
    V_hat    = J(R_eval, 1, .)

    // Pilot bandwidth constant
    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else C_c = 2.34

    x_sd = sqrt(variance(X))
    // Interquartile range
    {
        real colvector X_sorted
        real scalar q25, q75, idx_r, j_r, h_r
        X_sorted = sort(X, 1)
        idx_r = 1 + (N - 1) * 0.25
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q25 = X_sorted[N]
        else q25 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        idx_r = 1 + (N - 1) * 0.75
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q75 = X_sorted[N]
        else q75 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        x_iq = q75 - q25
    }

    x_min = min(X)
    x_max = max(X)
    range_X = x_max - x_min

    // Sort data
    sort_idx = order(X, 1)
    X_s = X[sort_idx]
    Y_s = Y[sort_idx]

    // Compute duplicate counts and IDs
    dups = J(N, 1, 0)
    dupsid = J(N, 1, 0)
    for (j = 1; j <= N; j++) {
        dups[j] = sum(X_s :== X_s[j])
    }
    j = 1
    while (j <= N) {
        {
            real scalar k_dup
            for (k_dup = 0; k_dup < dups[j]; k_dup++) {
                dupsid[j + k_dup] = k_dup + 1
            }
            j = j + dups[j]
        }
    }

    // Four-stage pilot bandwidth cascade per evaluation point
    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        real scalar c_eval_r

        c_eval_r = eval[r_idx]

        // Maximum bandwidth
        bw_max_r = max((abs(c_eval_r - x_min), abs(c_eval_r - x_max)))

        // Pilot bandwidth constant
        c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
        c_bw = min((c_bw, bw_max_r))

        // Minimum bandwidth (bwcheck)
        {
            real colvector abs_dist
            abs_dist = sort(abs(X_s :- c_eval_r), 1)
            if (21 <= N) {
                bw_min = abs_dist[21]
            }
            else {
                bw_min = abs_dist[N]
            }
        }
        c_bw = max((c_bw, bw_min))

        // Stage 1: Pilot bandwidth for bias estimation
        bw_mp2 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 1, q + 1, q + 2,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (bw_mp2 == .) bw_mp2 = range_X
        bw_mp2 = min((bw_mp2, bw_max_r))
        bw_mp2 = max((bw_mp2, bw_min))

        // Stage 2: Pilot bandwidth for higher-order bias
        bw_mp3 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 2, q + 2, q + 3,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (bw_mp3 == .) bw_mp3 = range_X
        bw_mp3 = min((bw_mp3, bw_max_r))
        bw_mp3 = max((bw_mp3, bw_min))

        // Stage 3: Pilot bandwidth for MSE
        b_mse_dpi = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q, p + 1, q + 1,
                    c_bw, bw_mp2, bw_mp3, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (b_mse_dpi == .) b_mse_dpi = range_X
        b_mse_dpi = min((b_mse_dpi, bw_max_r))
        b_mse_dpi = max((b_mse_dpi, bw_min))

        // Stage 4: Extract variance and bias estimates
        h_val = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    p, deriv, q,
                    c_bw, b_mse_dpi, bw_mp2, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)

        // Store estimates
        if (V_tmp != . & B1_tmp != .) {
            B_hat_sq[r_idx] = B1_tmp^2
            V_hat[r_idx]    = V_tmp
        }
    }
}


// --- 2.4a2 Extended pilot estimates for even case ---
// Computes pre-asymptotic quantities for cases where (p - deriv) is even.
// Returns raw variance, bias components, and optimal bandwidth per point.
void _didhetero_bw_pilot_even(real colvector Y, real colvector X,
                                         real colvector eval, real scalar p,
                                         real scalar deriv,
                                         string scalar kernel,
                                         real colvector V_raw,
                                         real colvector B1_raw,
                                         real colvector B2_raw,
                                         real colvector h_mse_vec)
{
    real scalar N, R_eval, q, C_c, x_iq, x_sd
    real scalar c_bw, bw_max_r, bw_min
    real scalar bw_mp2, bw_mp3, b_mse_dpi, h_val
    real scalar V_tmp, B1_tmp, B2_tmp, R_tmp
    real scalar r_idx, j
    real colvector X_s, Y_s, dups, dupsid
    real colvector sort_idx
    real scalar range_X, x_min, x_max

    N = rows(X)
    R_eval = rows(eval)
    q = p + 1

    V_raw    = J(R_eval, 1, .)
    B1_raw   = J(R_eval, 1, .)
    B2_raw   = J(R_eval, 1, .)
    h_mse_vec = J(R_eval, 1, .)

    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else C_c = 2.34

    x_sd = sqrt(variance(X))
    {
        real colvector X_sorted
        real scalar q25, q75, idx_r, j_r, h_r
        X_sorted = sort(X, 1)
        idx_r = 1 + (N - 1) * 0.25
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q25 = X_sorted[N]
        else q25 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        idx_r = 1 + (N - 1) * 0.75
        j_r = floor(idx_r)
        h_r = idx_r - j_r
        if (j_r >= N) q75 = X_sorted[N]
        else q75 = (1 - h_r) * X_sorted[j_r] + h_r * X_sorted[j_r + 1]
        x_iq = q75 - q25
    }

    x_min = min(X)
    x_max = max(X)
    range_X = x_max - x_min

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
        {
            real scalar k_dup
            for (k_dup = 0; k_dup < dups[j]; k_dup++) {
                dupsid[j + k_dup] = k_dup + 1
            }
            j = j + dups[j]
        }
    }

    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        real scalar c_eval_r

        c_eval_r = eval[r_idx]
        bw_max_r = max((abs(c_eval_r - x_min), abs(c_eval_r - x_max)))

        c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
        c_bw = min((c_bw, bw_max_r))

        {
            real colvector abs_dist
            abs_dist = sort(abs(X_s :- c_eval_r), 1)
            if (21 <= N) bw_min = abs_dist[21]
            else bw_min = abs_dist[N]
        }
        c_bw = max((c_bw, bw_min))

        // Stage 1
        bw_mp2 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 1, q + 1, q + 2,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        bw_mp2 = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, 0,
                    N, range_X, q + 1, q + 1, 0)
        if (bw_mp2 == .) bw_mp2 = range_X
        bw_mp2 = min((bw_mp2, bw_max_r))
        bw_mp2 = max((bw_mp2, bw_min))

        // Stage 2
        bw_mp3 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 2, q + 2, q + 3,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        bw_mp3 = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, 0,
                    N, range_X, q + 2, q + 2, 0)
        if (bw_mp3 == .) bw_mp3 = range_X
        bw_mp3 = min((bw_mp3, bw_max_r))
        bw_mp3 = max((bw_mp3, bw_min))

        // Stage 3
        b_mse_dpi = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q, p + 1, q + 1,
                    c_bw, bw_mp2, bw_mp3, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        b_mse_dpi = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, R_tmp,
                    N, range_X, q, p + 1, 1)
        if (b_mse_dpi == .) b_mse_dpi = range_X
        b_mse_dpi = min((b_mse_dpi, bw_max_r))
        b_mse_dpi = max((b_mse_dpi, bw_min))

        // Stage 4: extract V, B1, B2 and compute h_mse_dpi
        h_val = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    p, deriv, q,
                    c_bw, b_mse_dpi, bw_mp2, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)

        // For even case: need raw V, B1, B2, and the optimize-based h_mse_dpi.
        // IMSE-DPI formula: V_h = V and B_h = (B1 + h_mse_dpi*B2)^2.
        if (V_tmp != . & B1_tmp != .) {
            V_raw[r_idx]     = V_tmp
            B1_raw[r_idx]    = B1_tmp
            B2_raw[r_idx]    = (B2_tmp != . ? B2_tmp : 0)
            h_val = _didhetero_bw_even_optimize(V_tmp, B1_tmp, B2_tmp, R_tmp,
                        N, range_X, p, deriv, 1)
            if (h_val == .) h_val = range_X
            h_val = min((h_val, bw_max_r))
            h_val = max((h_val, bw_min))
            h_mse_vec[r_idx] = h_val
        }
    }
}


// =====================================================================
// Section 3: Bandwidth undersmoothing condition check and adjustment
// Assumption 4(iii): C*n^{-1/2+eps} <= h <= C*n^{-1/9-eps}
// =====================================================================

// --- 3.0 Adjust undersmoothing: clip bandwidth to theoretical bounds ---
// Returns the (possibly adjusted) bandwidth vector.
// When do_adjust=1, bandwidths outside bounds are clamped.
// When do_adjust=0, only warnings are issued (current behavior).
//
// Returns:
//   bw_vec   - (modified in place) K x 1 adjusted bandwidths
//   real scalar indicating whether any adjustment occurred (1=yes, 0=no)
//
real scalar _didhetero_bw_adjust_undersmooth(real colvector bw_vec,
                                              real scalar n,
                                              real scalar do_adjust)
{
    real scalar K, eps, h_lower, h_upper
    real scalar id_gt, h_val, h_new, any_adjusted

    K = rows(bw_vec)
    if (K == 0) return(0)

    eps = 0.01  // Theoretical epsilon from Assumption 4(iii)
    h_lower = n^(-0.5 + eps)        // n^{-0.49}: lower bound
    h_upper = n^(-1.0/9.0 - eps)    // n^{-1/9 - 0.01}: upper bound
    any_adjusted = 0

    for (id_gt = 1; id_gt <= K; id_gt++) {
        h_val = bw_vec[id_gt]
        if (h_val >= . | h_val <= 0) continue

        if (h_val > h_upper) {
            if (do_adjust) {
                h_new = h_upper
                printf("{txt}Note: bandwidth adjusted from %9.6f to %9.6f ", h_val, h_new)
                printf("(undersmoothing upper bound, Assumption 4(iii))\n")
                printf("{txt}  Adjustment ratio: %5.1f%%%% reduction\n",
                       (h_val - h_new)/h_val * 100)
                bw_vec[id_gt] = h_new
                any_adjusted = 1
            }
        }

        if (h_val < h_lower) {
            if (do_adjust) {
                h_new = h_lower
                printf("{txt}Note: bandwidth adjusted from %9.6f to %9.6f ", h_val, h_new)
                printf("(undersmoothing lower bound, Assumption 4(iii))\n")
                printf("{txt}  Adjustment ratio: %5.1f%%%% increase\n",
                       (h_new - h_val)/h_val * 100)
                bw_vec[id_gt] = h_new
                any_adjusted = 1
            }
        }
    }

    return(any_adjusted)
}


// --- 3.1 Check undersmoothing condition for a vector of bandwidths ---
// Prints a summary warning if any bandwidths violate the theoretical bounds.
//
// Args:
//   bw_vec      - K x 1 vector of selected bandwidths (one per (g,t) pair)
//   n           - sample size
//   bwselect    - bandwidth selection method string (for context in message)
//
void _didhetero_bw_check_undersmooth(real colvector bw_vec, real scalar n,
                                      string scalar bwselect)
{
    real scalar K, eps, h_lower, h_upper
    real scalar n_below, n_above, n_nh_low
    real scalar h_min, h_max, nh_min
    real scalar id_gt, h_val, nh_eff

    K = rows(bw_vec)
    if (K == 0) return

    // Skip check for manual bandwidth (user has full control)
    // Still check — user may want the diagnostic
    // if (bwselect == "manual") return

    eps = 0.01  // Theoretical epsilon from Assumption 4(iii)
    h_lower = n^(-0.5 + eps)      // n^{-0.49}: lower bound
    h_upper = n^(-1/9 - eps)      // n^{-1/9 - 0.01}: upper bound (undersmoothing)

    n_below = 0
    n_above = 0
    n_nh_low = 0
    h_min = .
    h_max = 0
    nh_min = .

    for (id_gt = 1; id_gt <= K; id_gt++) {
        h_val = bw_vec[id_gt]
        if (h_val >= . | h_val <= 0) continue

        // Track range
        if (h_min >= . | h_val < h_min) h_min = h_val
        if (h_val > h_max) h_max = h_val

        // Effective sample size
        nh_eff = n * h_val
        if (nh_min >= . | nh_eff < nh_min) nh_min = nh_eff

        // Check bounds
        if (h_val < h_lower) n_below = n_below + 1
        if (h_val > h_upper) n_above = n_above + 1
        if (nh_eff < 30)     n_nh_low = n_nh_low + 1
    }

    // --- Output summary warnings ---

    // Check 1: Bandwidths below lower bound (effective sample too small)
    if (n_below > 0) {
        printf("{txt}Warning: %g out of %g (g,t) pairs have bandwidth below lower bound.\n",
               n_below, K)
        printf("{txt}  h < n^(-0.49) = %9.6f. Effective sample size may be too small.\n",
               h_lower)
    }

    // Check 2: Bandwidths above upper bound (undersmoothing violated)
    if (n_above > 0) {
        printf("{txt}Warning: %g out of %g (g,t) pairs have bandwidths violating undersmoothing condition.\n",
               n_above, K)
        printf("{txt}  Range of selected bandwidths: [%9.6f, %9.6f]\n", h_min, h_max)
        printf("{txt}  Theoretical upper bound for n=%g: n^(-1/9-eps) = %9.6f\n",
               n, h_upper)
        printf("{txt}  Assumption 4(iii) may be violated; analytical confidence bands\n")
        printf("{txt}  may have incorrect coverage. Consider bootstrap inference.\n")
    }

    // Check 3: Extremely low effective sample size
    if (n_nh_low > 0) {
        printf("{txt}Warning: %g out of %g (g,t) pairs have effective sample size nh < 30.\n",
               n_nh_low, K)
        printf("{txt}  Minimum nh = %9.1f. Asymptotic approximations may be unreliable.\n",
               nh_min)
    }
}


end
