mata:

// =============================================================================
// didhetero_bwselect.mata
// Bandwidth selection + density estimation wrappers
//
// Implements three bandwidth selection algorithms:
//   A. lpbwselect MSE-DPI  - for local polynomial regression
//   B. kdbwselect MSE-DPI  - for kernel density estimation
//   C. lpdensity MSE-DPI   - for density derivative estimation
//
// References:
//   Calonico, Cattaneo, Farrell (2019). JSS, 91(8): 1-33.
//   Cattaneo, Jansson, Ma (2020). JASA, 115(531): 1449-1455.
//   Cattaneo, Jansson, Ma (2022). JSS, 101(2): 1-25.
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
    _error("invalid kernel or m")
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
// Pre-asymptotic MSE-DPI bandwidth selection following
// Calonico, Cattaneo, Farrell (2019, JSS).
//
// Internal helper functions:
//   _didhetero_lpbwselect_mse   - pointwise MSE-DPI bandwidth
//   _didhetero_lprobust_bw      - core bandwidth computation
//   _didhetero_lprobust_res     - nearest-neighbor residuals
//   _didhetero_lprobust_vce     - sandwich variance estimator
// =====================================================================

// --- 2.1 _didhetero_lprobust_res: NN residuals on SORTED data ---
// Computes nearest-neighbor residuals for vce="nn"
// X must be sorted, dups/dupsid precomputed
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
        // y.J = sum(y[ind.J]) - y[pos]
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

// --- 2.2 _didhetero_lprobust_vce: sandwich variance (no cluster) ---
// Returns crossprod(c(res) * RX) = (res .* RX)' * (res .* RX)
real matrix _didhetero_lprobust_vce(real matrix RX, real colvector res)
{
    return(cross(res :* RX, res :* RX))
}

// --- Trapezoidal integration on sorted grid with missing value handling ---
// Used by IMSE1/IMSE2 bandwidth selectors (Stories 3.3/3.4)
// 
// Args:
//   grid   - R x 1 sorted (ascending) evaluation points
//   values - R x 1 function values at grid points
//
// Returns:
//   scalar: approximate integral via trapezoidal rule
//   Returns missing (.) if R < 2 or all intervals contain missing values
//
// Missing value handling:
//   Skips intervals where grid[r], grid[r+1], values[r], or values[r+1] is missing
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

// --- 2.3 _didhetero_lprobust_bw: core bandwidth computation ---
// Implements the pre-asymptotic MSE-DPI bandwidth formula for vce="nn", no cluster
//
// Args:
//   Y, X      - full data (sorted by X)
//   c_eval    - scalar evaluation point
//   o         - polynomial order for variance estimation
//   nu        - derivative order for variance extraction
//   o_B       - polynomial order for bias estimation
//   h_V       - bandwidth for variance part
//   h_B1      - bandwidth for bias B1 part
//   h_B2      - bandwidth for bias B2 part
//   scale     - bwregul parameter (0 or 1)
//   nnmatch   - number of NN matches
//   kernel    - "epa" or "gau"
//   dups      - duplicate counts (precomputed on sorted X)
//   dupsid    - duplicate IDs (precomputed on sorted X)
//
// Returns: bw (scalar bandwidth), or . on failure
// Also returns V, B1, B2, R via pointer arguments
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

    // Unscaled basis: R_V[,j] = (eX - c)^(j-1)
    R_V = J(n_V, o + 1, .)
    for (j = 1; j <= o + 1; j++) {
        R_V[., j] = (eX :- c_eval) :^ (j - 1)
    }

    // invG_V = (R_V' * diag(eW) * R_V)^{-1} via Cholesky decomposition
    invG_V = cholinv(cross(R_V :* sqrt(eW), R_V :* sqrt(eW)))
    if (hasmissing(invG_V)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }

    // NN residuals on the variance subsample
    dups_V = select(dups, ind_V)
    dupsid_V = select(dupsid, ind_V)
    res_V = _didhetero_lprobust_res(eX, eY, nnmatch, dups_V, dupsid_V)

    // V_V = (invG_V * crossprod(res_V .* RW) * invG_V)[nu+1, nu+1]
    RW = R_V :* eW
    vce_mat = _didhetero_lprobust_vce(RW, res_V)
    V_V = (invG_V * vce_mat * invG_V)[nu + 1, nu + 1]

    // ===== BIAS CONSTANTS (uses h_V data) =====
    // Hp[j] = h_V^(j-1)
    Hp = J(o + 1, 1, .)
    for (j = 1; j <= o + 1; j++) {
        Hp[j] = h_V^(j - 1)
    }

    // v1 = R_V' * diag(eW) * ((eX-c)/h_V)^(o+1)
    v1 = cross(R_V, eW :* (((eX :- c_eval) / h_V) :^ (o + 1)))
    // v2 = R_V' * diag(eW) * ((eX-c)/h_V)^(o+2)
    v2 = cross(R_V, eW :* (((eX :- c_eval) / h_V) :^ (o + 2)))

    // BConst1 = (Hp .* (invG_V * v1))[nu+1]
    BConst1 = (Hp :* (invG_V * v1))[nu + 1]
    // BConst2 = (Hp .* (invG_V * v2))[nu+1]
    BConst2 = (Hp :* (invG_V * v2))[nu + 1]

    // ===== BIAS B1 (uses h_B1, kernel weights NOT divided by h) =====
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
    invG_B1 = cholinv(cross(R_B1 :* sqrt(eW_B), R_B1 :* sqrt(eW_B)))
    if (hasmissing(invG_B1)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }
    beta_B1 = invG_B1 * cross(R_B1, eW_B :* eY_B)

    // BWreg computation (if scale > 0)
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

    // ===== BIAS B2 (uses h_B2, kernel weights NOT divided by h) =====
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
    invG_B2 = cholinv(cross(R_B2 :* sqrt(eW_B), R_B2 :* sqrt(eW_B)))
    if (hasmissing(invG_B2)) {
        V_out = .; B1_out = .; B2_out = .; R_out = 0
        return(.)
    }
    beta_B2 = invG_B2 * cross(R_B2, eW_B :* eY_B)

    // ===== ASSEMBLE =====
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


// --- 2.4a _didhetero_bw_pilot_estimates: shared pre-asymptotic pilot ---
// Computes pre-asymptotic bias^2 and variance terms at each evaluation
// point using the 4-stage pilot bandwidth cascade.
//
// Both MSE-DPI and IMSE-DPI share this function:
//   - MSE-DPI applies the pointwise formula to B_hat_sq, V_hat
//   - IMSE-DPI integrates B_hat_sq, V_hat then applies the global formula
//
// Stages 1-3 compute pilot bandwidths,
// Stage 4 extracts the final V and B1 components.
//
// Args:
//   Y        - n x 1 dependent variable
//   X        - n x 1 independent variable
//   eval     - R x 1 evaluation points
//   p        - polynomial order
//   deriv    - derivative order (nu)
//   kernel   - "epa" or "gau"
//   B_hat_sq - output: R x 1, pre-asymptotic B_hat_1^2(z_r)
//   V_hat    - output: R x 1, pre-asymptotic V_hat(z_r)
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

    // --- Initialize outputs ---
    B_hat_sq = J(R_eval, 1, .)
    V_hat    = J(R_eval, 1, .)

    // --- Pilot bandwidth constant ---
    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else C_c = 2.34

    x_sd = sqrt(variance(X))
    // IQR: Hyndman & Fan type-7 quantile interpolation
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

    // --- Sort data by X for NN residuals ---
    sort_idx = order(X, 1)
    X_s = X[sort_idx]
    Y_s = Y[sort_idx]

    // --- Compute dups and dupsid on sorted data ---
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

    // --- Per eval point: 4-stage cascade ---
    for (r_idx = 1; r_idx <= R_eval; r_idx++) {
        real scalar c_eval_r

        c_eval_r = eval[r_idx]

        // bw.max = max(|eval - x.min|, |eval - x.max|)
        bw_max_r = max((abs(c_eval_r - x_min), abs(c_eval_r - x_max)))

        // c.bw = C.c * min(sd(x), IQR/1.349) * N^(-1/5)
        c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
        c_bw = min((c_bw, bw_max_r))

        // bwcheck: bw.min = sort(|x - eval|)[21]
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

        // === Stage 1: bw.mp2 ===
        bw_mp2 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 1, q + 1, q + 2,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (bw_mp2 == .) bw_mp2 = range_X
        bw_mp2 = min((bw_mp2, bw_max_r))
        bw_mp2 = max((bw_mp2, bw_min))

        // === Stage 2: bw.mp3 ===
        bw_mp3 = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q + 2, q + 2, q + 3,
                    c_bw, range_X, range_X, 0,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (bw_mp3 == .) bw_mp3 = range_X
        bw_mp3 = min((bw_mp3, bw_max_r))
        bw_mp3 = max((bw_mp3, bw_min))

        // === Stage 3: b.mse.dpi ===
        b_mse_dpi = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    q, p + 1, q + 1,
                    c_bw, bw_mp2, bw_mp3, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)
        if (b_mse_dpi == .) b_mse_dpi = range_X
        b_mse_dpi = min((b_mse_dpi, bw_max_r))
        b_mse_dpi = max((b_mse_dpi, bw_min))

        // === Stage 4: extract pre-asymptotic quantities ===
        h_val = _didhetero_lprobust_bw(Y_s, X_s, c_eval_r,
                    p, deriv, q,
                    c_bw, b_mse_dpi, bw_mp2, 1,
                    3, kernel, dups, dupsid,
                    V_tmp, B1_tmp, B2_tmp, R_tmp)

        // Store pre-asymptotic terms (missing if cascade failed)
        if (V_tmp != . & B1_tmp != .) {
            B_hat_sq[r_idx] = B1_tmp^2
            V_hat[r_idx]    = V_tmp
        }
    }
}


// --- 2.4a2 _didhetero_bw_pilot_even: extended pilot for even case ---
// Same 4-stage cascade as _didhetero_bw_pilot_estimates, but also returns
// the raw B1, B2, and h_mse_dpi per point needed for the even branch of
// lpbwselect.imse.dpi.
//
// For even case (p-deriv)%%2==0, the IMSE-DPI formula uses:
//   V_h = V  (raw, no rV scaling)
//   B_h = (B1 + h_mse_dpi * B2)^2
// So we need B1, B2, h_mse per point.
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


// --- 2.4b _dh_lpbw_mse_details: exact pointwise MSE-DPI details ---
// Computes the pointwise MSE-DPI bandwidth for one eval point.
// Returns the pointwise h_mse_dpi and the V/B objects that IMSE-DPI averages.
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

    // --- Pilot bandwidth constant ---
    if (kernel == "epa") C_c = 2.34
    else if (kernel == "gau") C_c = 1.06
    else C_c = 2.34

    x_sd = sqrt(variance(X))
    // IQR: Hyndman & Fan type-7 quantile interpolation
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

    // --- Sort data by X for NN residuals ---
    sort_idx = order(X, 1)
    X_s = X[sort_idx]
    Y_s = Y[sort_idx]

    // --- Compute dups and dupsid on sorted data ---
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

    // bw.max = max(|eval - x.min|, |eval - x.max|)
    bw_max_r = max((abs(eval_pt - x_min), abs(eval_pt - x_max)))

    // c.bw = C.c * min(sd(x), IQR/1.349) * N^(-1/5)
    c_bw = C_c * min((x_sd, x_iq / 1.349)) * N^(-1 / 5)
    c_bw = min((c_bw, bw_max_r))

    // bwcheck: bw.min = sort(|x - eval|)[21]
    {
        real colvector abs_dist
        abs_dist = sort(abs(X_s :- eval_pt), 1)
        if (21 <= N) bw_min = abs_dist[21]
        else bw_min = abs_dist[N]
    }
    c_bw = max((c_bw, bw_min))

    // === Stage 1: bw.mp2 ===
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

    // === Stage 2: bw.mp3 ===
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

    // === Stage 3: b.mse.dpi ===
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

    // === Stage 4: h_mse_dpi (final bandwidth h) ===
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


// --- 2.4c _didhetero_lpbwselect_mse: main MSE-DPI bandwidth selector ---
// Computes pointwise MSE-DPI bandwidth for each evaluation point.
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


// --- 2.5 _didhetero_lpbwselect_imse: IMSE-DPI common bandwidth selector ---
// Returns a SCALAR bandwidth that minimizes integrated MSE over eval grid.
// Implements IMSE-DPI bandwidth selection (both even and odd branches)
// following Calonico, Cattaneo, Farrell (2019, JSS).
//
// Algorithm (IMSE-DPI bandwidth selection):
//   even = mod(p - deriv, 2) == 0
//
//   ODD branch (even==FALSE):
//     V_h = rV * V,  B_h = rB * B1^2
//     h = (mean(V_h) / (N * mean(B_h)))^(1/(2p+3))
//
//   EVEN branch (even==TRUE, interior==FALSE):
//     V_h = V (raw),  B_h = (B1 + h_mse*B2)^2
//     h = optimize(|H^(2p+2-2*deriv)*mean(B_h) + mean(V_h)/(N*H^(1+2*deriv))|)
//
//   Note: bwcheck/bwregul are NOT applied to the final IMSE result.
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

    // IMSE-DPI averages pilot quantities over a 30-point equally
    // spaced grid on [min(X), max(X)], ignoring the supplied eval.
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


// =====================================================================
// Section 3: kdbwselect MSE-DPI bandwidth for kernel density estimation
// ROT-based bandwidth for kernel density estimation (v=0, p=2)
// Reference: Calonico, Cattaneo, Farrell (2019, JSS)
// =====================================================================

real colvector _didhetero_kdbwselect_mse(real colvector X,
                                         real colvector eval,
                                         string scalar kernel)
{
    real scalar n, ng, i, j, sd_X, h_rot, b_rot, p, deriv
    real scalar f_b, f_h, B_i, V_i, h_mse_i, bw_min
    real colvector h_mse, u_b, K_eq, u_h, K_h
    real colvector sorted_dist

    n = rows(X)
    ng = rows(eval)
    p = 2
    deriv = 0
    sd_X = sqrt(variance(X))

    // ROT pilot bandwidths
    // C_h=2.34 for epa, C_b=3.49 for epa
    if (kernel == "epa") {
        h_rot = sd_X * 2.34 * n^(-1 / (1 + 2 * p))
        b_rot = sd_X * 3.49 * n^(-1 / (1 + 2 * (p + 2) + 2 * p))
    }
    else {
        h_rot = sd_X * 1.06 * n^(-1 / (1 + 2 * p))
        b_rot = sd_X * 1.00 * n^(-1 / (1 + 2 * (p + 2) + 2 * p))
    }

    h_mse = J(ng, 1, .)

    for (i = 1; i <= ng; i++) {
        // bwcheck: minimum bandwidth
        sorted_dist = sort(abs(X :- eval[i]), 1)
        if (21 <= n) {
            bw_min = sorted_dist[21]
        }
        else {
            bw_min = sorted_dist[n]
        }

        // Bias pilot: equivalent kernel for (v=p+2=4, r=p=2)
        // K_eq(u) = (105/16)*(6u^2 - 5u^4 - 1) for |u|<=1
        u_b = (X :- eval[i]) / b_rot
        K_eq = J(n, 1, 0)
        for (j = 1; j <= n; j++) {
            if (abs(u_b[j]) <= 1) {
                K_eq[j] = (105 / 16) * (6 * u_b[j]^2 - 5 * u_b[j]^4 - 1)
            }
        }
        f_b = mean(K_eq) / b_rot^(1 + p)

        // Density pilot: standard kernel (v=p=2, r=deriv=0)
        // For epa: K(u) = 0.75*(1-u^2) for |u|<=1
        u_h = (X :- eval[i]) / h_rot
        K_h = didhetero_kernel_eval(u_h, kernel)
        f_h = mean(K_h) / h_rot

        // k_v and R_v for epa, v=2, r=0: k_v=0.1, R_v=0.6
        // B = f_b * k_v, V = f_h * R_v
        B_i = f_b * 0.1
        V_i = f_h * 0.6

        // MSE formula: h_mse = ((1+2*r)*V/(2*v*N*B^2))^(1/(1+2*v+2*r))
        // with v=p=2, r=deriv=0: ((1+0)*V/(2*2*N*B^2))^(1/(1+4+0)) = (V/(4NB^2))^(1/5)
        if (abs(B_i) > 1e-20 & V_i > 0) {
            h_mse_i = ((1 + 2 * deriv) * V_i / (2 * p * n * B_i^2))^(1 / (1 + 2 * p + 2 * deriv))
        }
        else {
            h_mse_i = .
        }

        // bwcheck
        if (h_mse_i < . & h_mse_i > 0) {
            if (h_mse_i < bw_min) h_mse_i = bw_min
        }
        else {
            h_mse_i = bw_min
        }

        h_mse[i] = h_mse_i
    }

    return(h_mse)
}

// didhetero_kdrobust: density estimation using kdbwselect + simple KDE
// Implements kernel density estimation with MSE-DPI bandwidth selection.
real colvector didhetero_kdrobust(real colvector X, real colvector eval,
                                  string scalar kernel)
{
    real colvector h_bw
    h_bw = _didhetero_kdbwselect_mse(X, eval, kernel)
    return(didhetero_kde(X, eval, kernel, h_bw))
}


// =====================================================================
// Section 4: Local polynomial density derivative estimation
// Implements bw_IROT, bw_MSE, and density estimation functions
// Reference: Cattaneo, Jansson, Ma (2022, JSS)
// =====================================================================

// --- 4.1 Analytical kernel matrix generators ---

// Sgenerate: S[i,j] = integral(u^(i+j-2) * K(u), -1, 1)
real matrix _didhetero_lpdensity_Sgenerate(real scalar p, string scalar kernel)
{
    real matrix S
    real scalar i, j, a
    S = J(p + 1, p + 1, 0)
    for (i = 1; i <= p + 1; i++) {
        for (j = 1; j <= p + 1; j++) {
            a = i + j - 2
            if (mod(a, 2) == 1) {
                S[i, j] = 0
            }
            else {
                S[i, j] = _didhetero_kernel_int(a, 1, kernel)
            }
        }
    }
    return(S)
}

// Tgenerate: T[i,j] = integral(u^(i+j-2) * K(u)^2, -1, 1)
real matrix _didhetero_lpdensity_Tgenerate(real scalar p, string scalar kernel)
{
    real matrix T
    real scalar i, j, a
    T = J(p + 1, p + 1, 0)
    for (i = 1; i <= p + 1; i++) {
        for (j = 1; j <= p + 1; j++) {
            a = i + j - 2
            if (mod(a, 2) == 1) {
                T[i, j] = 0
            }
            else {
                T[i, j] = _didhetero_kernel_int(a, 2, kernel)
            }
        }
    }
    return(T)
}

// Cgenerate: C[i] = integral(u^(i+k-1) * K(u), -1, 1)
real colvector _didhetero_lpdensity_Cgenerate(real scalar k, real scalar p,
                                               string scalar kernel)
{
    real colvector C
    real scalar i, a
    C = J(p + 1, 1, 0)
    for (i = 1; i <= p + 1; i++) {
        a = i + k - 1
        if (mod(a, 2) == 1) {
            C[i] = 0
        }
        else {
            C[i] = _didhetero_kernel_int(a, 1, kernel)
        }
    }
    return(C)
}

// Ggenerate: G[i,j] = double integral for influence function variance
// G[i,j] = int_-1^1 int_-1^y x^i * y^(j-1) * K(x)*K(y) dx dy
//         + int_-1^1 int_y^1 x^(i-1) * y^j * K(x)*K(y) dx dy
// Uses midpoint rule with M points per dimension
real matrix _didhetero_lpdensity_Ggenerate(real scalar p, string scalar kernel)
{
    real matrix G
    real scalar i, j, M, dy, dx, y_m, x_m, Ky, Kx
    real scalar sum1, sum2, m1, m2

    M = 400
    G = J(p + 1, p + 1, 0)

    dy = 2.0 / M

    for (i = 1; i <= p + 1; i++) {
        for (j = 1; j <= p + 1; j++) {
            sum1 = 0
            sum2 = 0
            for (m2 = 1; m2 <= M; m2++) {
                y_m = -1 + (m2 - 0.5) * dy
                if (kernel == "epa") {
                    Ky = 0.75 * (1 - y_m^2)
                }
                else {
                    Ky = normalden(y_m)
                }
                if (abs(Ky) < 1e-30) continue

                // First integral: x from -1 to y_m
                if (y_m > -1) {
                    dx = (y_m - (-1)) / M
                    for (m1 = 1; m1 <= M; m1++) {
                        x_m = -1 + (m1 - 0.5) * dx
                        if (kernel == "epa") {
                            Kx = 0.75 * (1 - x_m^2)
                        }
                        else {
                            Kx = normalden(x_m)
                        }
                        sum1 = sum1 + x_m^i * y_m^(j - 1) * Kx * Ky * dx * dy
                    }
                }

                // Second integral: x from y_m to 1
                if (y_m < 1) {
                    dx = (1 - y_m) / M
                    for (m1 = 1; m1 <= M; m1++) {
                        x_m = y_m + (m1 - 0.5) * dx
                        if (kernel == "epa") {
                            Kx = 0.75 * (1 - x_m^2)
                        }
                        else {
                            Kx = normalden(x_m)
                        }
                        sum2 = sum2 + x_m^(i - 1) * y_m^j * Kx * Ky * dx * dy
                    }
                }
            }
            G[i, j] = sum1 + sum2
        }
    }
    return(G)
}


// --- 4.2 Normal PDF derivative via Hermite polynomials ---
// Computes d^v/dx^v [phi(x; mu, sd)] where phi is the normal PDF
// Uses the relation: phi^(v)(x) = (-1)^v / sd^v * He_v(z) * phi(z) / sd
// where z = (x-mu)/sd and He_v is the probabilist's Hermite polynomial
// He_0(z) = 1, He_1(z) = z, He_n(z) = z*He_{n-1}(z) - (n-1)*He_{n-2}(z)
real scalar _didhetero_normal_pdf_deriv(real scalar x, real scalar mu,
                                        real scalar sd_val, real scalar v)
{
    real scalar z, phi_z, He_prev, He_curr, He_next, k

    z = (x - mu) / sd_val
    phi_z = normalden(z) / sd_val

    if (v == 0) return(phi_z)

    // Hermite polynomial recurrence
    He_prev = 1
    He_curr = z
    if (v == 1) return((-1)^1 / sd_val * He_curr * phi_z)

    for (k = 2; k <= v; k++) {
        He_next = z * He_curr - (k - 1) * He_prev
        He_prev = He_curr
        He_curr = He_next
    }
    return((-1)^v / sd_val^v * He_curr * phi_z)
}

// --- 4.3 bw_IROT: Integrated Rule-of-Thumb bandwidth ---
// Returns a SCALAR bandwidth (integrated over grid points)
// Following Cattaneo, Jansson, Ma (2022, JSS)
real scalar _didhetero_lpdensity_bw_IROT(real colvector data,
                                          real colvector grid,
                                          real scalar p, real scalar v,
                                          string scalar kernel,
                                          real scalar nLocalMin,
                                          real scalar nUniqueMin)
{
    real scalar n, ng, nUnique, center_temp, scale_temp
    real scalar mean_hat, sd_hat, j, k
    real scalar phi_val, phi_deriv_val, phi_p, phi_p1
    real scalar bias1, bias2, h_opt, range_data, h_max
    real matrix bias_dgp, S, C1, C2, G, S2
    real colvector sd_dgp, dataUnique, sorted_abs, data_local, grid_local
    real matrix S_inv
    real scalar const_bias1, const_bias2, const_var
    real scalar sum_bias_sq, sum_var_sq, a_lo, a_hi
    real scalar a_m1, a_m2, f_m1, f_m2, gr

    // Mata passes arguments by reference. Work on local copies so repeated
    // pilot calls inside bw_MSE do not keep re-standardizing the caller's data.
    data_local = data
    grid_local = grid

    n = rows(data_local)
    ng = rows(grid_local)

    // Get unique values
    dataUnique = uniqrows(data_local)
    nUnique = rows(dataUnique)

    // Standardize
    center_temp = mean(data_local)
    scale_temp = sqrt(variance(data_local))
    data_local = (data_local :- center_temp) / scale_temp
    dataUnique = (dataUnique :- center_temp) / scale_temp
    grid_local = (grid_local :- center_temp) / scale_temp

    // Normal reference model.
    // The IROT rule computes the post-standardization scale with denominator n,
    // not n-1, so use the empirical second moment rather than variance().
    mean_hat = mean(data_local)
    sd_hat = sqrt(mean((data_local :- mean_hat) :^ 2))

    // Kernel matrices
    S  = _didhetero_lpdensity_Sgenerate(p, kernel)
    C1 = _didhetero_lpdensity_Cgenerate(p + 1, p, kernel)
    C2 = _didhetero_lpdensity_Cgenerate(p + 2, p, kernel)
    S2 = _didhetero_lpdensity_Tgenerate(p, kernel)
    G  = _didhetero_lpdensity_Ggenerate(p, kernel)

    S_inv = cholinv(S)

    // Kernel constants for bias and variance
    const_bias1 = (S_inv * C1)[v + 1]
    const_bias2 = (S_inv * C2)[v + 1]

    // DGP constants per grid point
    bias_dgp = J(ng, 2, .)
    sd_dgp = J(ng, 1, .)

    for (j = 1; j <= ng; j++) {
        // phi(grid[j]; mean_hat, sd_hat) = normal PDF at grid point
        phi_val = normalden(grid_local[j], mean_hat, sd_hat)

        // phi'(grid[j]) = first derivative of normal PDF
        phi_deriv_val = _didhetero_normal_pdf_deriv(grid_local[j], mean_hat, sd_hat, 1)

        // The IROT source differentiates temp_3/temp_4 only p and p+1 times.
        // This follows the pre-asymptotic pilot formula for numerical equivalence.
        phi_p = _didhetero_normal_pdf_deriv(grid_local[j], mean_hat, sd_hat, p)
        phi_p1 = _didhetero_normal_pdf_deriv(grid_local[j], mean_hat, sd_hat, p + 1)

        // bias_dgp[j,1] = temp_3 / (p+1)! * v!
        bias_dgp[j, 1] = phi_p / factorial(p + 1) * factorial(v)
        // bias_dgp[j,2] = temp_4 / (p+2)! * v! + bias_dgp[j,1] * phi'/phi
        bias_dgp[j, 2] = phi_p1 / factorial(p + 2) * factorial(v)
        if (abs(phi_val) > 1e-30) {
            bias_dgp[j, 2] = bias_dgp[j, 2] + bias_dgp[j, 1] * phi_deriv_val / phi_val
        }

        // Apply kernel constants
        bias_dgp[j, 1] = bias_dgp[j, 1] * const_bias1
        bias_dgp[j, 2] = bias_dgp[j, 2] * const_bias2

        // Variance DGP
        if (v > 0) {
            sd_dgp[j] = factorial(v) * sqrt(phi_val / n)
            sd_dgp[j] = sd_dgp[j] * sqrt(abs((S_inv * G * S_inv)[v + 1, v + 1]))
        }
        else {
            // v==0 case (not used in our application but included for completeness)
            phi_val = max((phi_val, 1e-30))
            sd_dgp[j] = sqrt(normal(grid_local[j] * sd_hat + mean_hat) * (1 - normal(grid_local[j] * sd_hat + mean_hat)) / phi_val / (0.5 * n^2))
            sd_dgp[j] = sd_dgp[j] * sqrt(abs((S_inv * S2 * S_inv)[v + 1, v + 1]))
        }
    }

    // Golden section search to minimize integrated MSE
    // For v > 0: f(a) = a^(2p+2-2v) * sum(bias1 + a*bias2)^2 + sum(sd^2) / a^(2v-1)
    range_data = max(data_local) - min(data_local)
    gr = (sqrt(5) + 1) / 2

    a_lo = 1e-10
    a_hi = range_data

    for (k = 1; k <= 100; k++) {
        a_m1 = a_hi - (a_hi - a_lo) / gr
        a_m2 = a_lo + (a_hi - a_lo) / gr

        if (v > 0) {
            sum_bias_sq = 0
            sum_var_sq = 0
            for (j = 1; j <= ng; j++) {
                sum_bias_sq = sum_bias_sq + (bias_dgp[j, 1] + a_m1 * bias_dgp[j, 2])^2
                sum_var_sq = sum_var_sq + sd_dgp[j]^2
            }
            f_m1 = a_m1^(2 * p + 2 - 2 * v) * sum_bias_sq + sum_var_sq / a_m1^(2 * v - 1)

            sum_bias_sq = 0
            for (j = 1; j <= ng; j++) {
                sum_bias_sq = sum_bias_sq + (bias_dgp[j, 1] + a_m2 * bias_dgp[j, 2])^2
            }
            f_m2 = a_m2^(2 * p + 2 - 2 * v) * sum_bias_sq + sum_var_sq / a_m2^(2 * v - 1)
        }
        else {
            sum_bias_sq = 0
            sum_var_sq = 0
            for (j = 1; j <= ng; j++) {
                sum_bias_sq = sum_bias_sq + (bias_dgp[j, 1] + a_m1 * bias_dgp[j, 2])^2
                sum_var_sq = sum_var_sq + sd_dgp[j]^2
            }
            f_m1 = a_m1^(2 * p + 2) * sum_bias_sq + sum_var_sq / a_m1

            sum_bias_sq = 0
            for (j = 1; j <= ng; j++) {
                sum_bias_sq = sum_bias_sq + (bias_dgp[j, 1] + a_m2 * bias_dgp[j, 2])^2
            }
            f_m2 = a_m2^(2 * p + 2) * sum_bias_sq + sum_var_sq / a_m2
        }

        if (f_m1 < f_m2) {
            a_hi = a_m2
        }
        else {
            a_lo = a_m1
        }

        if (abs(a_hi - a_lo) < 1e-12) break
    }

    h_opt = (a_lo + a_hi) / 2

    // Handle NA
    if (h_opt == . | h_opt <= 0) {
        h_opt = 0
        for (j = 1; j <= ng; j++) {
                sorted_abs = sort(abs(data_local :- grid_local[j]), 1)
            if (nLocalMin <= n) {
                if (sorted_abs[min((n, max((nLocalMin, 20 + p + 1))))] > h_opt) {
                    h_opt = sorted_abs[min((n, max((nLocalMin, 20 + p + 1))))]
                }
            }
        }
    }

    // Regularize
    for (j = 1; j <= ng; j++) {
        sorted_abs = sort(abs(data_local :- grid_local[j]), 1)
        if (nLocalMin > 0 & nLocalMin <= n) {
            if (sorted_abs[min((n, nLocalMin))] > h_opt) {
                h_opt = sorted_abs[min((n, nLocalMin))]
            }
        }
        if (nUniqueMin > 0 & nUniqueMin <= nUnique) {
            sorted_abs = sort(abs(dataUnique :- grid_local[j]), 1)
            if (sorted_abs[min((nUnique, nUniqueMin))] > h_opt) {
                h_opt = sorted_abs[min((nUnique, nUniqueMin))]
            }
        }
    }

    // Cap at max range
    h_max = max((abs(max(dataUnique) - min(grid_local)), abs(min(dataUnique) - max(grid_local))))
    if (h_opt > h_max) h_opt = h_max

    // Rescale back
    h_opt = h_opt * scale_temp

    return(h_opt)
}


// --- 4.4 bw_MSE: MSE-optimal pointwise bandwidth ---
// Returns ng x 1 vector of bandwidths
// Following Cattaneo, Jansson, Ma (2022, JSS)
real colvector _didhetero_lpdensity_bw_MSE(real colvector data_in,
                                            real colvector grid_in,
                                            real scalar p, real scalar v,
                                            string scalar kernel)
{
    real scalar n, ng, nUnique, center_temp, scale_temp
    real scalar j, jj, k, i_idx
    real scalar h1, hp1, hp2
    real colvector data, grid, Fn, dataUnique
    real colvector freqUnique, indexUnique_vec
    real colvector h
    real matrix dgp_hat, const_hat
    real scalar nLocalMin_h1, nLocalMin_hp1, nLocalMin_hp2
    real scalar nUniqueMin_h1, nUniqueMin_hp1, nUniqueMin_hp2
    // Working variables for per-grid-point computation
    real colvector index_temp_vec, Xh_temp, Kh_temp, Y_temp
    real matrix Xh_p_temp
    real matrix S_hat, S_hat_inv, G_hat, G_col
    real colvector C_p_hat, C_p1_hat
    real scalar temp_val
    real matrix beta_temp
    // Golden section variables
    real scalar a_lo, a_hi, a_m1, a_m2, f_m1, f_m2, gr
    real scalar sum_bias_sq, sum_var_sq
    real colvector sorted_abs, sorted_abs_u
    real scalar bw_min_local, bw_min_unique
    real scalar range_data, h_max_j
    // Influence function variables
    real colvector F_Xh_p_Kh_row
    real matrix Xh_p_Kh_temp, cumsum_rev
    real scalar n_eff
    real matrix G_full

    // --- Sort data ---
    data = sort(data_in, 1)
    grid = grid_in  // don't sort grid, keep original order
    n = rows(data)
    ng = rows(grid)

    // --- Unique values with frequencies ---
    dataUnique = uniqrows(data)
    nUnique = rows(dataUnique)

    // Build frequency and index vectors for mass points
    freqUnique = J(nUnique, 1, 0)
    indexUnique_vec = J(nUnique, 1, 0)
    k = 1
    for (i_idx = 1; i_idx <= n; i_idx++) {
        if (i_idx == n) {
            freqUnique[k] = freqUnique[k] + 1
            indexUnique_vec[k] = i_idx
        }
        else if (data[i_idx] != data[i_idx + 1]) {
            freqUnique[k] = freqUnique[k] + 1
            indexUnique_vec[k] = i_idx
            k = k + 1
        }
        else {
            freqUnique[k] = freqUnique[k] + 1
        }
    }

    // --- Standardize ---
    center_temp = mean(data)
    scale_temp = sqrt(variance(data))
    data = (data :- center_temp) / scale_temp
    dataUnique = (dataUnique :- center_temp) / scale_temp
    grid = (grid :- center_temp) / scale_temp

    // --- ECDF with mass points ---
    // Fn = cumulative proportion, with ties getting the same value
    Fn = (1::n) / n
    // Apply mass point correction: all tied values get the max rank proportion
    for (k = 1; k <= nUnique; k++) {
        if (freqUnique[k] > 1) {
            temp_val = Fn[indexUnique_vec[k]]
            if (k == 1) {
                for (i_idx = 1; i_idx <= indexUnique_vec[k]; i_idx++) {
                    Fn[i_idx] = temp_val
                }
            }
            else {
                for (i_idx = indexUnique_vec[k - 1] + 1; i_idx <= indexUnique_vec[k]; i_idx++) {
                    Fn[i_idx] = temp_val
                }
            }
        }
    }

    // --- IROT pilot bandwidths ---
    nLocalMin_h1 = 20 + 2 + 1
    nUniqueMin_h1 = 20 + 2 + 1
    h1 = _didhetero_lpdensity_bw_IROT(data, grid, 2, 1, kernel, nLocalMin_h1, nUniqueMin_h1)
    // h1 is already in standardized scale (bw_IROT standardizes internally then rescales).
    // bw_MSE standardizes data first, then passes the standardized data to bw_IROT
    // which standardizes AGAIN internally. This double-standardization is intentional
    // and matches the density estimation convention.

    nLocalMin_hp1 = 20 + p + 2 + 1
    nUniqueMin_hp1 = 20 + p + 2 + 1
    hp1 = _didhetero_lpdensity_bw_IROT(data, grid, p + 2, p + 1, kernel, nLocalMin_hp1, nUniqueMin_hp1)

    nLocalMin_hp2 = 20 + p + 3 + 1
    nUniqueMin_hp2 = 20 + p + 3 + 1
    hp2 = _didhetero_lpdensity_bw_IROT(data, grid, p + 3, p + 2, kernel, nLocalMin_hp2, nUniqueMin_hp2)

    // --- Per grid point: estimate dgp and constants ---
    dgp_hat = J(ng, 2, .)
    const_hat = J(ng, 3, .)
    h = J(ng, 1, .)
    range_data = max(data) - min(data)
    gr = (sqrt(5) + 1) / 2

    for (j = 1; j <= ng; j++) {

        // === Estimate F_{p+2} via scaled-basis CDF-LPR with hp2 ===
        index_temp_vec = J(n, 1, 0)
        n_eff = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (abs(data[i_idx] - grid[j]) <= hp2) {
                index_temp_vec[i_idx] = 1
                n_eff = n_eff + 1
            }
        }
        if (n_eff < p + 4) continue

        Xh_temp = J(n_eff, 1, .)
        Y_temp = J(n_eff, 1, .)
        k = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (index_temp_vec[i_idx] == 1) {
                k = k + 1
                Xh_temp[k] = (data[i_idx] - grid[j]) / hp2
                Y_temp[k] = Fn[i_idx]
            }
        }

        // Scaled polynomial basis: columns 0..(p+3)
        Xh_p_temp = J(n_eff, p + 4, .)
        for (k = 0; k <= p + 3; k++) {
            Xh_p_temp[., k + 1] = Xh_temp :^ k
        }

        // Kernel weights / hp2
        Kh_temp = 0.75 * (1 :- Xh_temp :^ 2) / hp2

        // WLS: beta = (X'KX)^{-1} X'KY
        // Result: beta[p+3] / hp2^(p+2) = F_{p+2} estimate
        beta_temp = cholinv(cross(Xh_p_temp, Kh_temp, Xh_p_temp)) * cross(Xh_p_temp, Kh_temp, Y_temp)
        if (hasmissing(beta_temp)) continue
        dgp_hat[j, 2] = beta_temp[p + 3] / hp2^(p + 2)

        // === Estimate F_{p+1} via scaled-basis CDF-LPR with hp1 ===
        index_temp_vec = J(n, 1, 0)
        n_eff = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (abs(data[i_idx] - grid[j]) <= hp1) {
                index_temp_vec[i_idx] = 1
                n_eff = n_eff + 1
            }
        }
        if (n_eff < p + 3) continue

        Xh_temp = J(n_eff, 1, .)
        Y_temp = J(n_eff, 1, .)
        k = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (index_temp_vec[i_idx] == 1) {
                k = k + 1
                Xh_temp[k] = (data[i_idx] - grid[j]) / hp1
                Y_temp[k] = Fn[i_idx]
            }
        }

        Xh_p_temp = J(n_eff, p + 3, .)
        for (k = 0; k <= p + 2; k++) {
            Xh_p_temp[., k + 1] = Xh_temp :^ k
        }

        Kh_temp = 0.75 * (1 :- Xh_temp :^ 2) / hp1

        beta_temp = cholinv(cross(Xh_p_temp, Kh_temp, Xh_p_temp)) * cross(Xh_p_temp, Kh_temp, Y_temp)
        if (hasmissing(beta_temp)) continue
        dgp_hat[j, 1] = beta_temp[p + 2] / hp1^(p + 1)

        // === Pre-asymptotic matrices with h1 ===
        index_temp_vec = J(n, 1, 0)
        n_eff = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (abs(data[i_idx] - grid[j]) <= h1) {
                index_temp_vec[i_idx] = 1
                n_eff = n_eff + 1
            }
        }
        if (n_eff < p + 2) continue

        Xh_temp = J(n_eff, 1, .)
        k = 0
        for (i_idx = 1; i_idx <= n; i_idx++) {
            if (index_temp_vec[i_idx] == 1) {
                k = k + 1
                Xh_temp[k] = (data[i_idx] - grid[j]) / h1
            }
        }

        Kh_temp = 0.75 * (1 :- Xh_temp :^ 2) / h1

        // C_p_hat: (p+1) x 1, C_p_hat[i] = (1/n) * sum(Xh^(p+i) * Kh)
        C_p_hat = J(p + 1, 1, 0)
        for (k = 1; k <= p + 1; k++) {
            C_p_hat[k] = sum(Xh_temp :^ (p + k) :* Kh_temp) / n
        }

        // C_p1_hat: (p+1) x 1, C_p1_hat[i] = (1/n) * sum(Xh^(p+1+i) * Kh)
        C_p1_hat = J(p + 1, 1, 0)
        for (k = 1; k <= p + 1; k++) {
            C_p1_hat[k] = sum(Xh_temp :^ (p + 1 + k) :* Kh_temp) / n
        }

        // S_hat: (p+1) x (p+1)
        Xh_p_temp = J(n_eff, p + 1, .)
        for (k = 0; k <= p; k++) {
            Xh_p_temp[., k + 1] = Xh_temp :^ k
        }
        S_hat = cross(Xh_p_temp, Kh_temp, Xh_p_temp) / n
        S_hat_inv = cholinv(S_hat)
        if (hasmissing(S_hat_inv)) continue

        // G_hat via influence function approach (for v > 0)
        if (v > 0) {
            // Use full sample, influence function approach with mass points
            // Recompute Xh and Kh for ALL unique data points
            Xh_temp = (dataUnique :- grid[j]) / h1
            Xh_p_temp = J(nUnique, p + 1, .)
            for (k = 0; k <= p; k++) {
                Xh_p_temp[., k + 1] = Xh_temp :^ k
            }

            // Kh with index_temp masking
            Kh_temp = J(nUnique, 1, 0)
            for (i_idx = 1; i_idx <= nUnique; i_idx++) {
                if (abs(dataUnique[i_idx] - grid[j]) <= h1) {
                    Kh_temp[i_idx] = 0.75 * (1 - Xh_temp[i_idx]^2) / h1
                }
            }

            // Xh_p_Kh_temp = Xh_p_temp .* Kh_temp (each row scaled by Kh)
            Xh_p_Kh_temp = Xh_p_temp :* Kh_temp

            // F_Xh_p_Kh = (1/n) * Fn_unique' * Xh_p_Kh_temp : 1 x (p+1)
            // Get Fn at unique points
            Y_temp = Fn[indexUnique_vec]
            F_Xh_p_Kh_row = (Y_temp' * Xh_p_Kh_temp / n)'  // (p+1) x 1

            // Build full n x (p+1) influence matrix for G_hat
            G_full = J(n, p + 1, 0)
            for (jj = 1; jj <= p + 1; jj++) {
                cumsum_rev = J(nUnique, 1, 0)
                cumsum_rev[nUnique] = Xh_p_Kh_temp[nUnique, jj]
                for (i_idx = nUnique - 1; i_idx >= 1; i_idx--) {
                    cumsum_rev[i_idx] = cumsum_rev[i_idx + 1] + Xh_p_Kh_temp[i_idx, jj]
                }
                cumsum_rev = cumsum_rev / n

                for (i_idx = 1; i_idx <= nUnique; i_idx++) {
                    temp_val = cumsum_rev[i_idx] - F_Xh_p_Kh_row[jj]
                    if (i_idx == 1) {
                        for (k = 1; k <= indexUnique_vec[i_idx]; k++) {
                            G_full[k, jj] = temp_val
                        }
                    }
                    else {
                        for (k = indexUnique_vec[i_idx - 1] + 1; k <= indexUnique_vec[i_idx]; k++) {
                            G_full[k, jj] = temp_val
                        }
                    }
                }
            }
            G_hat = cross(G_full, G_full) / n
        }
        else {
            // v == 0: G_hat = T matrix (kernel squared)
            Xh_temp = J(n_eff, 1, .)
            k = 0
            for (i_idx = 1; i_idx <= n; i_idx++) {
                if (index_temp_vec[i_idx] == 1) {
                    k = k + 1
                    Xh_temp[k] = (data[i_idx] - grid[j]) / h1
                }
            }
            Kh_temp = 0.75 * (1 :- Xh_temp :^ 2) / h1
            Xh_p_temp = J(n_eff, p + 1, .)
            for (k = 0; k <= p; k++) {
                Xh_p_temp[., k + 1] = Xh_temp :^ k
            }
            G_hat = cross(Xh_p_temp, Kh_temp :^ 2, Xh_p_temp) / n
        }

        // === Assemble constants ===
        const_hat[j, 1] = factorial(v) * (S_hat_inv * C_p_hat)[v + 1]
        const_hat[j, 2] = factorial(v) * (S_hat_inv * C_p1_hat)[v + 1]

        if (v > 0) {
            const_hat[j, 3] = factorial(v) * sqrt(abs((S_hat_inv * G_hat * S_hat_inv)[v + 1, v + 1]) / (n * h1))
        }
        else {
            temp_val = min((max((mean(data :<= grid[j]), 1 / n)), 1 - 1 / n))
            const_hat[j, 3] = factorial(v) * sqrt(abs((S_hat_inv * G_hat * S_hat_inv)[v + 1, v + 1] / (0.5 * n^2) * h1 * temp_val * (1 - temp_val)))
        }

        // === Golden section search for optimal h[j] ===
        a_lo = 1e-10
        a_hi = range_data

        for (k = 1; k <= 100; k++) {
            a_m1 = a_hi - (a_hi - a_lo) / gr
            a_m2 = a_lo + (a_hi - a_lo) / gr

            if (v > 0) {
                f_m1 = a_m1^(2 * p + 2 - 2 * v) * (dgp_hat[j, 1] * const_hat[j, 1] + a_m1 * dgp_hat[j, 2] * const_hat[j, 2])^2 + const_hat[j, 3]^2 / a_m1^(2 * v - 1)
                f_m2 = a_m2^(2 * p + 2 - 2 * v) * (dgp_hat[j, 1] * const_hat[j, 1] + a_m2 * dgp_hat[j, 2] * const_hat[j, 2])^2 + const_hat[j, 3]^2 / a_m2^(2 * v - 1)
            }
            else {
                f_m1 = a_m1^(2 * p + 2) * (dgp_hat[j, 1] * const_hat[j, 1] + a_m1 * dgp_hat[j, 2] * const_hat[j, 2])^2 + const_hat[j, 3]^2 / a_m1
                f_m2 = a_m2^(2 * p + 2) * (dgp_hat[j, 1] * const_hat[j, 1] + a_m2 * dgp_hat[j, 2] * const_hat[j, 2])^2 + const_hat[j, 3]^2 / a_m2
            }

            if (f_m1 < f_m2) {
                a_hi = a_m2
            }
            else {
                a_lo = a_m1
            }

            if (abs(a_hi - a_lo) < 1e-12) break
        }

        h[j] = (a_lo + a_hi) / 2
    }

    // --- Post-processing: regularize ---
    for (j = 1; j <= ng; j++) {
        sorted_abs = sort(abs(data :- grid[j]), 1)
        if (h[j] == .) {
            h[j] = sorted_abs[min((n, max((20 + p + 1, 20 + p + 1))))]
        }
        // nLocalMin regularization
        bw_min_local = sorted_abs[min((n, 20 + p + 1))]
        if (h[j] < bw_min_local) h[j] = bw_min_local
        // nUniqueMin regularization
        sorted_abs_u = sort(abs(dataUnique :- grid[j]), 1)
        bw_min_unique = sorted_abs_u[min((nUnique, 20 + p + 1))]
        if (h[j] < bw_min_unique) h[j] = bw_min_unique
        // Cap at max distance
        h_max_j = max(abs(dataUnique :- grid[j]))
        if (h[j] > h_max_j) h[j] = h_max_j
    }

    // --- Rescale back ---
    h = h * scale_temp

    return(h)
}


// =====================================================================
// Section 5: IMSE1 / IMSE2 / US1 bandwidth selection
// (didhetero paper-specific algorithms)
//
// Reference: Imai, Qin, Yanagi (2025), Sections 4.2.5-4.2.6
// =====================================================================

// --- 5.1 Kernel integral constants ---
// Returns the 7 kernel integral constants used by IMSE1/IMSE2/US1:
//   I_{l,K^m} = integral(u^l * K^m(u) du)
// for (l,m) in {(2,1),(4,1),(6,1),(0,2),(2,2),(4,2),(6,2)}
//
// These are the SAME values computed by _didhetero_kernel_int but
// bundled for convenience and validated against the paper.
void _didhetero_kernel_constants(string scalar kernel,
                                  real scalar I_2_K1,
                                  real scalar I_4_K1,
                                  real scalar I_6_K1,
                                  real scalar I_0_K2,
                                  real scalar I_2_K2,
                                  real scalar I_4_K2,
                                  real scalar I_6_K2)
{
    if (kernel == "epa") {
        I_2_K1 = 0.2
        I_4_K1 = 3 / 35            // 0.08571428571...
        I_6_K1 = 1 / 21            // 0.04761904762...
        I_0_K2 = 0.6
        I_2_K2 = 3 / 35            // 0.08571428571...
        I_4_K2 = 1 / 35            // 0.02857142857...
        I_6_K2 = 1 / 77            // 0.01298701299...
    }
    else if (kernel == "gau") {
        I_2_K1 = 1
        I_4_K1 = 3
        I_6_K1 = 15
        I_0_K2 = 1 / (2 * sqrt(c("pi")))   // 0.28209479177...
        I_2_K2 = 1 / (4 * sqrt(c("pi")))   // 0.14104739589...
        I_4_K2 = 3 / (8 * sqrt(c("pi")))   // 0.21157109383...
        I_6_K2 = 15 / (16 * sqrt(c("pi"))) // 0.52892773458...
    }
    else {
        _error("unknown kernel: " + kernel)
    }
}


// --- 5.2 Shared variance estimation for IMSE1/IMSE2/US1 ---
// Estimates the conditional variance term mathcal_V at each evaluation point.
//
// Five-step process per evaluation point r:
//   1. IMSE-DPI common bandwidth (p=1, deriv=0) on 100-point support grid
//   2. LLR fitted values at all observation points (eval=Z)
//   3. Residuals: U_hat = y_r - mu_hat_0
//   4. MSE-DPI bandwidth for sigma^2 = U_hat^2 (p=1, deriv=0)
//   5. LPR estimate of sigma^2 at zeval[r], then V[r] = const_V * sigma^2 / f(z_r)
//
// Args:
//   B_g_t    - n x R matrix, each column is the influence function
//   Z        - n x 1 covariate vector
//   zeval    - R x 1 sorted evaluation points
//   kd0_Z    - R x 1 kernel density estimates at zeval
//   const_V  - scalar variance constant (C_V1 for IMSE1, C_V2 for IMSE2)
//   kernel   - "epa" or "gau"
//   mathcal_V - output: R x 1 variance estimates (modified by reference)
void _didhetero_bwselect_var_est(real matrix B_g_t,
                                             real colvector Z,
                                             real colvector zeval,
                                             real colvector kd0_Z,
                                             real scalar const_V,
                                             string scalar kernel,
                                             real colvector mathcal_V)
{
    real scalar n, R_eval, r
    real colvector y_r
    real scalar mu_B_0_bw
    real colvector mu_B_0, U_hat, U_hat_sq
    real colvector sigma2_bw_vec
    real scalar sigma2_bw, sigma2
    real colvector Z_supp

    n = rows(Z)
    R_eval = rows(zeval)

    // Initialize output
    mathcal_V = J(R_eval, 1, .)

    // Support grid for IMSE-DPI bandwidth (100 equally spaced points)
    Z_supp = rangen(min(Z), max(Z), 100)

    // === Loop over evaluation points ===
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // Skip if density is non-positive
        if (kd0_Z[r] >= . | kd0_Z[r] <= 0) continue

        // Step 1: IMSE-DPI common bandwidth (p=1, deriv=0)
        mu_B_0_bw = _didhetero_lpbwselect_imse(y_r, Z, Z_supp,
                        1, 0, kernel)

        if (mu_B_0_bw >= . | mu_B_0_bw <= 0) continue

        // Step 2: LLR at all observation points (eval=Z)
        mu_B_0 = didhetero_lpr(y_r, Z, Z, 1, 0, kernel, mu_B_0_bw)

        if (rows(mu_B_0) != n) continue

        // Step 3: Residuals
        U_hat = y_r - mu_B_0

        // Step 4: MSE-DPI bandwidth for sigma^2 (p=1, deriv=0)
        U_hat_sq = U_hat :^ 2
        sigma2_bw_vec = _didhetero_lpbwselect_mse(U_hat_sq, Z,
                            zeval[r..r], 1, 0, kernel)
        sigma2_bw = sigma2_bw_vec[1]

        if (sigma2_bw >= . | sigma2_bw <= 0) continue

        // Step 5: Estimate conditional variance
        {
            real colvector sigma2_vec
            sigma2_vec = didhetero_lpr(U_hat_sq, Z, zeval[r..r],
                             1, 0, kernel, sigma2_bw)
            sigma2 = sigma2_vec[1]
        }

        if (sigma2 >= . | sigma2 <= 0) continue

        // mathcal_V[r] = const_V * sigma^2 / f(z_r)
        mathcal_V[r] = const_V * sigma2 / kd0_Z[r]
    }
}


// --- 5.3a IMSE1/US1 shared internal implementation ---
// Parameterized IMSE1 bandwidth selection with configurable n exponent.
// IMSE1 uses n^(-1/5), US1 uses n^(-2/7). All other logic is identical.
//
// Algorithm (Paper Section 4.2.5):
//   For each evaluation point r = 1..R:
//     Bias:  MSE-DPI bw for 2nd deriv → LPR 2nd deriv → B[r] = mu2 * I_2_K1/2
//     Var:   Shared variance estimation via _didhetero_bwselect_var_est
//   Then: h = (int_V / (4 * int_B^2))^(1/5) * n^(n_exponent)
//
// Args:
//   B_g_t      - n x R matrix, each column is the influence function at one eval point
//   Z          - n x 1 covariate vector
//   zeval      - R x 1 sorted evaluation points
//   kd0_Z      - R x 1 kernel density estimates at zeval
//   kernel     - "epa" or "gau"
//   n_exponent - exponent for n: -1/5 for IMSE1, -2/7 for US1
//
// Returns: scalar bandwidth (or . on failure)
real scalar _dh_bwselect_imse1_internal(real matrix B_g_t,
                                      real colvector Z,
                                      real colvector zeval,
                                      real colvector kd0_Z,
                                      string scalar kernel,
                                      real scalar n_exponent)
{
    real scalar n, R_eval, r
    real scalar I_2_K1, I_4_K1, I_6_K1, I_0_K2, I_2_K2, I_4_K2, I_6_K2
    real scalar const_V1
    real colvector mathcal_B, mathcal_V
    real colvector y_r, mu_B_2_bw_vec
    real scalar mu_B_2_bw, mu_B_2
    real scalar int_bias, int_var, h_opt, h_max

    n = rows(Z)
    R_eval = rows(zeval)

    // Edge case: need at least 2 eval points for trapz
    if (R_eval < 2) return(.)

    // Get kernel constants
    _didhetero_kernel_constants(kernel, I_2_K1, I_4_K1, I_6_K1,
                                I_0_K2, I_2_K2, I_4_K2, I_6_K2)
    const_V1 = I_0_K2

    // Initialize bias vector
    mathcal_B = J(R_eval, 1, .)

    // === Bias estimation loop ===
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // --- Bias estimation ---
        // MSE-DPI bandwidth for 2nd derivative (p=3, deriv=2)
        mu_B_2_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 3, 2, kernel)
        mu_B_2_bw = mu_B_2_bw_vec[1]

        if (mu_B_2_bw >= . | mu_B_2_bw <= 0) continue

        // Estimate 2nd derivative via LPR (p=3, deriv=2)
        {
            real colvector mu_B_2_vec
            mu_B_2_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             3, 2, kernel, mu_B_2_bw)
            mu_B_2 = mu_B_2_vec[1]
        }

        if (mu_B_2 >= .) continue

        // mathcal_B[r] = mu^(2) * I_{2,K} / 2
        mathcal_B[r] = mu_B_2 * I_2_K1 / 2
    }

    // === Variance estimation via shared function ===
    _didhetero_bwselect_var_est(B_g_t, Z, zeval, kd0_Z, const_V1,
                                           kernel, mathcal_V)

    // === Trapezoidal integration ===
    // Note: integrates B^2 (not B), matching the IMSE formula
    int_bias = _didhetero_trapz(zeval, mathcal_B :^ 2)
    int_var  = _didhetero_trapz(zeval, mathcal_V)

    // === Edge cases ===
    h_max = (max(Z) - min(Z)) / 2

    if (int_bias >= . | int_bias < 1e-20) return(h_max)
    if (int_var >= . | int_var <= 0) return(.)

    // === IMSE1/US1 formula ===
    // h = (int_V / (4 * int_B^2))^(1/5) * n^(n_exponent)
    // IMSE1: n_exponent = -1/5, US1: n_exponent = -2/7
    h_opt = (int_var / (4 * int_bias))^(1/5) * n^(n_exponent)

    if (h_opt >= . | h_opt <= 0) return(.)

    return(h_opt)
}


// --- 5.3b IMSE1 bandwidth selection (thin wrapper) ---
// Calls the shared internal function with n_exponent = -1/5
real scalar _didhetero_bwselect_imse1(real matrix B_g_t,
                                      real colvector Z,
                                      real colvector zeval,
                                      real colvector kd0_Z,
                                      string scalar kernel)
{
    return(_dh_bwselect_imse1_internal(B_g_t, Z, zeval, kd0_Z, kernel, -1/5))
}


// --- 5.3c US1 bandwidth selection (undersmoothing) ---
// US1 uses the same bias/variance estimation as IMSE1 but with
// n^(-2/7) instead of n^(-1/5), producing a smaller bandwidth
// that satisfies the uniform inference condition.
//
// Paper ref: Section 4.2.6 (Undersmoothing bandwidth)
real scalar _didhetero_bwselect_us1(real matrix B_g_t,
                                     real colvector Z,
                                     real colvector zeval,
                                     real colvector kd0_Z,
                                     string scalar kernel)
{
    return(_dh_bwselect_imse1_internal(B_g_t, Z, zeval, kd0_Z, kernel, -2/7))
}


// --- 5.4 IMSE2 bandwidth selection (LQR bias) ---
// Computes the IMSE-optimal bandwidth for local quadratic regression.
//
// Algorithm (Paper Section 4.2.5, LQR case):
//   For each evaluation point r = 1..R:
//     Bias:  3rd deriv (p=4, deriv=3) + 4th deriv (p=5, deriv=4)
//            → B[r] = (1/(24*f(z_r))) * (2*mu3*f'(z_r) + mu4*f(z_r)) * C_B_LQ
//     Var:   Shared variance estimation with const_V = C_V2
//   Then: h = (int_V / (8 * int_B^2))^(1/9) * n^(-1/9)
//
// Args:
//   B_g_t  - n x R matrix, each column is the influence function
//   Z      - n x 1 covariate vector
//   zeval  - R x 1 sorted evaluation points
//   kd0_Z  - R x 1 kernel density estimates at zeval
//   kd1_Z  - R x 1 kernel density DERIVATIVE estimates at zeval
//   kernel - "epa" or "gau"
//
// Returns: scalar bandwidth (or . on failure)
real scalar _didhetero_bwselect_imse2(real matrix B_g_t,
                                      real colvector Z,
                                      real colvector zeval,
                                      real colvector kd0_Z,
                                      real colvector kd1_Z,
                                      string scalar kernel)
{
    real scalar n, R_eval, r
    real scalar I_2_K1, I_4_K1, I_6_K1, I_0_K2, I_2_K2, I_4_K2, I_6_K2
    real scalar C_B_LQ, const_V2
    real scalar cb_num, cb_den, cv_num, cv_den
    real colvector mathcal_B, mathcal_V
    real colvector y_r
    real colvector mu_B_3_bw_vec, mu_B_4_bw_vec
    real scalar mu_B_3_bw, mu_B_4_bw
    real scalar mu_B_3, mu_B_4
    real scalar int_bias, int_var, h_opt, h_max

    n = rows(Z)
    R_eval = rows(zeval)

    // Edge case: need at least 2 eval points for trapz
    if (R_eval < 2) return(.)

    // Get kernel constants
    _didhetero_kernel_constants(kernel, I_2_K1, I_4_K1, I_6_K1,
                                I_0_K2, I_2_K2, I_4_K2, I_6_K2)

    // Compute C_{B,LQ} = (I_4_K1^2 - I_2_K1*I_6_K1) / (I_4_K1 - I_2_K1^2)
    // For epa: -1/21 ≈ -0.04762
    cb_num = I_4_K1^2 - I_2_K1 * I_6_K1
    cb_den = I_4_K1 - I_2_K1^2
    if (abs(cb_den) < 1e-30) return(.)
    C_B_LQ = cb_num / cb_den

    // Compute C_{V2} = (I_4_K1^2*I_0_K2 - 2*I_2_K1*I_4_K1*I_2_K2 + I_2_K1^2*I_4_K2)
    //                  / (I_4_K1 - I_2_K1^2)^2
    // For epa: 1.25
    cv_num = I_4_K1^2 * I_0_K2 - 2 * I_2_K1 * I_4_K1 * I_2_K2 + I_2_K1^2 * I_4_K2
    cv_den = (I_4_K1 - I_2_K1^2)^2
    if (abs(cv_den) < 1e-30) return(.)
    const_V2 = cv_num / cv_den

    // Initialize bias vector
    mathcal_B = J(R_eval, 1, .)

    // === Bias estimation loop (LQR bias) ===
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // Skip if density is non-positive
        if (kd0_Z[r] >= . | kd0_Z[r] <= 0) continue

        // --- 3rd derivative: MSE-DPI bw (p=4, deriv=3) → LPR ---
        mu_B_3_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 4, 3, kernel)
        mu_B_3_bw = mu_B_3_bw_vec[1]

        if (mu_B_3_bw >= . | mu_B_3_bw <= 0) continue

        {
            real colvector mu_B_3_vec
            mu_B_3_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             4, 3, kernel, mu_B_3_bw)
            mu_B_3 = mu_B_3_vec[1]
        }

        if (mu_B_3 >= .) continue

        // --- 4th derivative: MSE-DPI bw (p=5, deriv=4) → LPR ---
        mu_B_4_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 5, 4, kernel)
        mu_B_4_bw = mu_B_4_bw_vec[1]

        if (mu_B_4_bw >= . | mu_B_4_bw <= 0) continue

        {
            real colvector mu_B_4_vec
            mu_B_4_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             5, 4, kernel, mu_B_4_bw)
            mu_B_4 = mu_B_4_vec[1]
        }

        if (mu_B_4 >= .) continue

        // mathcal_B[r] = (1/(24*f(z_r))) * (2*mu3*f'(z_r) + mu4*f(z_r)) * C_B_LQ
        mathcal_B[r] = (1 / (24 * kd0_Z[r])) * (2 * mu_B_3 * kd1_Z[r] + mu_B_4 * kd0_Z[r]) * C_B_LQ
    }

    // === Variance estimation via shared function (with C_V2) ===
    _didhetero_bwselect_var_est(B_g_t, Z, zeval, kd0_Z, const_V2,
                                           kernel, mathcal_V)

    // === Trapezoidal integration ===
    int_bias = _didhetero_trapz(zeval, mathcal_B :^ 2)
    int_var  = _didhetero_trapz(zeval, mathcal_V)

    // === Edge cases ===
    h_max = (max(Z) - min(Z)) / 2

    if (int_bias >= . | int_bias < 1e-20) return(h_max)
    if (int_var >= . | int_var <= 0) return(.)

    // === IMSE2 formula ===
    // h = (int_V / (8 * int_B^2))^(1/9) * n^(-1/9)
    // Note: denominator is 8 (not 4), exponent is 1/9 (not 1/5)
    h_opt = (int_var / (8 * int_bias))^(1/9) * n^(-1/9)

    if (h_opt >= . | h_opt <= 0) return(.)

    return(h_opt)
}


// --- 5.5 Bandwidth selection dispatch function ---
// Routes bwselect string to the appropriate bandwidth selection algorithm.
//
// Args:
//   bwselect  - "IMSE1", "IMSE2", "US1", or "manual"
//   B_g_t     - n x R matrix, each column is the influence function
//   Z         - n x 1 covariate vector
//   zeval     - R x 1 sorted evaluation points
//   kd0_Z     - R x 1 kernel density estimates at zeval
//   kd1_Z     - R x 1 kernel density DERIVATIVE estimates at zeval
//               (only needed for IMSE2; ignored otherwise)
//   kernel    - "epa" or "gau"
//   bw_manual - (optional) manual bandwidth value; required when bwselect="manual"
//
// Returns: scalar bandwidth (or . on failure)
real scalar _didhetero_bwselect(string scalar bwselect,
                                real matrix B_g_t,
                                real colvector Z,
                                real colvector zeval,
                                real colvector kd0_Z,
                                real colvector kd1_Z,
                                string scalar kernel,
                                | real scalar bw_manual)
{
    if (bwselect == "IMSE1") {
        return(_didhetero_bwselect_imse1(B_g_t, Z, zeval, kd0_Z, kernel))
    }
    else if (bwselect == "IMSE2") {
        return(_didhetero_bwselect_imse2(B_g_t, Z, zeval, kd0_Z, kd1_Z, kernel))
    }
    else if (bwselect == "US1") {
        return(_didhetero_bwselect_us1(B_g_t, Z, zeval, kd0_Z, kernel))
    }
    else if (bwselect == "manual") {
        if (args() < 8) {
            _error(3001, "bwselect='manual' requires bw_manual argument")
        }
        if (bw_manual >= .) {
            _error(3001, "bwselect='manual' requires non-missing bw_manual")
        }
        return(bw_manual)
    }
    else {
        _error(3498, "unknown bwselect: " + bwselect)
    }
    return(.)
}


// --- 5.6 Bandwidth selection for all (g,t) pairs ---
// Outer loop that computes bandwidth for each (g,t) pair, handles manual
// bandwidth broadcasting, and applies uniformall (minimum) logic.
//
// Args:
//   bwselect     - "IMSE1" | "IMSE2" | "US1" | "manual"
//   B_g_t_ptrs   - 1 x K pointer rowvector, each pointing to n x R matrix
//   Z            - n x 1 covariate vector
//   zeval        - R x 1 sorted evaluation points
//   kd0_Z        - R x 1 kernel density estimates at zeval
//   kd1_Z        - R x 1 kernel density derivative estimates at zeval
//   kernel       - "epa" or "gau"
//   uniformall   - 1 to take min across all (g,t), 0 for independent
//   bw_manual    - (optional) scalar or K x 1 vector of manual bandwidths
//
// Returns: K x 1 vector of bandwidths
//
// Manual bandwidth handling:
//   - Scalar bw_manual: broadcast to all K (g,t) pairs
//   - Vector bw_manual (K x 1): each (g,t) uses its own value
//   - Length mismatch: error
//
// Uniformall logic (Paper Section 4.2.5):
//   - uniformall=1: all (g,t) pairs use min(bw_vec)
//   - Also applies when bwselect="manual" and bw is scalar
real colvector _didhetero_bwselect_all(
    string scalar bwselect,
    pointer(real matrix) rowvector B_g_t_ptrs,
    real colvector Z,
    real colvector zeval,
    real colvector kd0_Z,
    real colvector kd1_Z,
    string scalar kernel,
    real scalar uniformall,
    | real colvector bw_manual)
{
    real scalar K, id_gt
    real colvector bw_vec

    K = cols(B_g_t_ptrs)
    bw_vec = J(K, 1, .)

    // === Manual bandwidth handling ===
    if (bwselect == "manual") {
        if (args() < 9) {
            _error(3001, "bwselect='manual' requires bw_manual argument")
        }
        if (rows(bw_manual) == 1) {
            // Scalar: broadcast to all (g,t) pairs
            bw_vec = J(K, 1, bw_manual[1])
        }
        else if (rows(bw_manual) == K) {
            // Vector: one per (g,t) pair
            bw_vec = bw_manual
        }
        else {
            _error(3200, "bw must be scalar or vector of length " + strofreal(K))
        }
    }
    else {
        // === Algorithmic bandwidth selection ===
        for (id_gt = 1; id_gt <= K; id_gt++) {
            bw_vec[id_gt] = _didhetero_bwselect(bwselect,
                *B_g_t_ptrs[id_gt], Z, zeval, kd0_Z, kd1_Z, kernel)
        }
    }

    // === Uniformall: take minimum across all (g,t) pairs ===
    // Paper Section 4.2.5: common bandwidth for uniform inference
    if (uniformall) {
        bw_vec = J(K, 1, min(bw_vec))
    }

    return(bw_vec)
}


// --- 4.5 didhetero_lpdensity: density derivative estimation ---
// Density derivative estimation with MSE-DPI bandwidth selection.
// Uses CDF-based LPR with unscaled basis (our didhetero_lpr)
real colvector didhetero_lpdensity(real colvector X, real colvector eval,
                                    real scalar p, real scalar v,
                                    string scalar kernel)
{
    real colvector F_n, h_bw, lpr_est
    real scalar kernel_mapped_ok
    string scalar kernel_mapped

    // Map kernel name: "epanechnikov" is aliased to "epa"
    kernel_mapped = kernel
    if (kernel == "epanechnikov") kernel_mapped = "epa"

    // Step 1: Compute ECDF
    F_n = _didhetero_ecdf(X)

    // Step 2: MSE-DPI bandwidth selection
    h_bw = _didhetero_lpdensity_bw_MSE(X, eval, p, v, kernel_mapped)

    // Step 3: CDF-based LPR estimation
    // Our didhetero_lpr with unscaled basis returns v! * beta_v directly
    // which equals the v-th derivative of the density
    // (equivalent to the scaled basis approach in the density estimation literature)
    lpr_est = didhetero_lpr(F_n, X, eval, p, v, kernel_mapped, h_bw)

    return(lpr_est)
}


end
