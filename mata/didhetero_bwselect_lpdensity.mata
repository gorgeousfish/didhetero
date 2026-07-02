mata:

// =============================================================================
// didhetero_bwselect_lpdensity.mata
// IMSE/US bandwidth selection and local polynomial density estimation
//
// Provides:
//   - _didhetero_kernel_constants()   // Kernel integral constants for BW
//   - _didhetero_bwselect_var_est()   // Variance estimation for IMSE/US
//   - _didhetero_bwselect_imse1()     // IMSE1 bandwidth (LLR optimal)
//   - _didhetero_bwselect_us1()       // US1 undersmoothing bandwidth
//   - _didhetero_bwselect_imse2()     // IMSE2 bandwidth (LQR optimal)
//   - _didhetero_bwselect()           // Single (g,t) bandwidth dispatcher
//   - _didhetero_bwselect_all()       // All (g,t) pairs bandwidth selection
//   - didhetero_lpdensity()           // Local polynomial density estimation
//
// Requires:
//   - didhetero_kernel.mata       (didhetero_kernel_eval)
//   - didhetero_lpr.mata          (didhetero_lpr)
//   - didhetero_bwselect_lp.mata  (_didhetero_lpbwselect_mse, _didhetero_lpbwselect_imse,
//                                  _didhetero_nn_variance, _didhetero_bwcheck)
//   - didhetero_errors.mata       (_dh_error_input)
//
// Paper reference: Section 3.2, bandwidth selection methods
// =============================================================================


// =====================================================================
// Section 5: IMSE1 / IMSE2 / US1 bandwidth selection
// =====================================================================

// --- 5.1 Kernel integral constants ---
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
        _dh_error_input("bwselect", "unknown kernel: " + kernel)
    }
}


// --- 5.2 Variance estimation for IMSE1/IMSE2/US1 ---
// Estimates conditional variance at each evaluation point via five-step
// procedure: IMSE-DPI bandwidth for conditional mean, residuals, MSE-DPI
// bandwidth for conditional variance, and LPR estimation scaled by density.
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

    // Support grid (100 equispaced points)
    Z_supp = rangen(min(Z), max(Z), 100)

    // Loop over evaluation points
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // Skip if density is non-positive
        if (kd0_Z[r] >= . | kd0_Z[r] <= 0) continue

        // Bandwidth for conditional mean
        mu_B_0_bw = _didhetero_lpbwselect_imse(y_r, Z, Z_supp,
                        1, 0, kernel)

        if (mu_B_0_bw >= . | mu_B_0_bw <= 0) continue

        // Local linear regression at observation points
        mu_B_0 = didhetero_lpr(y_r, Z, Z, 1, 0, kernel, mu_B_0_bw)

        if (rows(mu_B_0) != n) continue

        // Residuals
        U_hat = y_r - mu_B_0

        // Bandwidth for conditional variance
        U_hat_sq = U_hat :^ 2
        sigma2_bw_vec = _didhetero_lpbwselect_mse(U_hat_sq, Z,
                            zeval[r..r], 1, 0, kernel)
        sigma2_bw = sigma2_bw_vec[1]

        if (sigma2_bw >= . | sigma2_bw <= 0) continue

        // Estimate conditional variance
        {
            real colvector sigma2_vec
            sigma2_vec = didhetero_lpr(U_hat_sq, Z, zeval[r..r],
                             1, 0, kernel, sigma2_bw)
            sigma2 = sigma2_vec[1]
        }

        if (sigma2 >= . | sigma2 <= 0) continue

            mathcal_V[r] = const_V * sigma2 / kd0_Z[r]
    }
}


// --- 5.3a IMSE1/US1 internal implementation ---
// Bandwidth selection minimizing approximate IMSE. IMSE1 uses n^(-1/5),
// US1 uses n^(-2/7) for undersmoothing.
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

    // Bias estimation
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // Bias estimation
        // Bandwidth for second derivative
        mu_B_2_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 3, 2, kernel)
        mu_B_2_bw = mu_B_2_bw_vec[1]

        if (mu_B_2_bw >= . | mu_B_2_bw <= 0) continue

        // Estimate second derivative
        {
            real colvector mu_B_2_vec
            mu_B_2_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             3, 2, kernel, mu_B_2_bw)
            mu_B_2 = mu_B_2_vec[1]
        }

        if (mu_B_2 >= .) continue

        mathcal_B[r] = mu_B_2 * I_2_K1 / 2
    }

    // Variance estimation
    _didhetero_bwselect_var_est(B_g_t, Z, zeval, kd0_Z, const_V1,
                                           kernel, mathcal_V)

    // Trapezoidal integration of squared bias and variance
    int_bias = _didhetero_trapz(zeval, mathcal_B :^ 2)
    int_var  = _didhetero_trapz(zeval, mathcal_V)

    // Edge case handling
    h_max = (max(Z) - min(Z)) / 2

    if (int_bias >= . | int_bias < 1e-20) return(h_max)
    if (int_var >= . | int_var <= 0) return(.)

    // Final bandwidth formula
    h_opt = (int_var / (4 * int_bias))^(1/5) * n^(n_exponent)

    if (h_opt >= . | h_opt <= 0) return(.)

    return(h_opt)
}


// --- 5.3b IMSE1 bandwidth selection ---
real scalar _didhetero_bwselect_imse1(real matrix B_g_t,
                                      real colvector Z,
                                      real colvector zeval,
                                      real colvector kd0_Z,
                                      string scalar kernel)
{
    return(_dh_bwselect_imse1_internal(B_g_t, Z, zeval, kd0_Z, kernel, -1/5))
}


// --- 5.3c US1 bandwidth selection (undersmoothing) ---
real scalar _didhetero_bwselect_us1(real matrix B_g_t,
                                     real colvector Z,
                                     real colvector zeval,
                                     real colvector kd0_Z,
                                     string scalar kernel)
{
    return(_dh_bwselect_imse1_internal(B_g_t, Z, zeval, kd0_Z, kernel, -2/7))
}


// --- 5.4 IMSE2 bandwidth selection (local quadratic regression) ---
// Computes IMSE-optimal bandwidth using third and fourth derivatives
// of the influence function.
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
    real scalar _verbose, _n_reg_skip

    n = rows(Z)
    R_eval = rows(zeval)

    // Verbose diagnostics (read from Stata local set by ado layer)
    _verbose = strtoreal(st_local("_dh_verbose"))
    if (_verbose >= .) _verbose = 0
    _n_reg_skip = 0

    // Edge case: need at least 2 eval points for trapz
    if (R_eval < 2) return(.)

    // Get kernel constants
    _didhetero_kernel_constants(kernel, I_2_K1, I_4_K1, I_6_K1,
                                I_0_K2, I_2_K2, I_4_K2, I_6_K2)

    // Bias constant C_{B,LQ}
    cb_num = I_4_K1^2 - I_2_K1 * I_6_K1
    cb_den = I_4_K1 - I_2_K1^2
    if (abs(cb_den) < 1e-30) return(.)
    C_B_LQ = cb_num / cb_den

    // Variance constant C_{V2}
    cv_num = I_4_K1^2 * I_0_K2 - 2 * I_2_K1 * I_4_K1 * I_2_K2 + I_2_K1^2 * I_4_K2
    cv_den = (I_4_K1 - I_2_K1^2)^2
    if (abs(cv_den) < 1e-30) return(.)
    const_V2 = cv_num / cv_den

    // Initialize bias vector
    mathcal_B = J(R_eval, 1, .)

    // Bias estimation (local quadratic regression)
    for (r = 1; r <= R_eval; r++) {

        y_r = B_g_t[., r]

        // Skip if density is non-positive
        if (kd0_Z[r] >= . | kd0_Z[r] <= 0) continue

        // Third derivative bandwidth and estimation
        mu_B_3_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 4, 3, kernel)
        mu_B_3_bw = mu_B_3_bw_vec[1]

        if (mu_B_3_bw >= . | mu_B_3_bw <= 0) {
            if (_verbose) {
                _n_reg_skip++
                printf("{txt}  [bwselect] IMSE2: 3rd deriv bw failed at z=%.4f (r=%g/%g)\n",
                       zeval[r], r, R_eval)
            }
            continue
        }

        {
            real colvector mu_B_3_vec
            mu_B_3_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             4, 3, kernel, mu_B_3_bw)
            mu_B_3 = mu_B_3_vec[1]
        }

        if (mu_B_3 >= .) {
            if (_verbose) {
                _n_reg_skip++
                printf("{txt}  [bwselect] IMSE2: 3rd deriv estimate failed at z=%.4f (bw=%.4f)\n",
                       zeval[r], mu_B_3_bw)
            }
            continue
        }

        // Fourth derivative bandwidth and estimation
        mu_B_4_bw_vec = _didhetero_lpbwselect_mse(y_r, Z,
                            zeval[r..r], 5, 4, kernel)
        mu_B_4_bw = mu_B_4_bw_vec[1]

        if (mu_B_4_bw >= . | mu_B_4_bw <= 0) {
            if (_verbose) {
                _n_reg_skip++
                printf("{txt}  [bwselect] IMSE2: 4th deriv bw failed at z=%.4f (r=%g/%g)\n",
                       zeval[r], r, R_eval)
            }
            continue
        }

        {
            real colvector mu_B_4_vec
            mu_B_4_vec = didhetero_lpr(y_r, Z, zeval[r..r],
                             5, 4, kernel, mu_B_4_bw)
            mu_B_4 = mu_B_4_vec[1]
        }

        if (mu_B_4 >= .) {
            if (_verbose) {
                _n_reg_skip++
                printf("{txt}  [bwselect] IMSE2: 4th deriv estimate failed at z=%.4f (bw=%.4f)\n",
                       zeval[r], mu_B_4_bw)
            }
            continue
        }

        mathcal_B[r] = (1 / (24 * kd0_Z[r])) * (2 * mu_B_3 * kd1_Z[r] + mu_B_4 * kd0_Z[r]) * C_B_LQ
    }

    // Diagnostic summary
    if (_verbose & _n_reg_skip > 0) {
        printf("{txt}  [bwselect] IMSE2: %g/%g eval points skipped (boundary/ill-conditioned)\n",
               _n_reg_skip, R_eval)
    }

    // Variance estimation
    _didhetero_bwselect_var_est(B_g_t, Z, zeval, kd0_Z, const_V2,
                                           kernel, mathcal_V)

    // Trapezoidal integration
    int_bias = _didhetero_trapz(zeval, mathcal_B :^ 2)
    int_var  = _didhetero_trapz(zeval, mathcal_V)

    // Edge case handling
    h_max = (max(Z) - min(Z)) / 2

    if (int_bias >= . | int_bias < 1e-20) return(h_max)
    if (int_var >= . | int_var <= 0) return(.)

    // Final bandwidth formula
    h_opt = (int_var / (8 * int_bias))^(1/9) * n^(-1/9)

    if (h_opt >= . | h_opt <= 0) return(.)

    return(h_opt)
}


// --- 5.5 Bandwidth selection dispatch function ---
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
            _dh_error_dimension("bwselect", "bwselect='manual' requires bw_manual argument")
        }
        if (bw_manual >= .) {
            _dh_error_dimension("bwselect", "bwselect='manual' requires non-missing bw_manual")
        }
        return(bw_manual)
    }
    else {
        _dh_error_estimation("bwselect", "unknown bwselect: " + bwselect)
    }
    return(.)
}


// --- 5.6 Bandwidth selection for all (g,t) pairs ---
// Computes bandwidth for each (g,t) pair with manual bandwidth broadcasting
// and uniformall (common minimum) logic.
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

    // Manual bandwidth handling
    if (bwselect == "manual") {
        if (args() < 9) {
            _dh_error_dimension("bwselect", "bwselect='manual' requires bw_manual argument")
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
            _dh_error_data("bwselect", "bw must be scalar or vector of length " + strofreal(K))
        }
    }
    else {
        // Algorithmic bandwidth selection
        for (id_gt = 1; id_gt <= K; id_gt++) {
            bw_vec[id_gt] = _didhetero_bwselect(bwselect,
                *B_g_t_ptrs[id_gt], Z, zeval, kd0_Z, kd1_Z, kernel)
        }
    }

    // Common bandwidth (minimum across all groups and periods)
    if (uniformall) {
        bw_vec = J(K, 1, min(bw_vec))
    }

    // Undersmoothing condition diagnostic (Assumption 4(iii))
    _didhetero_bw_check_undersmooth(bw_vec, rows(Z), bwselect)

    return(bw_vec)
}


// --- 4.5 Density derivative estimation ---
// Computes density derivative estimates using CDF-based local polynomial
// regression with MSE-DPI bandwidth selection.
real colvector didhetero_lpdensity(real colvector X, real colvector eval,
                                    real scalar p, real scalar v,
                                    string scalar kernel)
{
    real colvector F_n, h_bw, lpr_est
    real scalar kernel_mapped_ok
    string scalar kernel_mapped

    // Kernel name mapping
    kernel_mapped = kernel
    if (kernel == "epanechnikov") kernel_mapped = "epa"

    // Empirical CDF
    F_n = _didhetero_ecdf(X)

    // MSE-DPI bandwidth selection
    h_bw = _didhetero_lpdensity_bw_MSE(X, eval, p, v, kernel_mapped)

    // CDF-based local polynomial regression
    lpr_est = didhetero_lpr(F_n, X, eval, p, v, kernel_mapped, h_bw)

    return(lpr_est)
}


end
mata:

// =============================================================================
// didhetero_bwselect_lpdensity.mata
// Local polynomial density estimation and IMSE/US bandwidth selection
//
// Dependencies:
//   - didhetero_kernel.mata (must be compiled first)
//   - didhetero_lpr.mata (must be compiled first)
//   - didhetero_bwselect_lp.mata (must be compiled first)
//
// Contains:
//   Section 4: Local polynomial density derivative estimation
//     - _didhetero_lpdensity_Sgenerate, Tgenerate, Cgenerate, Ggenerate
//     - _didhetero_normal_pdf_deriv
//     - _didhetero_lpdensity_bw_IROT
//     - _didhetero_lpdensity_bw_MSE
//     - didhetero_lpdensity (public interface)
//   Section 5: IMSE1 / IMSE2 / US1 bandwidth selection
//     - _didhetero_kernel_constants
//     - _didhetero_bwselect_var_est
//     - _dh_bwselect_imse1_internal
//     - _didhetero_bwselect_imse1, _didhetero_bwselect_us1
//     - _didhetero_bwselect_imse2
//     - _didhetero_bwselect (router)
//     - _didhetero_bwselect_all
// =============================================================================

// =====================================================================
// Section 4: Local polynomial density derivative estimation
// Implements bw_IROT, bw_MSE, and density estimation functions.
// =====================================================================

// --- 4.1 Kernel matrix generators ---

// S matrix: kernel moment integrals
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

// T matrix: squared kernel moment integrals
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

// C vector: bias term kernel integrals
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

// G matrix: influence function variance (double integral via midpoint rule)
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

// --- 4.3 Integrated Rule-of-Thumb bandwidth ---
// Returns scalar bandwidth integrated over grid points.
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

    // Standardize data
    center_temp = mean(data_local)
    scale_temp = sqrt(variance(data_local))
    data_local = (data_local :- center_temp) / scale_temp
    dataUnique = (dataUnique :- center_temp) / scale_temp
    grid_local = (grid_local :- center_temp) / scale_temp

    // Normal reference model
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

    // Handle missing values
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


// --- 4.4 MSE-optimal pointwise bandwidth ---
// Returns ng x 1 vector of bandwidths.
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

    // Sort data
    data = sort(data_in, 1)
    grid = grid_in  // don't sort grid, keep original order
    n = rows(data)
    ng = rows(grid)

    // Unique values with frequencies
    dataUnique = uniqrows(data)
    nUnique = rows(dataUnique)

    // Frequency and index vectors for tied values
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

    // Standardize data
    center_temp = mean(data)
    scale_temp = sqrt(variance(data))
    data = (data :- center_temp) / scale_temp
    dataUnique = (dataUnique :- center_temp) / scale_temp
    grid = (grid :- center_temp) / scale_temp

    // Empirical CDF (with ties receiving the same cumulative proportion)
    Fn = (1::n) / n
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

    // IROT pilot bandwidths
    nLocalMin_h1 = 20 + 2 + 1
    nUniqueMin_h1 = 20 + 2 + 1
    h1 = _didhetero_lpdensity_bw_IROT(data, grid, 2, 1, kernel, nLocalMin_h1, nUniqueMin_h1)

    nLocalMin_hp1 = 20 + p + 2 + 1
    nUniqueMin_hp1 = 20 + p + 2 + 1
    hp1 = _didhetero_lpdensity_bw_IROT(data, grid, p + 2, p + 1, kernel, nLocalMin_hp1, nUniqueMin_hp1)

    nLocalMin_hp2 = 20 + p + 3 + 1
    nUniqueMin_hp2 = 20 + p + 3 + 1
    hp2 = _didhetero_lpdensity_bw_IROT(data, grid, p + 3, p + 2, kernel, nLocalMin_hp2, nUniqueMin_hp2)

    // Per grid point: estimate data generating process and constants
    dgp_hat = J(ng, 2, .)
    const_hat = J(ng, 3, .)
    h = J(ng, 1, .)
    range_data = max(data) - min(data)
    gr = (sqrt(5) + 1) / 2

    for (j = 1; j <= ng; j++) {

        // Estimate CDF derivative F_{p+2}
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

        // Polynomial basis
        Xh_p_temp = J(n_eff, p + 4, .)
        for (k = 0; k <= p + 3; k++) {
            Xh_p_temp[., k + 1] = Xh_temp :^ k
        }

        // Kernel weights
        Kh_temp = 0.75 * (1 :- Xh_temp :^ 2) / hp2

        // Weighted least squares
        beta_temp = cholinv(cross(Xh_p_temp, Kh_temp, Xh_p_temp)) * cross(Xh_p_temp, Kh_temp, Y_temp)
        if (hasmissing(beta_temp)) continue
        dgp_hat[j, 2] = beta_temp[p + 3] / hp2^(p + 2)

        // Estimate CDF derivative F_{p+1}
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

        // Pre-asymptotic matrices
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
            // Influence function approach with tied values
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

        // Assemble bias and variance constants
        const_hat[j, 1] = factorial(v) * (S_hat_inv * C_p_hat)[v + 1]
        const_hat[j, 2] = factorial(v) * (S_hat_inv * C_p1_hat)[v + 1]

        if (v > 0) {
            const_hat[j, 3] = factorial(v) * sqrt(abs((S_hat_inv * G_hat * S_hat_inv)[v + 1, v + 1]) / (n * h1))
        }
        else {
            temp_val = min((max((mean(data :<= grid[j]), 1 / n)), 1 - 1 / n))
            const_hat[j, 3] = factorial(v) * sqrt(abs((S_hat_inv * G_hat * S_hat_inv)[v + 1, v + 1] / (0.5 * n^2) * h1 * temp_val * (1 - temp_val)))
        }

        // Golden section search for optimal bandwidth
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

    // Post-processing regularization
    for (j = 1; j <= ng; j++) {
        sorted_abs = sort(abs(data :- grid[j]), 1)
        if (h[j] == .) {
            h[j] = sorted_abs[min((n, max((20 + p + 1, 20 + p + 1))))]
        }
        // Minimum local observations regularization
        bw_min_local = sorted_abs[min((n, 20 + p + 1))]
        if (h[j] < bw_min_local) h[j] = bw_min_local
        // Minimum unique observations regularization
        sorted_abs_u = sort(abs(dataUnique :- grid[j]), 1)
        bw_min_unique = sorted_abs_u[min((nUnique, 20 + p + 1))]
        if (h[j] < bw_min_unique) h[j] = bw_min_unique
        // Cap at max distance
        h_max_j = max(abs(dataUnique :- grid[j]))
        if (h[j] > h_max_j) h[j] = h_max_j
    }

    // Rescale to original units
    h = h * scale_temp

    return(h)
}


end
