mata:

// =============================================================================
// didhetero_bwselect_kde.mata
// Kernel density bandwidth selection (MSE-DPI)
//
// Provides:
//   - _didhetero_kdbwselect_mse()  // KDE MSE-DPI bandwidth
//   - didhetero_kdrobust()         // Public KDE interface with automatic bandwidth
//
// Requires:
//   - didhetero_kernel.mata        (didhetero_kernel_eval)
//   - didhetero_bwselect_lp.mata   (_didhetero_bwcheck)
// =============================================================================

// =====================================================================
// Section 3: kdbwselect MSE-DPI bandwidth for kernel density estimation
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

// Kernel density estimation with MSE-DPI bandwidth selection
real colvector didhetero_kdrobust(real colvector X, real colvector eval,
                                  string scalar kernel)
{
    real colvector h_bw
    h_bw = _didhetero_kdbwselect_mse(X, eval, kernel)
    return(didhetero_kde(X, eval, kernel, h_bw))
}


end
