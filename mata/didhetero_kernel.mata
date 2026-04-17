mata:

// =============================================================================
// didhetero_kernel.mata
// Kernel functions and precomputed integral constants
//
// This file provides two core functions:
//   1. didhetero_kernel_eval()   - Evaluate kernel function at given points
//   2. didhetero_kernel_consts() - Initialize precomputed kernel integral constants
//
// Supported kernels:
//   "epa" - Epanechnikov: K(u) = 0.75 * (1 - u^2) * I(|u| <= 1)
//   "gau" - Gaussian:     K(u) = phi(u) (standard normal density)
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//     - Section 4.2.2 Eq. 18-19 (LQR bias/variance constants)
//     - Section 4.2.6 / Appendix A (LLR bias/variance constants)
//     - Eq. 21-22 (lambda, analytical UCB)
//   Section 4.2.2 Eq. 18-19, Section 4.2.6, Appendix A
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_kernel_eval()
// Evaluate kernel function at given points.
//
// Args:
//   u      - real matrix of standardized distances (any dimension)
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//
// Returns:
//   real matrix of kernel values, same dimension as u
//
// Paper ref: Gaussian and Epanechnikov kernel functions
// -----------------------------------------------------------------------------
real matrix didhetero_kernel_eval(real matrix u, string scalar kernel)
{
    if (kernel == "epa") {
        // K(u) = 0.75 * (1 - u^2) * I(|u| <= 1)
        // Vectorized via element-wise operators
        return((abs(u) :<= 1) :* 0.75 :* (1 :- u:^2))
    }
    else if (kernel == "gau") {
        // K(u) = phi(u) = standard normal density
        return(normalden(u))
    }
    else {
        _error("invalid kernel: " + kernel)
    }
}

// -----------------------------------------------------------------------------
// didhetero_kernel_consts()
// Initialize kernel integral constants struct.
//
// Computes base integral moments and derives LLR/LQR bias and variance
// constants at runtime from those base integrals.
//
// Base integrals:
//   I_{l,K^m} = int u^l [K(u)]^m du
//
// Derived constants:
//   const_B1 = I_{2,K} / 2                          (LLR bias)
//   const_V1 = I_{0,K^2}                             (LLR variance)
//   const_B2 = (I_{4,K}^2 - I_{2,K}*I_{6,K})        (LQR bias)
//              / (I_{4,K} - I_{2,K}^2)
//   const_V2 = (I_{4,K}^2*I_{0,K^2}                  (LQR variance)
//               - 2*I_{2,K}*I_{4,K}*I_{2,K^2}
//               + I_{2,K}^2*I_{4,K^2})
//              / (I_{4,K} - I_{2,K}^2)^2
//
// Args:
//   kernel - "epa" (Epanechnikov) or "gau" (Gaussian)
//
// Returns:
//   DidHeteroKernelConsts struct with all fields populated
//
// Paper refs: Section 4.2.2 Eq. 18-19, Section 4.2.6, Appendix A
// cf. Imai, Qin, and Yanagi (2025), Sections 4.2.2 and 4.2.6
// -----------------------------------------------------------------------------
struct DidHeteroKernelConsts scalar didhetero_kernel_consts(string scalar kernel)
{
    struct DidHeteroKernelConsts scalar kc
    
    kc.name = kernel
    
    if (kernel == "epa") {
        // Epanechnikov: K(u) = 0.75*(1-u^2), |u|<=1
        // I_{l,K^1} = 3/((l+1)(l+3)) for even l
        kc.I_2_K  = 1/5       // 3/(3*5)
        kc.I_4_K  = 3/35      // 3/(5*7)
        kc.I_6_K  = 1/21      // 3/(7*9)
        // I_{l,K^2} via analytic integration of (3/4)^2*(1-u^2)^2
        kc.I_0_K2 = 3/5
        kc.I_2_K2 = 3/35      // coincides with I_4_K (mathematical coincidence)
        kc.I_4_K2 = 1/35
        kc.I_6_K2 = 1/77
        // lambda = -int(K*K''du)/int(K^2 du) = (3/2)/(3/5) = 5/2
        kc.lambda  = 5/2
    }
    else if (kernel == "gau") {
        // Gaussian: K(u) = phi(u)
        // I_{l,K^1} = E[U^l] for U~N(0,1) = (l-1)!! for even l
        kc.I_2_K  = 1
        kc.I_4_K  = 3
        kc.I_6_K  = 15
        // I_{l,K^2} = l!/(2^l*(l/2)!) / (2*sqrt(pi))
        kc.I_0_K2 = 1 / (2 * sqrt(c("pi")))
        kc.I_2_K2 = 1 / (4 * sqrt(c("pi")))
        kc.I_4_K2 = 3 / (8 * sqrt(c("pi")))
        kc.I_6_K2 = 15 / (16 * sqrt(c("pi")))
        // lambda = -int(K * K'' du) / int(K^2 du) = 1/2
        kc.lambda  = 1/2
    }
    else {
        _error("invalid kernel: " + kernel)
    }
    
    // --- Derived constants (computed from base integrals, not hardcoded) ---
    
    // LLR constants (Paper Section 4.2.6 / Appendix A, p=1 special case)
    kc.const_B1 = kc.I_2_K / 2                    // Bias: h^2 * (I_{2,K}/2) * mu_B''(z)
    kc.const_V1 = kc.I_0_K2                        // Var: I_{0,K^2}/(nh) * sigma_B^2/f_Z
    
    // LQR constants (Paper Section 4.2.2, Eq. 18-19)
    real scalar denom
    denom = kc.I_4_K - kc.I_2_K^2
    
    // Eq. 18: const_B2 = (I_{4,K}^2 - I_{2,K}*I_{6,K}) / (I_{4,K} - I_{2,K}^2)
    kc.const_B2 = (kc.I_4_K^2 - kc.I_2_K * kc.I_6_K) / denom
    
    // Eq. 19: const_V2 = (I_{4,K}^2*I_{0,K^2} - 2*I_{2,K}*I_{4,K}*I_{2,K^2}
    //                     + I_{2,K}^2*I_{4,K^2}) / (I_{4,K} - I_{2,K}^2)^2
    kc.const_V2 = (kc.I_4_K^2 * kc.I_0_K2 ///
                   - 2 * kc.I_2_K * kc.I_4_K * kc.I_2_K2 ///
                   + kc.I_2_K^2 * kc.I_4_K2) / denom^2
    
    return(kc)
}

end
