mata:

// =============================================================================
// Mata struct definitions for heterogeneous treatment effect estimation
//
// Data structures:
//   1. DidHeteroData           - Panel data container for estimation
//   2. DidHeteroParamResults   - Parametric estimation results (GPS + OR)
//   3. DidHeteroStage1Results  - Full Stage 1 results (parametric + KDE)
//   4. DidHeteroKernelConsts   - Precomputed kernel integral constants
//   5. BootPrecomp             - Bootstrap pre-computed invariants
//   6. DidHeteroEstResult      - Single (g,t) pair estimation result
//   7. DidHeteroAggteResult    - Aggregated parameter estimation result
//   8. DidHeteroCattResult     - Reentrant CATT estimator return structure
//   9. AggtResult              - aggte_gt command output (all eval points)
// =============================================================================

// -----------------------------------------------------------------------------
// DidHeteroData
// Encapsulates panel data and evaluation grids for estimation.
// Populated during data preparation.
// -----------------------------------------------------------------------------
struct DidHeteroData {

    // === Panel data ===
    real matrix    Y_wide      // n x T outcome matrix in wide format
    real colvector G           // n x 1 group variable (first treatment period)
                               // G[i] = 0 indicates never-treated units
    real colvector Z           // n x 1 continuous covariate (heterogeneity focus)
    real colvector id          // n x 1 panel unit identifier
    
    // === Dimensions ===
    real scalar    n           // number of cross-sectional units
    real scalar    T_num       // number of time periods
                               // Named T_num to avoid conflict with Mata's transpose operator T
    real colvector t_vals      // T_num x 1 sorted unique time period values
    real scalar    period1     // first period value = t_vals[1]

    // === Treatment structure ===
    real scalar    gbar        // max treatment period for comparison group construction
                               // = . (missing) if never-treated group exists (encodes +Inf)
                               // = max(unique(G[G>0])) otherwise
    real colvector geval       // valid group evaluation points
                               // = intersect(supp_g, supp_g + anticipation) \ {0}
                               // When anticipation=0: geval = supp_g \ {0}
    real colvector teval       // valid time evaluation points
                               // = intersect(supp_t, supp_t - anticipation) \ {period1}
                               // When anticipation=0: teval = supp_t \ {period1}

    // === Evaluation grid ===
    real colvector zeval       // R x 1 user-specified or auto-generated Z evaluation points
                               // Points where CATT results are reported
    real scalar    num_zeval   // number of evaluation points R = length(zeval)
    real colvector Z_supp      // 100 x 1 equally-spaced grid on [min(Z), max(Z)]
                               // Used for lpbwselect(eval=Z_supp, bwselect="imse-dpi")
                               // to select a common bandwidth (scalar)
                               // Not the same as zeval

    // === Support sets ===
    real colvector supp_g       // sorted unique G values (including 0)
    real colvector supp_t       // sorted unique time period values

    // === (g,t) evaluation pairs ===
    real matrix    gteval       // K x 2 matrix of valid (g, t) pairs
    real scalar    num_gteval   // number of (g, t) pairs = rows(gteval)

    // === Covariate matrix ===
    real matrix    X            // n x k covariate matrix (with intercept)

    // === Configuration parameters (passed from Ado layer) ===
    string scalar  control_group // "nevertreated" or "notyettreated"
    real scalar    anticipation  // anticipation parameter (>= 0, integer)
    real scalar    porder        // polynomial order (1 or 2)
    string scalar  kernel        // "epa" or "gau"
    real scalar    alp           // significance level (0, 1)
    real scalar    biters        // bootstrap iterations (>= 1)
    real scalar    uniformall    // uniform across all (g,t) pairs flag

    // === Compatibility fields retained for option/data plumbing ===
    string scalar  base_period      // "universal" or "varying"
    string scalar  true_base_period // stored copy of base_period
    real scalar    print_details    // 1 = print per-(g,t) status, 0 = quiet
    real scalar    deriv            // derivative order for LPR (default 0)
    real scalar    panel            // 1 = panel data, 0 = repeated cross-section
    string scalar  est_method       // estimation method identifier
    string scalar  weightsname      // name of weights variable ("" if none)
    real scalar    bstrap           // 1 = perform bootstrap, 0 = skip
    real scalar    cband            // 1 = compute confidence bands
    real scalar    nboot            // number of bootstrap iterations (alias for biters)
    real scalar    bw_scalar        // scalar bandwidth (user-specified or selected)

    // === Kernel-derived constants (from DidHeteroKernelConsts) ===
    real scalar    const_V       // selected variance constant (by porder)
    real scalar    const_B1      // LLR bias constant
    real scalar    const_B2      // LQR bias constant
    real scalar    lambda        // kernel lambda for analytical UCB

    // === Core estimation arrays ===
    pointer(real matrix) rowvector A_g_t    // 1 x num_gteval, each n x num_zeval
    pointer(real matrix) rowvector B_g_t    // 1 x num_gteval, each n x num_zeval
    real matrix    G_g           // n x num_gteval, group indicators
    real matrix    mu_G_g        // num_zeval x num_gteval
    real matrix    mu_E_g_t      // num_zeval x num_gteval
    real matrix    mu_F_g_t      // num_zeval x num_gteval

    // === Stage 1 density estimates ===
    real colvector kd0_Z        // num_zeval x 1, kernel density estimates

    // === Standard errors and analytical uniform confidence bands ===
    real matrix    se           // num_zeval x num_gteval, standard errors
    real matrix    mathcal_V    // num_zeval x num_gteval, variance estimates
    real matrix    ci1_lower    // num_zeval x num_gteval, analytical UCB lower bound
    real matrix    ci1_upper    // num_zeval x num_gteval, analytical UCB upper bound
    real colvector c_hat        // num_gteval x 1, analytical critical values

    // === Bootstrap uniform confidence bands ===
    real matrix    ci2_lower    // num_zeval x num_gteval, bootstrap UCB lower bound
    real matrix    ci2_upper    // num_zeval x num_gteval, bootstrap UCB upper bound
    real colvector c_check_bs   // num_gteval x 1, bootstrap critical values
}

// -----------------------------------------------------------------------------
// DidHeteroParamResults
// Intermediate return structure from parametric estimation (GPS + OR).
// Contains GPS and OR estimation results for all (g,t) pairs.
// -----------------------------------------------------------------------------
struct DidHeteroParamResults {
    real matrix    gps_mat      // GPS estimates: nevertreated (id,g,est) 3 cols,
                                // notyettreated (id,g,t,est) 4 cols
    real matrix    gps_coef     // Logit coefficients: nevertreated (g,coef1,...),
                                // notyettreated (g,t,coef1,...)
    real matrix    or_mat       // OR estimates: always (id,g,t,est) 4 cols
    real matrix    or_coef      // OLS coefficients: always (g,t,coef1,...)
    string scalar  ctrl_type    // "nevertreated" or "notyettreated"
}

// -----------------------------------------------------------------------------
// DidHeteroStage1Results
// Full Stage 1 results: parametric estimation + kernel density estimation.
// Returned by didhetero_stage1_dispatch().
// -----------------------------------------------------------------------------
struct DidHeteroStage1Results {
    real matrix    gps_mat      // GPS estimates: (id, g, [t,] est)
    real matrix    gps_coef     // Logit coefficients: (g, [t,] coef1, ...)
    real matrix    or_mat       // OR estimates: (id, g, t, est)
    real matrix    or_coef      // OLS coefficients: (g, t, coef1, ...)
    string scalar  ctrl_type    // "nevertreated" or "notyettreated"
    real colvector Z_supp       // 100-point support grid on [min(Z), max(Z)]
    real colvector kd0_Z        // density estimates at zeval
    real colvector kd1_Z        // density derivative estimates at zeval
}

// -----------------------------------------------------------------------------
// DidHeteroKernelConsts
// Stores precomputed kernel integral moments and derived constants.
// Initialized once per estimation call based on kernel choice.
//
// Supported kernels:
//   "epa" - Epanechnikov: K(u) = 0.75*(1-u^2) * (|u|<=1)
//   "gau" - Gaussian:     K(u) = phi(u)
//
// Integral notation:
//   I_{j,f} = int u^j * f(u) du
// -----------------------------------------------------------------------------
struct DidHeteroKernelConsts {

    // --- Kernel identifier ---
    string scalar  name        // kernel name: "epa" or "gau"

    // --- Kernel integral moments I_{j,K} = int u^j K(u) du ---
    real scalar    I_2_K       // I_{2,K}:  epa = 1/5,   gau = 1
    real scalar    I_4_K       // I_{4,K}:  epa = 3/35,  gau = 3
    real scalar    I_6_K       // I_{6,K}:  epa = 1/21,  gau = 15

    // --- Squared-kernel integral moments I_{j,K^2} = int u^j K(u)^2 du ---
    real scalar    I_0_K2      // I_{0,K^2}: epa = 3/5,   gau = 1/(2*sqrt(pi))
    real scalar    I_2_K2      // I_{2,K^2}: epa = 3/35,  gau = 1/(4*sqrt(pi))
    real scalar    I_4_K2      // I_{4,K^2}: epa = 1/35,  gau = 3/(8*sqrt(pi))
    real scalar    I_6_K2      // I_{6,K^2}: epa = 1/77,  gau = 15/(16*sqrt(pi))

    // --- Analytical UCB constant ---
    // lambda = -int(K * K'' du) / int(K^2 du)
    real scalar    lambda

    // --- Variance constants ---
    // const_V1: LLR variance constant = I_{0,K^2}
    real scalar    const_V1

    // const_V2: LQR variance constant
    real scalar    const_V2

    // --- Bias constants ---
    // const_B1: LLR bias constant = I_{2,K} / 2
    real scalar    const_B1

    // const_B2: LQR bias constant
    real scalar    const_B2
}

// -----------------------------------------------------------------------------
// BootPrecomp
// Pre-computed invariants for bootstrap loop optimization.
// Stores kernel values, design matrices, and partial matrix products
// that do not change across bootstrap iterations.
// -----------------------------------------------------------------------------
struct BootPrecomp {

    // Kernel values: num_gteval × num_zeval matrix of pointers
    // Each pointer -> n × 1 real colvector K_h(Z_i - z_r)
    pointer(real colvector) matrix kernel_vals

    // Design matrices: num_gteval × num_zeval matrix of pointers
    // Each pointer -> n × (p+1) real matrix R_{z_r}
    pointer(real matrix) matrix design_mats

    // R'diag(K)R: num_gteval × num_zeval matrix of pointers
    // Each pointer -> (p+1) × (p+1) real matrix
    pointer(real matrix) matrix RtKR

    // R'diag(K)A: num_gteval × num_zeval matrix of pointers
    // Each pointer -> (p+1) × 1 real colvector
    pointer(real colvector) matrix RtKA

    // Dimensions for validation
    real scalar num_gteval
    real scalar num_zeval
    real scalar n
    real scalar porder
}

// -----------------------------------------------------------------------------
// DidHeteroEstResult
// Stores estimation results for a single (g,t) pair.
// All vector fields are num_zeval x 1 (one entry per z evaluation point).
// -----------------------------------------------------------------------------
struct DidHeteroEstResult {

    // --- Point estimates and inference ---
    real colvector est         // num_zeval x 1 point estimates
    real colvector se          // num_zeval x 1 standard errors

    // --- Analytical uniform confidence bands ---
    real colvector ci1_L       // num_zeval x 1 analytical UCB lower bound
    real colvector ci1_U       // num_zeval x 1 analytical UCB upper bound

    // --- Bootstrap uniform confidence bands ---
    real colvector ci2_L       // num_zeval x 1 bootstrap UCB lower bound
    real colvector ci2_U       // num_zeval x 1 bootstrap UCB upper bound

    // --- Influence function matrix ---
    real matrix    B_gt_mat    // n x num_zeval influence function matrix
                               // Column r corresponds to evaluation point z_r

    // --- Bandwidth ---
    real scalar    bw          // bandwidth used for this (g,t) pair
}

// -----------------------------------------------------------------------------
// DidHeteroAggteResult
// Stores aggregated parameter estimation results.
//
// Aggregation types:
//   dynamic  - event-time aggregation (eval_pts = event times)
//   group    - group aggregation       (eval_pts = group values)
//   calendar - calendar aggregation    (eval_pts = calendar times)
//   simple   - simple aggregation      (eval_pts = single value)
//
// All matrix fields are num_zeval x num_eval.
// -----------------------------------------------------------------------------
struct DidHeteroAggteResult {

    // --- Evaluation points ---
    real colvector eval_pts    // num_eval x 1 aggregation evaluation points
                               // dynamic: event times
                               // group: group values
                               // calendar: calendar times
                               // simple: single value

    // --- Point estimates and inference ---
    real matrix    est         // num_zeval x num_eval point estimates
    real matrix    se          // num_zeval x num_eval standard errors

    // --- Analytical uniform confidence band (UCB) ---
    real matrix    ci1_L       // num_zeval x num_eval analytical UCB lower
    real matrix    ci1_U       // num_zeval x num_eval analytical UCB upper

    // --- Bootstrap uniform confidence band (UCB) ---
    real matrix    ci2_L       // num_zeval x num_eval bootstrap UCB lower
    real matrix    ci2_U       // num_zeval x num_eval bootstrap UCB upper

    // --- Bandwidth ---
    real colvector bw          // num_eval x 1 bandwidth per eval point
}

// -----------------------------------------------------------------------------
// DidHeteroCattResult
// Return structure from didhetero_catt_core() — the reentrant CATT estimator.
// Contains all outputs needed by the aggregation re-estimation pass.
//
// Dimensions:
//   K = num_gteval (number of (g,t) pairs estimated)
//   n = sample size
//   R = num_zeval (number of z evaluation points)
// -----------------------------------------------------------------------------
struct DidHeteroCattResult {
    pointer(real matrix) rowvector A_g_t    // 1 x K, each -> n x R DR score matrix
    pointer(real matrix) rowvector B_g_t    // 1 x K, each -> n x R influence function
    real matrix    G_g           // n x K group indicator matrix
    real matrix    mu_G_g        // R x K conditional group density
    real matrix    mu_E_g_t      // R x K conditional mean of E (R*Y_diff)
    real matrix    mu_F_g_t      // R x K conditional mean of F (G*Y_diff)
    real matrix    catt_est      // R x K CATT point estimates (DR estimates)
    real colvector bw_vec        // K x 1 bandwidths used
    real scalar    n             // sample size
    real scalar    num_gteval    // number of (g,t) pairs
    real scalar    num_zeval     // number of z evaluation points
}

// -----------------------------------------------------------------------------
// AggtResult
// Aggregated parameter estimation results for the aggte_gt command.
// This is the top-level return structure from didhetero_aggte_main(),
// containing results across ALL eval points (the full command output).
//
// Orientation: num_eval (rows) x num_zeval (cols) for matrix fields.
//
// Distinct from DidHeteroAggteResult which stores per-eval-point results
// with num_zeval x num_eval orientation.
// -----------------------------------------------------------------------------
struct AggtResult {
    real matrix    aggte_est    // num_eval x num_zeval: point estimates
    real matrix    aggte_se     // num_eval x num_zeval: standard errors
    real matrix    ci1_lower    // num_eval x num_zeval: analytic UCB lower
    real matrix    ci1_upper    // num_eval x num_zeval: analytic UCB upper
    real matrix    ci2_lower    // num_eval x num_zeval: bootstrap UCB lower
    real matrix    ci2_upper    // num_eval x num_zeval: bootstrap UCB upper
    real colvector aggte_bw     // num_eval x 1: aggte bandwidths
    real colvector eval_info    // num_eval x 1: eval point values
    real colvector zeval        // num_zeval x 1: z evaluation points
    string scalar  type         // aggregation type: "dynamic"|"group"|"calendar"|"simple"
}

// -----------------------------------------------------------------------------
// Version function
// -----------------------------------------------------------------------------
string scalar didhetero_version()
{
    return("0.1.0")
}

end
