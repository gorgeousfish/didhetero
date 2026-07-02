mata:

// =============================================================================
// didhetero_utils_init.mata
// Data initialization utility functions for didhetero-stata
//
// Provides:
//   - didhetero_gen_z_supp()         // Generate 100-point support grid
//   - didhetero_init_arrays()        // Initialize core estimation arrays
//   - didhetero_get_uniformall()     // Return uniformall flag from global struct
//   - _didhetero_read_panel_id()     // Read Stata panel id into numeric vector
//   - didhetero_prepare_data()       // Main data preparation (Stata -> struct)
//   - didhetero_init_core_arrays()   // Initialize core arrays variant
//   - didhetero_init_from_ado()      // Initialize data from ado-layer locals
//   - didhetero_init_param_results() // Initialize empty ParamResults struct
//   - didhetero_init_stage1_results()// Initialize empty Stage1Results struct
//
// Requires:
//   - didhetero_types.mata              (DidHeteroData, DidHeteroParamResults,
//                                        DidHeteroStage1Results)
//   - didhetero_utils_formula.mata      (didhetero_selectindex,
//                                        didhetero_parse_xformula_locals)
//   - didhetero_utils_numerical.mata    (_didhetero_intersect, _didhetero_setdiff)
//   - didhetero_utils_domain.mata       (didhetero_build_gteval)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_gen_z_supp()
// Generate uniform support grid on [min(Z), max(Z)] with 100 points
//
// Args:
//   Z - real colvector of continuous treatment values
//
// Returns:
//   real colvector of 100 equally-spaced points from min(Z) to max(Z)
// -----------------------------------------------------------------------------
real colvector didhetero_gen_z_supp(real colvector Z)
{
    return(rangen(min(Z), max(Z), 100))
}

// -----------------------------------------------------------------------------
// didhetero_init_arrays()
// Initialize the 6 core estimation arrays in DidHeteroData.
// Called after num_gteval is determined.
//
// Args:
//   d          - DidHeteroData struct (modified in-place)
//   n          - number of cross-sectional units
//   num_zeval  - number of Z evaluation points
//   num_gteval - number of (g,t) evaluation pairs
//
// Side effects:
//   Populates d.A_g_t, d.B_g_t, d.G_g, d.mu_G_g, d.mu_E_g_t, d.mu_F_g_t
//   A_g_t, B_g_t are pointer arrays (1 x num_gteval), each n x num_zeval
//   G_g is n x num_gteval
//   mu_G_g, mu_E_g_t, mu_F_g_t are num_zeval x num_gteval
// -----------------------------------------------------------------------------
void didhetero_init_arrays(struct DidHeteroData scalar d,
                           real scalar n,
                           real scalar num_zeval,
                           real scalar num_gteval)
{
    real scalar k

    // A_g_t, B_g_t: pointer rowvector, each element is n x num_zeval matrix of missing
    d.A_g_t = J(1, num_gteval, NULL)
    d.B_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.A_g_t[k] = &(J(n, num_zeval, .))
        d.B_g_t[k] = &(J(n, num_zeval, .))
    }

    // G_g: n x num_gteval group indicator matrix
    d.G_g = J(n, num_gteval, .)

    // mu_G_g: num_zeval x num_gteval conditional group density
    d.mu_G_g = J(num_zeval, num_gteval, .)

    // mu_E_g_t, mu_F_g_t: num_zeval x num_gteval
    d.mu_E_g_t = J(num_zeval, num_gteval, .)
    d.mu_F_g_t = J(num_zeval, num_gteval, .)
}

// -----------------------------------------------------------------------------
// didhetero_get_uniformall()
// Return the effective uniformall flag stored in the global data struct.
// -----------------------------------------------------------------------------
real scalar didhetero_get_uniformall()
{
    external struct DidHeteroData scalar _dh_data
    return(_dh_data.uniformall)
}

// -----------------------------------------------------------------------------
// _didhetero_read_panel_id()
// Read a Stata panel id variable into a numeric column vector.
//
// Numeric ids are read directly. String ids are converted into stable
// consecutive numeric codes using the current row order. The caller guarantees
// the data are pre-sorted by (id, time), so equal ids are contiguous and this
// scan preserves panel blocks exactly.
// -----------------------------------------------------------------------------
real colvector _didhetero_read_panel_id(string scalar idvar)
{
    real scalar n_obs, i, id_code
    real colvector id_long
    string matrix id_str

    if (!st_isstrvar(idvar)) {
        return(st_data(., idvar))
    }

    id_str = st_sdata(., idvar)
    n_obs = rows(id_str)
    id_long = J(n_obs, 1, .)

    if (n_obs == 0) {
        return(id_long)
    }

    id_code = 1
    id_long[1] = id_code

    for (i = 2; i <= n_obs; i++) {
        if (id_str[i, 1] != id_str[i - 1, 1]) {
            id_code++
        }
        id_long[i] = id_code
    }

    return(id_long)
}

// -----------------------------------------------------------------------------
// didhetero_prepare_data()
// Main data preparation: read Stata data, reshape, and fill DidHeteroData struct
//
// Args:
//   depvar   - string scalar, name of the dependent variable
//   idvar    - string scalar, name of the panel id variable
//   timevar  - string scalar, name of the time variable
//   groupvar - string scalar, name of the group variable
//   zvar     - string scalar, name of the continuous treatment variable
//   zeval    - real colvector, evaluation points for Z
//   xformula - string scalar, space-separated covariate names (empty = none)
//
// Returns:
//   struct DidHeteroData scalar with all fields populated
//
// Note: Data must be pre-sorted by (id, time) via _didhetero_validate.ado
// -----------------------------------------------------------------------------
struct DidHeteroData scalar didhetero_prepare_data(
    string scalar depvar,
    string scalar idvar,
    string scalar timevar,
    string scalar groupvar,
    string scalar zvar,
    real colvector zeval,
    string scalar xformula)
{
    struct DidHeteroData scalar d
    real colvector Y_long, id_long, time_long, G_long, Z_long
    real colvector id_unique, t_vals, G, Z, id, Z_supp
    real colvector period1_idx, x_long
    real matrix Y_wide, X, X_sub
    real scalar n, T_num, period1, num_zeval, i, row_start, row_end, j, has_intercept
    string rowvector tok
    string scalar xformula_has_intercept

    // === Step 1: Read data from Stata ===
    // Data already sorted by _didhetero_validate.ado
    Y_long    = st_data(., depvar)
    id_long   = _didhetero_read_panel_id(idvar)
    time_long = st_data(., timevar)
    G_long    = st_data(., groupvar)
    Z_long    = st_data(., zvar)

    // === Step 2: Basic dimensions ===
    // Number of cross-sectional units
    id_unique = uniqrows(id_long)
    n = rows(id_unique)

    // Sorted unique time periods
    t_vals = sort(uniqrows(time_long), 1)
    T_num  = rows(t_vals)

    // First time period
    period1 = t_vals[1]

    // === Step 3: Long -> Wide Y_wide ===
    // Data is sorted by (id, time), so each block of T_num rows = one unit
    Y_wide = J(n, T_num, .)
    for (i = 1; i <= n; i++) {
        row_start = (i - 1) * T_num + 1
        row_end   = i * T_num
        Y_wide[i, .] = Y_long[row_start..row_end]'
    }

    // === Step 4: G vector (n x 1) ===
    // Extract one value per id (time-invariant, take first period)
    // Indices: 1, 1+T_num, 1+2*T_num, ..., 1+(n-1)*T_num
    // range(1, n*T_num, T_num) rounds up point count; use explicit arithmetic
    period1_idx = (0::(n-1)) :* T_num :+ 1
    G = G_long[period1_idx]

    // === Step 5: Z vector (n x 1) ===
    Z = Z_long[period1_idx]

    // === Step 6: id vector (n x 1) ===
    id = id_unique

    // === Step 7: zeval sort and range validation ===
    zeval = sort(zeval, 1)
    num_zeval = rows(zeval)
    if (min(zeval) < min(Z) | max(zeval) > max(Z)) {
        _error("zeval values must be within the range of Z: [" +
               strofreal(min(Z)) + ", " + strofreal(max(Z)) + "]")
    }

    // === Step 8: Z_supp generation ===
    Z_supp = didhetero_gen_z_supp(Z)

    // === Step 9: Covariate matrix X ===
    // The default design matrix is [intercept, formula_vars].
    // Standard formula syntax can suppress the intercept via
    // `~ 0 + ...` or `~ ... - 1`.
    xformula_has_intercept = st_local("xformula_has_intercept")
    if (xformula_has_intercept == "") has_intercept = 1
    else {
        has_intercept = strtoreal(xformula_has_intercept)
        if (missing(has_intercept)) has_intercept = 1
    }

    if (xformula == "") {
        // No extra covariates: [intercept, Z] -> n x 2
        X = J(n, 1, 1), Z
    }
    else {
        // Collapse duplicate formula terms into unique columns before
        // extracting first-period covariates.
        tok = didhetero_unique_tokens(xformula)
        X_sub = J(n, cols(tok), .)
        for (j = 1; j <= cols(tok); j++) {
            x_long = st_data(., tok[j])
            X_sub[., j] = x_long[period1_idx]
        }
        // Do NOT prepend Z separately — it is included via xformula if needed.
        // Explicit no-intercept formulas keep only the requested covariate
        // columns.
        if (has_intercept) X = J(n, 1, 1), X_sub
        else X = X_sub
    }

    // === Step 10: Fill struct ===
    d.n         = n
    d.T_num     = T_num
    d.Y_wide    = Y_wide
    d.G         = G
    d.Z         = Z
    d.id        = id
    d.t_vals    = t_vals
    d.period1   = period1
    d.zeval     = zeval
    d.num_zeval = num_zeval
    d.Z_supp    = Z_supp
    d.X         = X

    return(d)
}

// -----------------------------------------------------------------------------
// didhetero_init_core_arrays()
// Initialize core estimation arrays after num_gteval is known.
// Called after gteval computation.
//
// Args:
//   d          - DidHeteroData struct (modified in place via pointer)
//   num_gteval - number of (g,t) evaluation pairs
//
// Side effects:
//   Populates d.A_g_t, d.B_g_t, d.G_g, d.mu_G_g, d.mu_E_g_t, d.mu_F_g_t
// -----------------------------------------------------------------------------
void didhetero_init_core_arrays(struct DidHeteroData scalar d,
                                real scalar num_gteval)
{
    real scalar k

    // A_g_t: 1 x num_gteval pointer array, each -> n x num_zeval DR score matrix
    d.A_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.A_g_t[k] = &(J(d.n, d.num_zeval, .))
    }

    // B_g_t: 1 x num_gteval pointer array, each -> n x num_zeval influence function
    d.B_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.B_g_t[k] = &(J(d.n, d.num_zeval, .))
    }

    // G_g: n x num_gteval group indicator matrix
    d.G_g = J(d.n, num_gteval, .)

    // mu_G_g: num_zeval x num_gteval conditional group density
    d.mu_G_g = J(d.num_zeval, num_gteval, .)

    // mu_E_g_t, mu_F_g_t: num_zeval x num_gteval
    d.mu_E_g_t = J(d.num_zeval, num_gteval, .)
    d.mu_F_g_t = J(d.num_zeval, num_gteval, .)
}

// -----------------------------------------------------------------------------
// didhetero_init_from_ado()
// Wrapper called from catt_gt.ado to initialize data and kernel structs.
// Reads Stata locals set by catt_gt.ado, calls prepare_data and kernel_consts,
// and stores results as Mata external globals.
//
// Side effects:
//   Creates external globals _dh_data (DidHeteroData) and _dh_kc (DidHeteroKernelConsts)
// -----------------------------------------------------------------------------
void didhetero_init_from_ado()
{
    external struct DidHeteroData scalar _dh_data
    external struct DidHeteroKernelConsts scalar _dh_kc
    
    string rowvector _zeval_tokens
    real rowvector _zeval_vec
    real colvector _zeval_col
    real scalar _i, _porder

    // Parse zeval from Stata local into Mata vector
    _zeval_tokens = tokens(st_local("zeval"))
    _zeval_vec = J(1, cols(_zeval_tokens), .)
    for (_i = 1; _i <= cols(_zeval_tokens); _i++) {
        _zeval_vec[_i] = strtoreal(_zeval_tokens[_i])
    }
    _zeval_col = _zeval_vec'
    
    // Call data preparation function
    _dh_data = didhetero_prepare_data(  
        st_local("depvar"),             
        st_local("id"),                 
        st_local("time"),               
        st_local("group"),              
        st_local("z"),                  
        _zeval_col,                     
        st_local("xformula")            
    )
    
    // Store configuration parameters into struct
    _dh_data.control_group = st_local("control")
    _dh_data.anticipation  = strtoreal(st_local("anticip"))
    _dh_data.porder        = strtoreal(st_local("porder"))
    _dh_data.kernel        = st_local("kernel")
    _dh_data.alp           = strtoreal(st_local("alp"))
    _dh_data.biters        = strtoreal(st_local("biters"))
    _dh_data.uniformall    = strtoreal(st_local("uniform"))
    
    // Initialize kernel constants
    _dh_kc = didhetero_kernel_consts(st_local("kernel"))
    
    // Copy key derived constants to DidHeteroData struct
    _dh_data.const_B1 = _dh_kc.const_B1
    _dh_data.const_B2 = _dh_kc.const_B2
    _dh_data.lambda   = _dh_kc.lambda
    
    // Select const_V based on porder
    _porder = strtoreal(st_local("porder"))
    _dh_data.const_V = (_porder == 1) * _dh_kc.const_V1 + (_porder == 2) * _dh_kc.const_V2

    // Store RBC flag (0 by default; overridden in didhetero_run_from_ado if set)
    _dh_data.rbc = 0
}

// ---------------------------------------------------------------------------
// didhetero_init_param_results()
// Factory function for DidHeteroParamResults struct.
// ---------------------------------------------------------------------------
struct DidHeteroParamResults scalar didhetero_init_param_results()
{
    struct DidHeteroParamResults scalar r
    r.gps_mat  = J(0, 0, .)
    r.gps_coef = J(0, 0, .)
    r.or_mat   = J(0, 0, .)
    r.or_coef  = J(0, 0, .)
    r.ctrl_type = ""
    return(r)
}

// ---------------------------------------------------------------------------
// didhetero_init_stage1_results()
// Factory function for DidHeteroStage1Results struct.
// ---------------------------------------------------------------------------
struct DidHeteroStage1Results scalar didhetero_init_stage1_results()
{
    struct DidHeteroStage1Results scalar r
    r.gps_mat  = J(0, 0, .)
    r.gps_coef = J(0, 0, .)
    r.or_mat   = J(0, 0, .)
    r.or_coef  = J(0, 0, .)
    r.ctrl_type = ""
    r.Z_supp   = J(0, 1, .)
    r.kd0_Z    = J(0, 1, .)
    r.kd1_Z    = J(0, 1, .)
    r.kde_trimmed = J(0, 1, .)
    return(r)
}

end
