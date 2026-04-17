*! aggte_gt.ado
*! Aggregated Group-Time Treatment Effect Estimation
*! 
*! Computes doubly robust estimates and uniform confidence bands
*! for summary parameters that aggregate group-time conditional
*! average treatment effects (CATT) given a continuous covariate.
*! 
*! Syntax:
*!   aggte_gt, [type(string) eval(numlist) bstrap(string) biters(integer)
*!              porder(integer) bwselect(string) bw(numlist)
*!              uniformall(string) seed(integer)]

program define aggte_gt, eclass
    version 16.0

    // -------------------------------------------------------------
    // Syntax parsing
    // -------------------------------------------------------------
    syntax , [TYpe(string) EVAL(numlist) ///
              BSTRap(string) BITers(integer 1000) ///
              POrder(integer 2) ///
              BWSelect(string) BW(numlist) ///
              UNIFormall(string) ///
              SEed(integer -1)]

    // -------------------------------------------------------------
    // Normalize upstream estimation results
    // -------------------------------------------------------------
    _aggte_normalize_upstream_e
    quietly _dh_ensure_backend

    // -------------------------------------------------------------
    // Parameter defaults
    // -------------------------------------------------------------
    if "`type'" == "" local type "dynamic"
    if "`bstrap'" == "" local bstrap "true"
    local kernel = e(kernel)
    if "`kernel'" == "" local kernel "gau"
    if "`bwselect'" == "" local bwselect "IMSE1"
    if "`uniformall'" == "" local uniformall "true"

    if "`kernel'" == "gaussian" local kernel "gau"
    if "`kernel'" == "epanechnikov" local kernel "epa"
    local bstrap = lower(trim("`bstrap'"))
    local uniformall = lower(trim("`uniformall'"))

    capture confirm scalar e(alp)
    if _rc {
        local alp 0.05
    }
    else {
        local alp = e(alp)
        if missing(`alp') local alp 0.05
    }

    // -------------------------------------------------------------
    // Parameter validation
    // -------------------------------------------------------------

    // Validate type
    if !inlist("`type'", "dynamic", "group", "calendar", "simple") {
        di as error "aggte_gt: type() must be one of: dynamic, group, calendar, simple"
        di as error "  received: `type'"
        exit 198
    }

    // simple aggregation has no eval dimension
    if "`type'" == "simple" & "`eval'" != "" {
        di as error "aggte_gt: eval() is not allowed when type(simple)"
        di as error "  type(simple) aggregates over all post-treatment (g,t) pairs and has no eval dimension"
        exit 198
    }

    // Validate porder
    if !inlist(`porder', 1, 2) {
        di as error "aggte_gt: porder() must be 1 or 2"
        di as error "  received: `porder'"
        exit 198
    }

    // Validate bwselect
    if !inlist("`bwselect'", "IMSE1", "IMSE2", "US1", "manual") {
        di as error "bwselect must be 'IMSE1', 'IMSE2', 'US1', or 'manual'"
        exit 198
    }

    // Validate alp
    if `alp' <= 0 | `alp' >= 1 {
        di as error "aggte_gt: alp() must be in (0, 1)"
        di as error "  received: `alp'"
        exit 198
    }

    // Validate kernel
    if !inlist("`kernel'", "gau", "epa") {
        di as error "aggte_gt: kernel() must be gau or epa"
        di as error "  long forms gaussian and epanechnikov are also accepted"
        di as error "  received: `kernel'"
        exit 198
    }

    // Validate Boolean string options
    if !inlist("`bstrap'", "true", "false") {
        di as error "aggte_gt: bstrap() must be true or false"
        di as error "  received: `bstrap'"
        exit 198
    }

    if !inlist("`uniformall'", "true", "false") {
        di as error "aggte_gt: uniformall() must be true or false"
        di as error "  received: `uniformall'"
        exit 198
    }

    // Validate seed domain
    if (`seed' < -1) {
        di as error "aggte_gt: seed() must be -1 or a nonnegative integer"
        di as error "  seed(-1) leaves the current RNG state unchanged"
        di as error "  received: `seed'"
        exit 198
    }

    // Validate biters when bootstrap is enabled
    if "`bstrap'" == "true" & `biters' <= 0 {
        di as error "When bstrap = TRUE, biters must be a positive number."
        exit 198
    }

    // -------------------------------------------------------------
    // e() precondition checks
    // -------------------------------------------------------------

    // Check e(cmd) source
    if ("`e(cmd)'" != "catt_gt") & ("`e(cmd)'" != "didhetero") & ("`e(cmd)'" != "aggte_gt") {
        di as error "aggte_gt requires catt_gt or didhetero results in e()"
        di as error "  e(cmd) = `e(cmd)'"
        di as error "  Please run catt_gt or didhetero before aggte_gt."
        exit 301
    }

    // Check required matrices exist
    local _aggte_req_matrices "B_g_t G_g Z mu_G_g gteval catt_est catt_se zeval bw kd0_Z kd1_Z Z_supp"
    foreach _mat of local _aggte_req_matrices {
        capture confirm matrix e(`_mat')
        if _rc {
            di as error "aggte_gt: required matrix e(`_mat') not found in e() results"
            di as error "  This may indicate an incomplete catt_gt run."
            di as error "  Please re-run catt_gt before aggte_gt."
            exit 198
        }
    }

    // Check required scalar e(gbar)
    capture confirm scalar e(gbar)
    if _rc {
        di as error "aggte_gt: required scalar e(gbar) not found in e() results"
        di as error "  This may indicate an incomplete catt_gt run."
        di as error "  Please re-run catt_gt before aggte_gt."
        exit 198
    }

    // Snapshot upstream e() state and clear stale results
    tempname _agg_base_results _agg_base_estimate _agg_base_estimate_b _agg_base_gteval ///
             _agg_base_zeval _agg_base_bw _agg_base_c_hat _agg_base_c_check ///
             _agg_base_B_g_t _agg_base_G_g _agg_base_Z _agg_base_dh_Y_wide ///
             _agg_base_dh_G_unit _agg_base_dh_t_vals _agg_base_dh_gps_mat ///
             _agg_base_dh_or_mat _agg_base_mu_G_g _agg_base_catt_est ///
             _agg_base_catt_se _agg_base_kd0_Z _agg_base_kd1_Z ///
             _agg_base_b _agg_base_V ///
             _agg_base_Z_supp
    matrix `_agg_base_results' = e(results)
    matrix `_agg_base_estimate' = e(results)
    capture confirm matrix e(b)
    local _has_base_b = (_rc == 0)
    if `_has_base_b' {
        matrix `_agg_base_b' = e(b)
    }
    capture confirm matrix e(V)
    local _has_base_V = (_rc == 0)
    if `_has_base_V' {
        matrix `_agg_base_V' = e(V)
    }
    capture confirm matrix e(Estimate_b)
    local _has_base_estimate_b = (_rc == 0)
    if `_has_base_estimate_b' {
        matrix `_agg_base_estimate_b' = e(Estimate_b)
    }
    matrix `_agg_base_gteval' = e(gteval)
    matrix `_agg_base_zeval' = e(zeval)
    matrix `_agg_base_bw' = e(bw)
    matrix `_agg_base_c_hat' = e(c_hat)
    capture confirm matrix e(c_check)
    local _has_base_c_check = (_rc == 0)
    if `_has_base_c_check' {
        matrix `_agg_base_c_check' = e(c_check)
    }
    matrix `_agg_base_B_g_t' = e(B_g_t)
    matrix `_agg_base_G_g' = e(G_g)
    matrix `_agg_base_Z' = e(Z)
    matrix `_agg_base_dh_Y_wide' = e(dh_Y_wide)
    matrix `_agg_base_dh_G_unit' = e(dh_G_unit)
    matrix `_agg_base_dh_t_vals' = e(dh_t_vals)
    matrix `_agg_base_dh_gps_mat' = e(dh_gps_mat)
    matrix `_agg_base_dh_or_mat' = e(dh_or_mat)
    matrix `_agg_base_mu_G_g' = e(mu_G_g)
    matrix `_agg_base_catt_est' = e(catt_est)
    matrix `_agg_base_catt_se' = e(catt_se)
    matrix `_agg_base_kd0_Z' = e(kd0_Z)
    matrix `_agg_base_kd1_Z' = e(kd1_Z)
    matrix `_agg_base_Z_supp' = e(Z_supp)

    local _agg_base_gbar = e(gbar)
    capture confirm scalar e(gbar_isinf)
    if _rc == 0 {
        local _agg_base_gbar_isinf = e(gbar_isinf)
    }
    else {
        local _agg_base_gbar_isinf = missing(e(gbar))
    }
    local _agg_base_N = e(N)
    local _agg_base_num_gteval = e(num_gteval)
    local _agg_base_num_zeval = e(num_zeval)
    local _agg_base_T = e(T)
    local _agg_base_porder = e(porder)
    capture confirm scalar e(anticipation)
    if _rc == 0 local _agg_base_anticipation = e(anticipation)
    else local _agg_base_anticipation = e(anticip)
    local _agg_base_anticip = `_agg_base_anticipation'
    local _agg_base_alp = e(alp)
    local _agg_base_bstrap = e(bstrap)
    local _agg_base_biters = e(biters)
    capture confirm scalar e(seed_request)
    if _rc == 0 local _agg_base_seed_request = e(seed_request)
    else local _agg_base_seed_request = -1
    capture confirm scalar e(seed)
    if _rc == 0 local _agg_base_seed = e(seed)
    else local _agg_base_seed = .
    local _agg_base_uniformall = e(uniformall)
    local _agg_base_pretrend = e(pretrend)

    local _agg_base_source_cmd "`e(cmd)'"
    if "`_agg_base_source_cmd'" == "aggte_gt" {
        local _agg_base_source_cmd "`e(aggte_source_cmd)'"
        if "`_agg_base_source_cmd'" == "" {
            if "`e(control)'" != "" local _agg_base_source_cmd "didhetero"
            else local _agg_base_source_cmd "catt_gt"
        }
    }
    local _agg_base_depvar "`e(depvar)'"
    local _agg_base_idvar "`e(idvar)'"
    local _agg_base_timevar "`e(timevar)'"
    local _agg_base_groupvar "`e(groupvar)'"
    local _agg_base_zvar "`e(zvar)'"
    local _agg_base_kernel "`e(kernel)'"
    local _agg_base_bwselect "`e(bwselect)'"
    local _agg_base_control_group "`e(control_group)'"
    local _agg_base_control "`e(control)'"
    if "`_agg_base_control_group'" == "" local _agg_base_control_group "`_agg_base_control'"
    if "`_agg_base_control'" == "" local _agg_base_control "`_agg_base_control_group'"

    ereturn clear
    if `_has_base_b' & `_has_base_V' {
        ereturn post `_agg_base_b' `_agg_base_V', obs(`_agg_base_N')
    }
    ereturn matrix results = `_agg_base_results'
    ereturn matrix Estimate = `_agg_base_estimate'
    if `_has_base_estimate_b' {
        ereturn matrix Estimate_b = `_agg_base_estimate_b'
    }
    ereturn matrix gteval = `_agg_base_gteval'
    ereturn matrix zeval = `_agg_base_zeval'
    ereturn matrix bw = `_agg_base_bw'
    ereturn matrix c_hat = `_agg_base_c_hat'
    if `_has_base_c_check' {
        ereturn matrix c_check = `_agg_base_c_check'
    }
    ereturn matrix B_g_t = `_agg_base_B_g_t'
    ereturn matrix G_g = `_agg_base_G_g'
    ereturn matrix Z = `_agg_base_Z'
    ereturn matrix dh_Y_wide = `_agg_base_dh_Y_wide'
    ereturn matrix dh_G_unit = `_agg_base_dh_G_unit'
    ereturn matrix dh_t_vals = `_agg_base_dh_t_vals'
    ereturn matrix dh_gps_mat = `_agg_base_dh_gps_mat'
    ereturn matrix dh_or_mat = `_agg_base_dh_or_mat'
    ereturn matrix mu_G_g = `_agg_base_mu_G_g'
    ereturn matrix catt_est = `_agg_base_catt_est'
    ereturn matrix catt_se = `_agg_base_catt_se'
    ereturn matrix kd0_Z = `_agg_base_kd0_Z'
    ereturn matrix kd1_Z = `_agg_base_kd1_Z'
    ereturn matrix Z_supp = `_agg_base_Z_supp'

    ereturn scalar gbar = `_agg_base_gbar'
    ereturn scalar gbar_isinf = `_agg_base_gbar_isinf'
    ereturn scalar N = `_agg_base_N'
    ereturn scalar num_gteval = `_agg_base_num_gteval'
    ereturn scalar num_zeval = `_agg_base_num_zeval'
    ereturn scalar T = `_agg_base_T'
    ereturn scalar porder = `_agg_base_porder'
    ereturn scalar anticipation = `_agg_base_anticipation'
    ereturn scalar anticip = `_agg_base_anticip'
    ereturn scalar alp = `_agg_base_alp'
    ereturn scalar bstrap = `_agg_base_bstrap'
    ereturn scalar biters = `_agg_base_biters'
    ereturn scalar seed_request = `_agg_base_seed_request'
    ereturn scalar seed = `_agg_base_seed'
    ereturn scalar uniformall = `_agg_base_uniformall'
    ereturn scalar pretrend = `_agg_base_pretrend'

    ereturn local cmd "`_agg_base_source_cmd'"
    ereturn local depvar "`_agg_base_depvar'"
    ereturn local idvar "`_agg_base_idvar'"
    ereturn local timevar "`_agg_base_timevar'"
    ereturn local groupvar "`_agg_base_groupvar'"
    ereturn local zvar "`_agg_base_zvar'"
    ereturn local kernel "`_agg_base_kernel'"
    ereturn local bwselect "`_agg_base_bwselect'"
    ereturn local control_group "`_agg_base_control_group'"
    ereturn local control "`_agg_base_control'"

    // Reject non-dynamic aggregation with no post-treatment pairs
    if inlist("`type'", "group", "calendar", "simple") {
        mata: st_numscalar("__aggte_post_pairs", ///
            sum((st_matrix("e(gteval)")[., 2] :>= st_matrix("e(gteval)")[., 1]) :& ///
                (st_matrix("e(gteval)")[., 2] :< .) :& (st_matrix("e(gteval)")[., 1] :< .)))
        if (`=scalar(__aggte_post_pairs)') == 0 {
            di as error "aggte_gt: no post-treatment (g,t) pairs available for type(`type')"
            exit 198
        }
    }

    // -------------------------------------------------------------
    // Build eval vector
    // -------------------------------------------------------------

    // Convert user eval numlist to matrix
    tempname tmp_eval_user
    local _has_user_eval = 0
    if "`eval'" != "" {
        local _has_user_eval = 1
        local n_eval_user : word count `eval'
        matrix `tmp_eval_user' = J(`n_eval_user', 1, .)
        local _i = 0
        foreach _v of local eval {
            local _i = `_i' + 1
            matrix `tmp_eval_user'[`_i', 1] = `_v'
        }
        mata: st_numscalar("__aggte_eval_hasdup", ///
            rows(uniqrows(st_matrix("`tmp_eval_user'"))) < rows(st_matrix("`tmp_eval_user'")))
        if (`=scalar(__aggte_eval_hasdup)') {
            di as error "aggte_gt: eval() cannot contain duplicate values"
            exit 198
        }
    }

    // Build eval vector via Mata
    tempname tmp_eval_result
    if `_has_user_eval' == 1 {
        mata: st_matrix("`tmp_eval_result'",    ///
            _didhetero_aggte_build_eval(        ///
                st_matrix("e(gteval)"),          ///
                "`type'",                        ///
                st_matrix("`tmp_eval_user'")))
    }
    else {
        mata: st_matrix("`tmp_eval_result'",    ///
            _didhetero_aggte_build_eval(        ///
                st_matrix("e(gteval)"),          ///
                "`type'",                        ///
                J(0, 1, .)))
    }

    local num_eval = rowsof(`tmp_eval_result')

    // Validate automatic bandwidth requirements
    if lower("`bwselect'") != "manual" & `_agg_base_num_zeval' < 2 {
        mata: st_numscalar("__aggte_need_pass1_autobw", ///
            _aggte_requires_pass1_bw(                    ///
                st_matrix("e(gteval)"),                  ///
                st_numscalar("e(gbar)"),                 ///
                "`type'",                                ///
                st_matrix("`tmp_eval_result'")))
        if (`=scalar(__aggte_need_pass1_autobw)') {
            di as error "aggte_gt: automatic aggregation bandwidth requires at least two zeval points; use bwselect(manual) with bw()"
            exit 198
        }
    }

    // -------------------------------------------------------------
    // Parse manual bandwidth vector
    // -------------------------------------------------------------
    tempname tmp_bw_user
    local _has_user_bw = 0
    if lower("`bwselect'") == "manual" & "`bw'" == "" {
        di as error "bw must be specified manually when bwselect = 'manual'."
        exit 198
    }
    if "`bw'" != "" {
        local _has_user_bw = 1
        local n_bw : word count `bw'
        // bw() requires bwselect(manual)
        if lower("`bwselect'") != "manual" {
            di as error "aggte_gt: bw() can only be specified when bwselect(manual)"
            exit 198
        }
        // Parse numlist to matrix
        matrix `tmp_bw_user' = J(`n_bw', 1, .)
        local _i = 0
        foreach _v of local bw {
            local _i = `_i' + 1
            matrix `tmp_bw_user'[`_i', 1] = `_v'
        }
        // Verify positive values
        mata: st_numscalar("__bw_min_pos", min(st_matrix("`tmp_bw_user'")))
        if (`=scalar(__bw_min_pos)' <= 0) {
            di as error "aggte_gt: all bw() values must be positive"
            exit 198
        }
        // Manual bandwidth: scalar or vector matching eval count
        if (`n_bw' != 1 & `n_bw' != `num_eval') {
            di as error "bw must be a positive scalar or vector whose length equals to the number of eval."
            di as error "  received length(bw)=`n_bw', eval count=`num_eval'"
            exit 198
        }
    }

    // -------------------------------------------------------------
    // Convert string options to numeric flags
    // -------------------------------------------------------------
    local _bstrap_flag = ("`bstrap'" == "true")
    local _uniformall_flag = ("`uniformall'" == "true")

    local _kernel_mata = "`kernel'"

    // -------------------------------------------------------------
    // Call Mata orchestrator
    // -------------------------------------------------------------
    local _bw_matname ""
    if `_has_user_bw' == 1 local _bw_matname "`tmp_bw_user'"

    mata: _didhetero_aggte_ado_entry(       ///
        "`type'",                            ///
        "`tmp_eval_result'",                 ///
        `_bstrap_flag',                      ///
        `biters',                            ///
        `seed',                              ///
        `porder',                            ///
        "`_kernel_mata'",                    ///
        "`bwselect'",                        ///
        "`_bw_matname'",                     ///
        `_uniformall_flag',                  ///
        `alp',                               ///
        "`_agg_base_control_group'",         ///
        `_agg_base_anticipation')

    // -------------------------------------------------------------
    // Display results table
    // -------------------------------------------------------------
    capture mata: _didhetero_aggte_display_table( ///
        "`type'", `porder', "`kernel'",           ///
        "`bwselect'", `alp',                      ///
        `_bstrap_flag', `biters',                 ///
        `_uniformall_flag')

    // -------------------------------------------------------------
    // Store results in e()
    // -------------------------------------------------------------

    // Move matrices from Stata globals to e()
    tempname _Est _est _se _ci1l _ci1u _ci2l _ci2u _bwm _evalm _zevalm

    matrix `_Est'   = __aggte_Estimate
    matrix `_est'   = __aggte_est
    matrix `_se'    = __aggte_se
    matrix `_ci1l'  = __aggte_ci1_lower
    matrix `_ci1u'  = __aggte_ci1_upper
    matrix `_ci2l'  = __aggte_ci2_lower
    matrix `_ci2u'  = __aggte_ci2_upper
    matrix `_bwm'   = __aggte_bw
    matrix `_evalm' = __aggte_eval
    matrix `_zevalm' = __aggte_zeval
    matrix colnames `_Est' = eval z est se ci1_lower ci1_upper ci2_lower ci2_upper bw
    matrix colnames `_evalm' = eval

    // Clean up global matrices
    capture matrix drop __aggte_Estimate
    capture matrix drop __aggte_est
    capture matrix drop __aggte_se
    capture matrix drop __aggte_ci1_lower
    capture matrix drop __aggte_ci1_upper
    capture matrix drop __aggte_ci2_lower
    capture matrix drop __aggte_ci2_upper
    capture matrix drop __aggte_bw
    capture matrix drop __aggte_eval
    capture matrix drop __aggte_zeval

    // Append to e()
    ereturn matrix Estimate = `_Est'
    ereturn matrix aggte_est = `_est'
    ereturn matrix aggte_se = `_se'
    ereturn matrix aggte_ci1_lower = `_ci1l'
    ereturn matrix aggte_ci1_upper = `_ci1u'
    ereturn matrix aggte_ci2_lower = `_ci2l'
    ereturn matrix aggte_ci2_upper = `_ci2u'
    ereturn matrix aggte_bw = `_bwm'
    ereturn matrix aggte_eval = `_evalm'
    ereturn matrix aggte_zeval = `_zevalm'

    _didhetero_post_aggte_eclass
    tempname _aggte_estb
    matrix `_aggte_estb' = e(b)
    ereturn matrix Estimate_b = `_aggte_estb'

    local _aggte_effective_biters = cond(`_bstrap_flag', `biters', 0)
    local _aggte_effective_uniformall = `_uniformall_flag'
    local _aggte_effective_seed = .
    if (`_bstrap_flag' == 1) & (`seed' >= 0) {
        local _aggte_effective_seed = `seed'
    }
    // Single eval point: joint inference degenerates to z-only case
    if `num_eval' == 1 {
        local _aggte_effective_uniformall = 0
    }

    ereturn local depvar "`_agg_base_depvar'"
    ereturn local idvar "`_agg_base_idvar'"
    ereturn local timevar "`_agg_base_timevar'"
    ereturn local groupvar "`_agg_base_groupvar'"
    ereturn local zvar "`_agg_base_zvar'"
    ereturn local kernel "`kernel'"
    ereturn local bwselect "`bwselect'"
    ereturn local control_group "`_agg_base_control_group'"
    ereturn local control "`_agg_base_control'"
    ereturn scalar porder = `porder'
    ereturn scalar anticipation = `_agg_base_anticipation'
    ereturn scalar anticip = `_agg_base_anticip'
    ereturn scalar alp = `alp'
    ereturn scalar bstrap = `_bstrap_flag'
    ereturn scalar biters = `_aggte_effective_biters'
    ereturn scalar seed_request = `seed'
    ereturn scalar seed = `_aggte_effective_seed'
    ereturn scalar uniformall = `_aggte_effective_uniformall'
    ereturn scalar pretrend = `_agg_base_pretrend'

    ereturn scalar aggte_porder = `porder'
    ereturn scalar aggte_bstrap = `_bstrap_flag'
    ereturn scalar aggte_uniformall = `_aggte_effective_uniformall'
    ereturn scalar aggte_alp = `alp'
    ereturn scalar aggte_biters = `_aggte_effective_biters'
    ereturn scalar aggte_seed_request = `seed'
    ereturn scalar aggte_seed = `_aggte_effective_seed'
    ereturn scalar aggte_base_porder = `_agg_base_porder'
    ereturn scalar aggte_base_bstrap = `_agg_base_bstrap'
    ereturn scalar aggte_base_uniformall = `_agg_base_uniformall'
    ereturn scalar aggte_base_alp = `_agg_base_alp'
    ereturn scalar aggte_base_biters = `_agg_base_biters'
    ereturn scalar aggte_base_seed_request = `_agg_base_seed_request'
    ereturn scalar aggte_base_seed = `_agg_base_seed'

    ereturn local aggte_type "`type'"
    ereturn local aggte_source_cmd "`_agg_base_source_cmd'"
    ereturn local aggte_kernel "`kernel'"
    ereturn local aggte_bwselect "`bwselect'"
    ereturn local aggte_base_kernel "`_agg_base_kernel'"
    ereturn local aggte_base_bwselect "`_agg_base_bwselect'"
    ereturn local type "`type'"
    ereturn local cmd "aggte_gt"

end

program define _aggte_normalize_upstream_e, eclass
    version 16.0

    // Accept base results or previous aggte_gt result with upstream object
    if ("`e(cmd)'" != "catt_gt") & ("`e(cmd)'" != "didhetero") & ("`e(cmd)'" != "aggte_gt") {
        di as error "aggte_gt requires catt_gt or didhetero results in e()"
        di as error "  e(cmd) = `e(cmd)'"
        di as error "  Please run catt_gt or didhetero before aggte_gt."
        exit 301
    }

    if "`e(cmd)'" != "aggte_gt" {
        exit
    }

    tempname _agg_init_results _agg_init_estimate _agg_init_estimate_b _agg_init_gteval ///
             _agg_init_zeval _agg_init_bw _agg_init_c_hat _agg_init_c_check ///
             _agg_init_B_g_t _agg_init_G_g _agg_init_Z _agg_init_dh_Y_wide ///
             _agg_init_dh_G_unit _agg_init_dh_t_vals _agg_init_dh_gps_mat ///
             _agg_init_dh_or_mat _agg_init_mu_G_g _agg_init_catt_est ///
             _agg_init_catt_se _agg_init_kd0_Z _agg_init_kd1_Z ///
             _agg_init_Z_supp

    matrix `_agg_init_results' = e(results)
    matrix `_agg_init_estimate' = e(results)
    capture confirm matrix e(Estimate_b)
    local _agg_init_has_estimate_b = (_rc == 0)
    if `_agg_init_has_estimate_b' {
        matrix `_agg_init_estimate_b' = e(Estimate_b)
    }
    matrix `_agg_init_gteval' = e(gteval)
    matrix `_agg_init_zeval' = e(zeval)
    matrix `_agg_init_bw' = e(bw)
    matrix `_agg_init_c_hat' = e(c_hat)
    capture confirm matrix e(c_check)
    local _agg_init_has_c_check = (_rc == 0)
    if `_agg_init_has_c_check' {
        matrix `_agg_init_c_check' = e(c_check)
    }
    matrix `_agg_init_B_g_t' = e(B_g_t)
    matrix `_agg_init_G_g' = e(G_g)
    matrix `_agg_init_Z' = e(Z)
    matrix `_agg_init_dh_Y_wide' = e(dh_Y_wide)
    matrix `_agg_init_dh_G_unit' = e(dh_G_unit)
    matrix `_agg_init_dh_t_vals' = e(dh_t_vals)
    matrix `_agg_init_dh_gps_mat' = e(dh_gps_mat)
    matrix `_agg_init_dh_or_mat' = e(dh_or_mat)
    matrix `_agg_init_mu_G_g' = e(mu_G_g)
    matrix `_agg_init_catt_est' = e(catt_est)
    matrix `_agg_init_catt_se' = e(catt_se)
    matrix `_agg_init_kd0_Z' = e(kd0_Z)
    matrix `_agg_init_kd1_Z' = e(kd1_Z)
    matrix `_agg_init_Z_supp' = e(Z_supp)

    local _agg_init_gbar = e(gbar)
    capture confirm scalar e(gbar_isinf)
    if _rc == 0 {
        local _agg_init_gbar_isinf = e(gbar_isinf)
    }
    else {
        local _agg_init_gbar_isinf = missing(e(gbar))
    }
    local _agg_init_N = e(N)
    local _agg_init_num_gteval = e(num_gteval)
    local _agg_init_num_zeval = e(num_zeval)
    local _agg_init_T = e(T)
    capture confirm scalar e(aggte_base_porder)
    if _rc == 0 local _agg_init_porder = e(aggte_base_porder)
    else local _agg_init_porder = e(porder)
    capture confirm scalar e(anticipation)
    if _rc == 0 local _agg_init_anticipation = e(anticipation)
    else local _agg_init_anticipation = e(anticip)
    local _agg_init_anticip = `_agg_init_anticipation'
    capture confirm scalar e(aggte_base_alp)
    if _rc == 0 local _agg_init_alp = e(aggte_base_alp)
    else local _agg_init_alp = e(alp)
    capture confirm scalar e(aggte_base_bstrap)
    if _rc == 0 local _agg_init_bstrap = e(aggte_base_bstrap)
    else local _agg_init_bstrap = e(bstrap)
    capture confirm scalar e(aggte_base_biters)
    if _rc == 0 local _agg_init_biters = e(aggte_base_biters)
    else local _agg_init_biters = e(biters)
    capture confirm scalar e(aggte_base_seed_request)
    if _rc == 0 local _agg_init_seed_request = e(aggte_base_seed_request)
    else {
        capture confirm scalar e(seed_request)
        if _rc == 0 local _agg_init_seed_request = e(seed_request)
        else local _agg_init_seed_request = -1
    }
    capture confirm scalar e(aggte_base_seed)
    if _rc == 0 local _agg_init_seed = e(aggte_base_seed)
    else {
        capture confirm scalar e(seed)
        if _rc == 0 local _agg_init_seed = e(seed)
        else local _agg_init_seed = .
    }
    capture confirm scalar e(aggte_base_uniformall)
    if _rc == 0 local _agg_init_uniformall = e(aggte_base_uniformall)
    else local _agg_init_uniformall = e(uniformall)
    local _agg_init_pretrend = e(pretrend)

    local _agg_init_source_cmd "`e(aggte_source_cmd)'"
    if "`_agg_init_source_cmd'" == "" {
        if "`e(control)'" != "" local _agg_init_source_cmd "didhetero"
        else local _agg_init_source_cmd "catt_gt"
    }
    local _agg_init_depvar "`e(depvar)'"
    local _agg_init_idvar "`e(idvar)'"
    local _agg_init_timevar "`e(timevar)'"
    local _agg_init_groupvar "`e(groupvar)'"
    local _agg_init_zvar "`e(zvar)'"
    local _agg_init_kernel "`e(aggte_base_kernel)'"
    if "`_agg_init_kernel'" == "" local _agg_init_kernel "`e(kernel)'"
    local _agg_init_bwselect "`e(aggte_base_bwselect)'"
    if "`_agg_init_bwselect'" == "" local _agg_init_bwselect "`e(bwselect)'"
    local _agg_init_control_group "`e(control_group)'"
    local _agg_init_control "`e(control)'"
    if "`_agg_init_control_group'" == "" local _agg_init_control_group "`_agg_init_control'"
    if "`_agg_init_control'" == "" local _agg_init_control "`_agg_init_control_group'"

    ereturn clear
    ereturn matrix results = `_agg_init_results'
    ereturn matrix Estimate = `_agg_init_estimate'
    if `_agg_init_has_estimate_b' {
        ereturn matrix Estimate_b = `_agg_init_estimate_b'
    }
    ereturn matrix gteval = `_agg_init_gteval'
    ereturn matrix zeval = `_agg_init_zeval'
    ereturn matrix bw = `_agg_init_bw'
    ereturn matrix c_hat = `_agg_init_c_hat'
    if `_agg_init_has_c_check' {
        ereturn matrix c_check = `_agg_init_c_check'
    }
    ereturn matrix B_g_t = `_agg_init_B_g_t'
    ereturn matrix G_g = `_agg_init_G_g'
    ereturn matrix Z = `_agg_init_Z'
    ereturn matrix dh_Y_wide = `_agg_init_dh_Y_wide'
    ereturn matrix dh_G_unit = `_agg_init_dh_G_unit'
    ereturn matrix dh_t_vals = `_agg_init_dh_t_vals'
    ereturn matrix dh_gps_mat = `_agg_init_dh_gps_mat'
    ereturn matrix dh_or_mat = `_agg_init_dh_or_mat'
    ereturn matrix mu_G_g = `_agg_init_mu_G_g'
    ereturn matrix catt_est = `_agg_init_catt_est'
    ereturn matrix catt_se = `_agg_init_catt_se'
    ereturn matrix kd0_Z = `_agg_init_kd0_Z'
    ereturn matrix kd1_Z = `_agg_init_kd1_Z'
    ereturn matrix Z_supp = `_agg_init_Z_supp'

    ereturn scalar gbar = `_agg_init_gbar'
    ereturn scalar gbar_isinf = `_agg_init_gbar_isinf'
    ereturn scalar N = `_agg_init_N'
    ereturn scalar num_gteval = `_agg_init_num_gteval'
    ereturn scalar num_zeval = `_agg_init_num_zeval'
    ereturn scalar T = `_agg_init_T'
    ereturn scalar porder = `_agg_init_porder'
    ereturn scalar anticipation = `_agg_init_anticipation'
    ereturn scalar anticip = `_agg_init_anticip'
    ereturn scalar alp = `_agg_init_alp'
    ereturn scalar bstrap = `_agg_init_bstrap'
    ereturn scalar biters = `_agg_init_biters'
    ereturn scalar seed_request = `_agg_init_seed_request'
    ereturn scalar seed = `_agg_init_seed'
    ereturn scalar uniformall = `_agg_init_uniformall'
    ereturn scalar pretrend = `_agg_init_pretrend'

    ereturn local cmd "`_agg_init_source_cmd'"
    ereturn local depvar "`_agg_init_depvar'"
    ereturn local idvar "`_agg_init_idvar'"
    ereturn local timevar "`_agg_init_timevar'"
    ereturn local groupvar "`_agg_init_groupvar'"
    ereturn local zvar "`_agg_init_zvar'"
    ereturn local kernel "`_agg_init_kernel'"
    ereturn local bwselect "`_agg_init_bwselect'"
    ereturn local control_group "`_agg_init_control_group'"
    ereturn local control "`_agg_init_control'"
end
