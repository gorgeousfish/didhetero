*! catt_gt.ado
*! Conditional Average Treatment Effect on the Treated (CATT)
*! with continuous heterogeneity
*!
*! Core estimation command for group-time CATT(g,t,z) functions.
*!
*! Syntax:
*!   catt_gt depvar, id() time() group() z() zeval() [options]

program define catt_gt, eclass
    version 16.0

    // Parse syntax and validate options
    local _dh_raw_bstrap_found 0
    local _dh_raw_bstrap_value ""
    local _dh_raw_uniformall_found 0
    local _dh_raw_uniformall_value ""
    local _dh_rebuilt_0 ""
    local _dh_rest `"`0'"'

    while `"`_dh_rest'"' != "" {
        gettoken _dh_tok _dh_rest : _dh_rest, bind
        local _dh_tok = trim(`"`_dh_tok'"')
        if `"`_dh_tok'"' == "" {
            continue
        }

        local _dh_tok_lc = lower(`"`_dh_tok'"')
        if regexm(`"`_dh_tok_lc'"', "^bstrap[(](.*)[)]$") {
            if `_dh_raw_bstrap_found' {
                di as error "bstrap() may be specified only once"
                exit 198
            }

            local _dh_inner = trim(`"`=regexs(1)'"')
            local _dh_inner_len = length(`"`_dh_inner'"')
            if `_dh_inner_len' >= 2 & substr(`"`_dh_inner'"', 1, 1) == `"""' & substr(`"`_dh_inner'"', `_dh_inner_len', 1) == `"""' {
                local _dh_inner = substr(`"`_dh_inner'"', 2, `_dh_inner_len' - 2)
                local _dh_inner = trim(`"`_dh_inner'"')
            }

            local _dh_inner_lc = lower(`"`_dh_inner'"')
            if !inlist(`"`_dh_inner_lc'"', "true", "false") {
                di as error "bstrap() must be true or false"
                exit 198
            }

            local _dh_raw_bstrap_found 1
            local _dh_raw_bstrap_value `"`_dh_inner_lc'"'
            continue
        }

        if regexm(`"`_dh_tok_lc'"', "^uniformall[(](.*)[)]$") {
            if `_dh_raw_uniformall_found' {
                di as error "uniformall() may be specified only once"
                exit 198
            }

            local _dh_inner = trim(`"`=regexs(1)'"')
            local _dh_inner_len = length(`"`_dh_inner'"')
            if `_dh_inner_len' >= 2 & substr(`"`_dh_inner'"', 1, 1) == `"""' & substr(`"`_dh_inner'"', `_dh_inner_len', 1) == `"""' {
                local _dh_inner = substr(`"`_dh_inner'"', 2, `_dh_inner_len' - 2)
                local _dh_inner = trim(`"`_dh_inner'"')
            }

            local _dh_inner_lc = lower(`"`_dh_inner'"')
            if !inlist(`"`_dh_inner_lc'"', "true", "false") {
                di as error "uniformall() must be true or false"
                exit 198
            }

            local _dh_raw_uniformall_found 1
            local _dh_raw_uniformall_value `"`_dh_inner_lc'"'
            continue
        }

        local _dh_rebuilt_0 `"`_dh_rebuilt_0' `_dh_tok'"'
    }

    local 0 `"`_dh_rebuilt_0'"'

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        Id(varname)                                  ///
        Time(varname)                                ///
        Group(varname)                               ///
        Z(varname)                                   ///
        Zeval(numlist)                               ///
        [Xformula(string)]                           ///
        [GTeval(numlist)]                            ///
        [Porder(integer 2)]                          ///
        [Kernel(string)]                             ///
        [Control_group(string)]                      ///
        [Anticipation(integer 0)]                    ///
        [Alp(real -1)]                               ///
        [Level(cilevel)]                             ///
        [Biters(integer 1000)]                       ///
        [noBSTrap]                                   ///
        [noUNIFormall]                               ///
        [PREtrend]                                   ///
        [BWselect(string)]                           ///
        [BW(numlist)]                                ///
        [SEed(integer -1)]                           ///
        [VERBose]                                    ///
        [GPSStrict]                                  ///
        [KDETrim]                                    ///
        [GPSTrim(numlist)]                            ///
        [RBC]                                         ///
        [UNDERSmooth]                                    ///
        [PROFile]

    local depvar `varlist'

    // Resolve significance level: alp() takes priority if explicitly given;
    // otherwise compute from level() (which defaults to c(level)).
    if `alp' != -1 {
        // User explicitly specified alp(); validate range
        if `alp' <= 0 | `alp' >= 1 {
            di as error "alp() must be strictly between 0 and 1"
            exit 198
        }
    }
    else {
        // Derive from level option (cilevel defaults to c(level))
        local alp = 1 - `level'/100
    }

    // Verbose flag: default off (quiet mode)
    if "`verbose'" != "" {
        local _dh_verbose 1
    }
    else {
        local _dh_verbose 0
    }

    // GPS strict mode: error on non-convergence (default off)
    if "`gpsstrict'" != "" {
        local _dh_gps_strict 1
    }
    else {
        local _dh_gps_strict 0
    }

    // Profile flag: performance profiling (default off)
    if "`profile'" != "" {
        local _dh_profile 1
    }
    else {
        local _dh_profile 0
    }

    // Set string defaults
    if "`kernel'" == "" local kernel "gau"
    if "`control_group'" == "" local control_group "notyettreated"
    if "`bwselect'" == "" local bwselect "IMSE1"

    // Bootstrap default: ON
    if `_dh_raw_bstrap_found' & "`bstrap'" != "" {
        di as error "bstrap() cannot be combined with legacy bootstrap flags"
        exit 198
    }
    if `_dh_raw_bstrap_found' {
        if `"`_dh_raw_bstrap_value'"' == "true" {
            local bstrap "bstrap"
        }
        else {
            local bstrap ""
        }
    }
    else if "`bstrap'" == "nobstrap" {
        local bstrap ""
    }
    else {
        local bstrap "bstrap"
    }

    // Uniformall default: ON
    if `_dh_raw_uniformall_found' & "`uniformall'" != "" {
        di as error "uniformall() cannot be combined with legacy uniform flags"
        exit 198
    }
    if `_dh_raw_uniformall_found' {
        if `"`_dh_raw_uniformall_value'"' == "true" {
            local uniformall "uniformall"
        }
        else {
            local uniformall ""
        }
    }
    else if "`uniformall'" == "nouniformall" {
        local uniformall ""
    }
    else {
        local uniformall "uniformall"
    }

    // Validate seed: -1 preserves current RNG state; non-negative integers set seed
    if (`seed' < -1) {
        di as error "catt_gt: seed() must be -1 or a nonnegative integer"
        di as error "  seed(-1) leaves the current RNG state unchanged"
        di as error "  received: `seed'"
        exit 198
    }

    // Scalar bw() required when gteval() is omitted and bwselect = manual
    if "`bwselect'" == "manual" & "`bw'" != "" & "`gteval'" == "" {
        local n_bw_tokens : word count `bw'
        if `n_bw_tokens' != 1 {
            di as error "catt_gt requires scalar bw() when bwselect = 'manual' and gteval() is omitted"
            exit 198
        }
    }

    // Protect user data and apply [if] [in] restrictions
    preserve
    if "`if'`in'" != "" {
        quietly keep `if' `in'
    }

    capture noisily _dh_ensure_backend
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // Validate data and parameters
    _didhetero_validate `depvar',            ///
        id(`id')                              ///
        time(`time')                          ///
        group(`group')                        ///
        z(`z')                                ///
        zeval(`zeval')                        ///
        gteval(`gteval')                      ///
        xformula(`xformula')                  ///
        porder(`porder')                      ///
        kernel(`kernel')                      ///
        control_group(`control_group')        ///
        anticipation(`anticipation')          ///
        alp(`alp')                            ///
        biters(`biters')                      ///
        `bstrap'                              ///
        `uniformall'                          ///
        `pretrend'                            ///
        bwselect(`bwselect')                  ///
        bw(`bw')                              ///
        `kdetrim'                             ///
        gpstrim(`gpstrim')                    ///
        `rbc'                                 ///
        `undersmooth'

    // Retrieve validated parameters
    local depvar    `_dh_depvar'
    local id        `_dh_id'
    local time      `_dh_time'
    local group     `_dh_group'
    local z         `_dh_z'
    local zeval     `_dh_zeval'
    local xformula  `_dh_xformula'
    local xformula_display `_dh_xformula_display'
    local xformula_has_intercept `_dh_xformula_has_intercept'
    local porder    `_dh_porder'
    local kernel    `_dh_kernel'
    local control   `_dh_control'
    local anticip   `_dh_anticip'
    local alp       `_dh_alp'
    local biters    `_dh_biters'
    local bstrap    `_dh_bstrap'
    local uniform   `_dh_uniform'
    local pretrend  `_dh_pretrend'
    local kdetrim   `_dh_kdetrim'
    local bwselect  `_dh_bwselect'
    local bw        `_dh_bw'
    local gpstrim   `_dh_gpstrim'
    local rbc_flag  `_dh_rbc'
    local undersmooth_flag `_dh_undersmooth'
    local n_total   `_dh_n'

    // Set biters to 0 when bootstrap is disabled
    if `bstrap' == 0 {
        local biters = 0
    }

    local _dh_level = 100 * (1 - `alp')
    // Header is rendered from Mata (_didhetero_render_header) during the
    // estimation run, once the (g,t) pairs are known. Set the command name
    // so the renderer can label the output correctly.
    local _dh_cmdname "catt_gt"

    // Prepare Mata data structures and initialize kernel constants

    ereturn clear
    capture noisily mata: didhetero_init_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // Process user-specified gteval pairs
    if "`gteval'" != "" {
        local ngt_tokens : word count `gteval'
        if mod(`ngt_tokens', 2) != 0 {
            restore
            capture _est unhold `_dh_prev_est'
            di as error "gteval() must contain an even number of values (g1 t1 g2 t2 ...)"
            exit 198
        }
        local _ngt_pairs = `ngt_tokens' / 2
        local _gteval_user `gteval'

        if "`bwselect'" == "manual" & "`bw'" != "" {
            local _n_bw_tokens : word count `bw'
            if (`_n_bw_tokens' != 1) & (`_n_bw_tokens' != `_ngt_pairs') {
                restore
                capture _est unhold `_dh_prev_est'
                di as error "bw must be a positive scalar or vector whose length equals to the number of gteval."
                exit 198
            }
        }

        capture noisily mata: _didhetero_validate_user_gteval()
        local _dh_rc = _rc
        if `_dh_rc' {
            restore
            capture _est unhold `_dh_prev_est'
            exit `_dh_rc'
        }
        if "`_gteval_duplicate'" == "1" {
            restore
            capture _est unhold `_dh_prev_est'
            di as error ///
                "gteval() contains duplicate (g,t) pairs: `_gteval_duplicate_pairs'"
            exit 198
        }
        if "`_gteval_invalid'" == "1" {
            restore
            capture _est unhold `_dh_prev_est'
            di as error ///
                "gteval() contains pairs outside the identification domain implied by the observed sample, control_group(`control_group'), and anticipation(`anticipation'): `_gteval_invalid_pairs'"
            exit 198
        }

        capture noisily mata: _didhetero_set_user_gteval()
        local _dh_rc = _rc
        if `_dh_rc' {
            restore
            capture _est unhold `_dh_prev_est'
            exit `_dh_rc'
        }
    }

    // Execute estimation pipeline
    capture noisily mata: didhetero_run_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // Clean up stale bootstrap matrices when bootstrap is disabled
    if `bstrap' == 0 {
        capture matrix drop e(c_check)
    }

    // Restore original data and post results
    restore
    capture noisily _didhetero_post_eclass
    local _dh_rc = _rc
    if `_dh_rc' {
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }
    capture _est unhold `_dh_prev_est', not

    ereturn local cmd "catt_gt"
    ereturn local predict "catt_gt_predict"
    ereturn local estat_cmd "catt_gt_estat"
    ereturn local depvar     "`depvar'"
    ereturn local idvar      "`id'"
    ereturn local timevar    "`time'"
    ereturn local groupvar   "`group'"
    ereturn local zvar       "`z'"
    ereturn local kernel     "`kernel'"
    ereturn local control_group "`control'"
    ereturn local control    "`control'"
    ereturn local bwselect   "`bwselect'"
    ereturn scalar porder    = cond(`rbc_flag', 2, `porder')
    ereturn scalar anticipation = `anticip'
    ereturn scalar anticip   = `anticip'
    ereturn scalar alp       = `alp'
    ereturn scalar level     = 100 * (1 - `alp')
    ereturn scalar bstrap    = `bstrap'
    ereturn scalar biters    = `biters'
    local _dh_effective_seed = .
    if (`bstrap' == 1) & (`seed' >= 0) {
        local _dh_effective_seed = `seed'
    }
    ereturn scalar seed_request = `seed'
    ereturn scalar seed      = `_dh_effective_seed'
    if e(num_gteval) == 1 {
        local uniform 0
    }
    ereturn scalar uniformall = `uniform'
    ereturn scalar pretrend  = `pretrend'
    ereturn scalar rbc       = `rbc_flag'
    if `undersmooth_flag' == 1 {
        ereturn scalar undersmooth = 1
        if "`_dh_bw_adjusted'" != "" {
            ereturn scalar bw_adjusted = `_dh_bw_adjusted'
        }
        else {
            ereturn scalar bw_adjusted = 0
        }
    }

    // Display results table
    _didhetero_display

end
