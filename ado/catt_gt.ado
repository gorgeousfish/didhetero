*! catt_gt.ado
*! Version 0.1.0
*! Conditional Average Treatment Effect on the Treated (CATT)
*! with continuous heterogeneity
*!
*! Core estimation command for group-time CATT(g,t,z) functions.
*! Implements the three-step estimation procedure described in
*! Section 4 of Imai, Qin, and Yanagi (2025).
*!
*! Syntax:
*!   catt_gt depvar, id() time() group() z() zeval() [options]

program define catt_gt, eclass
    version 16.0

    // =========================================================================
    // Step 1: Parse syntax and call validation
    // =========================================================================
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
        [Alp(real 0.05)]                             ///
        [Biters(integer 1000)]                       ///
        [noBSTrap]                                   ///
        [noUNIFormall]                               ///
        [PREtrend]                                   ///
        [BWselect(string)]                           ///
        [BW(numlist)]                                ///
        [SEed(integer -1)]

    local depvar `varlist'

    // Set string defaults
    if "`kernel'" == "" local kernel "gau"
    if "`control_group'" == "" local control_group "notyettreated"
    if "`bwselect'" == "" local bwselect "IMSE1"

    // Bootstrap default: ON (paper default)
    // [noBSTrap] syntax: "nobstrap" → OFF, "bstrap" or "" → ON
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

    // Uniformall default: ON (paper default for joint uniform inference)
    // [noUNIFormall] syntax accepts legacy flag aliases uniformall/nouniformall.
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

    // catt_gt.sthlp documents seed(-1) as the only sentinel meaning "use the
    // current RNG state". Reject all other negative integers at the API
    // boundary so invalid seeds cannot silently slip through as metadata.
    if (`seed' < -1) {
        di as error "catt_gt: seed() must be -1 or a nonnegative integer"
        di as error "  seed(-1) leaves the current RNG state unchanged"
        di as error "  received: `seed'"
        exit 198
    }

    // When gteval() is omitted, catt_gt auto-builds the full (g,t) evaluation
    // set internally, so only a scalar manual bandwidth is a well-defined API
    // contract. Once users supply explicit gteval(), vector bw() is allowed
    // and validated against the number of (g,t) rows below.
    if "`bwselect'" == "manual" & "`bw'" != "" & "`gteval'" == "" {
        local n_bw_tokens : word count `bw'
        if `n_bw_tokens' != 1 {
            di as error "catt_gt requires scalar bw() when bwselect = 'manual' and gteval() is omitted"
            exit 198
        }
    }

    // =========================================================================
    // Step 2: Protect user data and honor if/in
    // =========================================================================
    // Protect user dataset from any modifications made during validation
    // and estimation (dropping, recasting, sorting).
    preserve
    // If the user specified [if] and/or [in], restrict the working sample
    // before validation so all checks and drops are scoped to the intended
    // subset. This mirrors the behavior of didhetero.ado.
    if "`if'`in'" != "" {
        quietly keep `if' `in'
    }

    capture noisily _dh_ensure_backend
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        exit `_dh_rc'
    }

    // =========================================================================
    // Step 3: Call validation subroutine on the (possibly restricted) sample
    // =========================================================================
    // Validation must operate on the already restricted working sample.
    // Reapplying the original [if] [in] here would double-filter [in].
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
        bw(`bw')

    // Retrieve validated parameters from _didhetero_validate
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
    local bwselect  `_dh_bwselect'
    local bw        `_dh_bw'
    local n_total   `_dh_n'

    // Override biters to 0 when bootstrap is off (matching didhetero.ado)
    if `bstrap' == 0 {
        local biters = 0
    }

    di as text ""
    di as text "catt_gt: Conditional ATT with continuous heterogeneity"
    di as text "  Observations (n): `n_total'"
    di as text "  Polynomial order: `porder'"
    di as text "  Kernel:           `kernel'"
    di as text "  Control group:    `control'"
    di as text "  Anticipation:     `anticip'"
    di as text "  Significance:     `alp'"
    if "`xformula_display'" != "" {
        di as text "  Covariate spec:   `xformula_display'"
    }
    di as text ""

    // =========================================================================
    // Step 4: Call Mata data preparation and initialize kernel constants
    // =========================================================================
    // Hold the caller's last estimates so any failure after ereturn clear can
    // restore the exact pre-call state without leaving stale partial results.
    tempname _dh_prev_est
    capture _est hold `_dh_prev_est', restore nullok

    ereturn clear
    // Note: didhetero_init_from_ado() is a compiled Mata function that handles
    // external struct declarations (not supported in inline mata blocks).
    // It reads Stata locals and creates external globals _dh_data and _dh_kc.
    capture noisily mata: didhetero_init_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // =========================================================================
    // Step 5: Handle user-specified gteval
    // =========================================================================
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

    // =========================================================================
    // Step 6: Run full estimation pipeline
    // =========================================================================
    // didhetero_run_from_ado() reads the following locals:
    //   pretrend, uniform, bstrap, bwselect, bw
    // These are already set from the validated parameters above.
    capture noisily mata: didhetero_run_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // Bootstrap-off runs do not post e(c_check). Drop any stale inherited
    // matrix before _didhetero_post_eclass snapshots the successful result.
    if `bstrap' == 0 {
        capture matrix drop e(c_check)
    }

    // Restore the original user dataset. e() results produced by Mata and
    // the locals we set below are unaffected by restore.
    restore
    capture noisily _didhetero_post_eclass
    local _dh_rc = _rc
    if `_dh_rc' {
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }
    capture _est unhold `_dh_prev_est', not

    ereturn local cmd "catt_gt"
    ereturn local depvar     "`depvar'"
    ereturn local idvar      "`id'"
    ereturn local timevar    "`time'"
    ereturn local groupvar   "`group'"
    ereturn local zvar       "`z'"
    ereturn local kernel     "`kernel'"
    ereturn local control_group "`control'"
    ereturn local control    "`control'"
    ereturn local bwselect   "`bwselect'"
    ereturn scalar porder    = `porder'
    ereturn scalar anticipation = `anticip'
    ereturn scalar anticip   = `anticip'
    ereturn scalar alp       = `alp'
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

    // =========================================================================
    // Step 7: Display results table
    // =========================================================================
    _didhetero_display

end
