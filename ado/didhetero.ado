*! didhetero.ado
*!
*! Doubly Robust Uniform Confidence Bands for Group-Time Conditional
*! Average Treatment Effects in Difference-in-Differences
*!
*! Syntax:
*!   didhetero depvar, id(varname) time(varname) group(varname) z(varname)
*!         zeval(numlist) [options]
*!
*! Required options:
*!   id(varname)        - panel unit identifier
*!   time(varname)      - time period variable
*!   group(varname)     - treatment group variable (0 = never-treated)
*!   z(varname)         - continuous covariate for heterogeneity
*!   zeval(numlist)     - evaluation points for z
*!
*! Optional:
*!   xformula(string)   - first-stage covariate specification
*!   gteval(numlist)    - (g,t) evaluation pairs (g1 t1 g2 t2 ...)
*!   porder(integer 2)  - polynomial order (1 or 2)
*!   kernel(string)     - kernel function: "gau" (default) or "epa"
*!   control_group(string) - "notyettreated" (default) or "nevertreated"
*!   anticipation(integer 0) - anticipation periods
*!   alp(real 0.05)     - significance level
*!   biters(integer 1000) - bootstrap iterations
*!   bstrap(true|false) - bootstrap toggle
*!   uniformall(true|false) - uniform inference toggle
*!   pretrend           - include pre-treatment periods
*!   bwselect(string)   - bandwidth selection: "IMSE1", "IMSE2", "US1", "manual"
*!   bw(numlist)        - manual bandwidth(s)
*!   seed(integer -1)   - random number generator seed
*!
*! Returns (eclass):
*!   e(results)       - estimation results matrix
*!   e(Estimate_b)    - vectorized point estimates
*!   e(gteval)        - (g,t) evaluation pairs
*!   e(zeval)         - z evaluation points
*!   e(bw)            - bandwidths per (g,t) pair
*!   e(c_hat)         - analytical critical values
*!   e(c_check)       - bootstrap critical values
*!   e(N)             - number of cross-sectional units
*!   e(T)             - number of time periods
*!   e(num_gteval)    - number of (g,t) pairs
*!   e(num_zeval)     - number of z evaluation points
*!   e(anticipation)  - anticipation periods
*!   e(control_group) - control group specification

program define didhetero, eclass
    version 16.0

    // Parse syntax
    local _dh_raw_bstrap_found 0
    local _dh_raw_bstrap_value ""
    local _dh_raw_uniformall_found 0
    local _dh_raw_uniformall_value ""
    local _dh_rebuilt_0 ""
    local _dh_rest `"`0'"'
    local _dh_in_options 0

    while `"`_dh_rest'"' != "" {
        gettoken _dh_tok _dh_rest : _dh_rest, bind
        local _dh_tok = trim(`"`_dh_tok'"')
        if `"`_dh_tok'"' == "" {
            continue
        }

        local _dh_tok_lc = lower(`"`_dh_tok'"')
        if `_dh_in_options' & `"`_dh_tok_lc'"' == "nobstrap" {
            di as err "didhetero does not allow nobstrap; use nobootstrap or bstrap(false)"
            exit 198
        }
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
        if !`_dh_in_options' & substr(`"`_dh_tok_lc'"', -1, 1) == "," {
            local _dh_in_options 1
        }
    }

    local 0 `"`_dh_rebuilt_0'"'

    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        Id(varname)                                  ///
        Time(varname)                                ///
        Group(varname)                               ///
        Z(varname)                                   ///
        Zeval(numlist)                               ///
        [Xformula(string)]                           ///
        [GTeval(numlist)]                             ///
        [Porder(integer 2)]                          ///
        [Kernel(string)]                             ///
        [Control_group(string)]                      ///
        [Anticipation(integer 0)]                    ///
        [Alp(real 0.05)]                             ///
        [Biters(integer 1000)]                       ///
        [BSTrap]                                     ///
        [noBOOTstrap]                                ///
        [noUNIFormall]                               ///
        [PREtrend]                                   ///
        [BWselect(string)]                           ///
        [BW(numlist)]                                ///
        [SEed(integer -1)]

    local depvar `varlist'

    // Resolve option defaults and conflicts

    // Validate bootstrap specifications.
    if `_dh_raw_bstrap_found' & ("`bstrap'" != "" | "`bootstrap'" != "") {
        di as err "bstrap() cannot be combined with legacy bootstrap flags"
        exit 198
    }
    // Reject legacy nobstrap flag.
    if "`bstrap'" == "nobstrap" {
        di as err "didhetero does not allow nobstrap; use nobootstrap or bstrap(false)"
        exit 198
    }
    if "`bstrap'" != "" & "`bootstrap'" == "nobootstrap" {
        di as err "bstrap and nobootstrap cannot be specified together"
        exit 198
    }

    // Set string defaults.
    if "`kernel'" == "" local kernel "gau"
    if "`control_group'" == "" local control_group "notyettreated"
    if "`bwselect'" == "" local bwselect "IMSE1"

    // Significance level default is 0.05.

    // Bootstrap default is ON.
    if `_dh_raw_bstrap_found' {
        local bstrap_flag = (`"`_dh_raw_bstrap_value'"' == "true")
    }
    else if "`bootstrap'" == "nobootstrap" {
        local bstrap_flag = 0
    }
    else if "`bstrap'" != "" {
        local bstrap_flag = 1
    }
    else {
        local bstrap_flag = 1
    }

    // Validate seed: -1 leaves RNG unchanged; other negative values rejected.
    if (`seed' < -1) {
        di as err "didhetero: seed() must be -1 or a nonnegative integer"
        di as err "  seed(-1) leaves the current RNG state unchanged"
        di as err "  received: `seed'"
        exit 198
    }

    // Uniform inference default is ON.
    if `_dh_raw_uniformall_found' & "`uniformall'" != "" {
        di as err "uniformall() cannot be combined with legacy uniform flags"
        exit 198
    }
    if `_dh_raw_uniformall_found' {
        if `"`_dh_raw_uniformall_value'"' == "true" {
            local uniformall "uniformall"
            local uniformall_flag = 1
        }
        else {
            local uniformall ""
            local uniformall_flag = 0
        }
    }
    else if "`uniformall'" == "nouniformall" {
        local uniformall ""
        local uniformall_flag = 0
    }
    else {
        local uniformall "uniformall"
        local uniformall_flag = 1
    }
    local pretrend_flag = ("`pretrend'" != "")

    // Preserve data and apply if/in
    preserve

    if "`if'`in'" != "" {
        quietly keep `if' `in'
    }

    capture noisily _dh_ensure_backend
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        exit `_dh_rc'
    }

    // Validate input parameters
    local _dh_validate_bstrap
    if `bstrap_flag' == 1 {
        local _dh_validate_bstrap "bstrap"
    }

    _didhetero_validate `depvar',                 ///
        id(`id')                                   ///
        time(`time')                               ///
        group(`group')                             ///
        z(`z')                                     ///
        zeval(`zeval')                             ///
        gteval(`gteval')                           ///
        xformula(`xformula')                       ///
        porder(`porder')                           ///
        kernel(`kernel')                           ///
        control_group(`control_group')             ///
        anticipation(`anticipation')               ///
        alp(`alp')                                 ///
        biters(`biters')                           ///
        `_dh_validate_bstrap'                      ///
        `uniformall'                               ///
        `pretrend'                                 ///
        bwselect(`bwselect')                       ///
        bw(`bw')

    // Retrieve validated parameters.
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
    local bstrap    `bstrap_flag'
    local uniform   `uniformall_flag'
    local pretrend  `pretrend_flag'
    local bwselect  `_dh_bwselect'
    local bw        `_dh_bw'
    local n_total   `_dh_n'

    // Set biters to 0 when bootstrap is off.
    if `bstrap' == 0 {
        local biters = 0
    }

    // Display estimation header
    di as text ""
    di as text "{hline 72}"
    di as text "Doubly Robust Uniform Confidence Bands for CATT"
    di as text "{hline 72}"
    di as text ""
    di as text "  Outcome variable:   `depvar'"
    di as text "  Panel id:           `id'"
    di as text "  Time variable:      `time'"
    di as text "  Group variable:     `group'"
    di as text "  Covariate (z):      `z'"
    if "`xformula_display'" != "" {
        di as text "  Covariate spec:     `xformula_display'"
    }
    di as text ""
    di as text "  Observations (n):   `n_total'"
    di as text "  Polynomial order:   `porder'"
    di as text "  Kernel:             `kernel'"
    di as text "  Control group:      `control'"
    di as text "  Anticipation:       `anticip'"
    di as text "  Significance level: `alp'"
    di as text "  BW selection:       `bwselect'"
    if `bstrap' == 1 {
        di as text "  Bootstrap:          ON (`biters' iterations)"
    }
    else {
        di as text "  Bootstrap:          OFF"
    }
    if `uniform' == 1 {
        di as text "  Uniform over:       (g, t, z)"
    }
    else {
        di as text "  Uniform over:       z only"
    }
    if `pretrend' == 1 {
        di as text "  Pre-trend test:     ON"
    }
    di as text ""

    // Initialize Mata data structures
    // Preserve previous estimates for recovery on failure.
    tempname _dh_prev_est
    capture _est hold `_dh_prev_est', restore nullok

    ereturn clear
    capture noisily mata: didhetero_init_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }

    // Handle user-specified (g,t) evaluation pairs
    if "`gteval'" != "" {
        // Parse gteval numlist: format is (g1 t1 g2 t2 ...)
        local ngt_tokens : word count `gteval'
        if mod(`ngt_tokens', 2) != 0 {
            restore
            capture _est unhold `_dh_prev_est'
            di as error "gteval() must contain an even number of values (g1 t1 g2 t2 ...)"
            exit 198
        }
        local _ngt_pairs = `ngt_tokens' / 2
        local _gteval_user `gteval'

        // Validate bandwidth specification for user-defined pairs.
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

    // Run estimation
    capture noisily mata: didhetero_run_from_ado()
    local _dh_rc = _rc
    if `_dh_rc' {
        restore
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }
    mata: st_local("_dh_effective_uniformall", strofreal(didhetero_get_uniformall(), "%9.0g"))

    // Remove stale bootstrap results when bootstrap is off.
    if `bstrap' == 0 {
        capture matrix drop e(c_check)
    }

    // Post estimation results
    restore
    capture noisily _didhetero_post_eclass
    local _dh_rc = _rc
    if `_dh_rc' {
        capture _est unhold `_dh_prev_est'
        exit `_dh_rc'
    }
    capture _est unhold `_dh_prev_est', not

    // Post eclass scalars and macros.
    ereturn local cmd        "didhetero"
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
    ereturn scalar uniformall = `_dh_effective_uniformall'
    ereturn scalar pretrend  = `pretrend'

    // Display results table
    _didhetero_display

end
