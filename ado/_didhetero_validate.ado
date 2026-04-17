*! _didhetero_validate.ado
*! Input validation for catt_gt command
*! Part of didhetero-stata package
*!
*! Validates all user inputs before entering Mata computation layer.
*! Called by catt_gt.ado after syntax parsing.
*!
*! Validates:
*!   1. Parameter ranges (porder, kernel, control_group, anticipation, alp, biters)
*!   2. Panel structure (unique id+time, balanced panel, time-invariance of G and Z)
*!   3. Variable types (numeric depvar/z, integer group/time)
*!   4. Missing values (detect and exclude entire id)
*!   5. Sample size warnings (n < 100, n_g < 30)
*!
*! On error: exits with appropriate error code and message
*! On success: leaves validated, sorted data in memory

program define _didhetero_validate
    version 16.0

    // =========================================================================
    // Syntax parsing
    // =========================================================================
    syntax varlist(min=1 max=1 numeric) [if] [in], ///
        Id(varname)                                  ///
        Time(varname)                                ///
        Group(varname)                               ///
        Z(varname)                                   ///
        Zeval(numlist)                               ///
        [GTeval(numlist)]                            ///
        [Xformula(string)]                           ///
        [Porder(integer 2)]                          ///
        [Kernel(string)]                             ///
        [Control_group(string)]                      ///
        [Anticipation(integer 0)]                    ///
        [Alp(real 0.05)]                             ///
        [Biters(integer 1000)]                       ///
        [BSTrap]                                     ///
        [UNIFormall]                                 ///
        [PREtrend]                                   ///
        [BWselect(string)]                           ///
        [BW(numlist)]

    local depvar `varlist'

    // Set string defaults (Stata syntax does not support string defaults)
    if "`kernel'" == "" local kernel "gau"
    if "`control_group'" == "" local control_group "notyettreated"
    if "`bwselect'" == "" local bwselect "IMSE1"

    local xformula_display `"`xformula'"'
    local xformula_has_intercept 1
    if `"`xformula'"' != "" {
        mata: didhetero_parse_xformula_locals()
        local xformula_display `"`dh_xformula_display'"'
        local xformula_has_intercept `dh_xformula_has_intercept'
        local xformula

        forvalues i = 1/`dh_xformula_n' {
            local kind_name dh_xformula_kind`i'
            local payload_name dh_xformula_payload`i'
            local term_name dh_xformula_term`i'
            local kind ``kind_name''
            local payload ``payload_name''
            local term ``term_name''

            if "`kind'" == "var" {
                local xformula `xformula' `payload'
            }
            else if "`kind'" == "expr" {
                local gen_ix = `i'
                local genvar "__dhxf_`gen_ix'"
                capture confirm new variable `genvar'
                while _rc {
                    local ++gen_ix
                    local genvar "__dhxf_`gen_ix'"
                    capture confirm new variable `genvar'
                }

                capture gen double `genvar' = `payload'
                if _rc {
                    di as error "invalid xformula() term: `term'"
                    exit 198
                }

                local xformula `xformula' `genvar'
            }
        }

        local xformula : list uniq xformula

        local _dh_has_z_in_xformula : list posof "`z'" in xformula
        if !`_dh_has_z_in_xformula' {
            di as error ///
                "xformula() must explicitly include the z() variable `z'; for the default specification, use xformula(`z')"
            exit 198
        }
    }

    // Parse boolean options
    local bstrap_flag = ("`bstrap'" != "")
    local uniformall_flag = ("`uniformall'" != "")
    local pretrend_flag = ("`pretrend'" != "")

    // pretrend only affects the automatically constructed (g,t) domain.
    // Explicit gteval() already fixes that domain, so the combination is invalid.
    if "`gteval'" != "" & `pretrend_flag' {
        di as error "pretrend cannot be combined with gteval(); pretrend is only applicable when gteval() is omitted"
        exit 198
    }

    // =========================================================================
    // Step 1: Parameter range validation (FR-01)
    // cf. Imai, Qin, and Yanagi (2025), parameter constraints
    // =========================================================================

    // FR-01-01: porder must be 1 or 2
    if !inlist(`porder', 1, 2) {
        di as error "porder must be 1 or 2"
        exit 198
    }

    // FR-01-02: kernel must be "gau" or "epa"
    if !inlist("`kernel'", "gau", "epa") {
        di as error "kernel must be 'gau' or 'epa'"
        exit 198
    }

    // FR-01-03: control_group must be "nevertreated" or "notyettreated"
    if !inlist("`control_group'", "nevertreated", "notyettreated") {
        di as error "control_group must be 'nevertreated' or 'notyettreated'"
        exit 198
    }

    // FR-01-04: anticipation must be a non-negative integer
    // Stata special: check for missing value (. >= 0 is true in Stata)
    if `anticipation' >= . {
        di as error "anticipation must be a non-negative integer"
        exit 198
    }
    if `anticipation' < 0 {
        di as error "anticipation must be a non-negative integer"
        exit 198
    }
    if mod(`anticipation', 1) != 0 {
        di as error "anticipation must be a non-negative integer"
        exit 198
    }

    // FR-01-05: alp must be strictly between 0 and 1
    // Stata special: check for missing value (. > 0 is true in Stata)
    if `alp' >= . {
        di as error "alp must be strictly between 0 and 1"
        exit 198
    }
    if `alp' <= 0 | `alp' >= 1 {
        di as error "alp must be strictly between 0 and 1"
        exit 198
    }

    // FR-01-06: biters must be >= 1 only when bootstrap is enabled.
    if `bstrap_flag' & `biters' < 1 {
        di as error "When bstrap = TRUE, biters must be a positive number"
        exit 198
    }

    // FR-01-08: bwselect validation
    if !inlist("`bwselect'", "IMSE1", "IMSE2", "US1", "manual") {
        di as error "bwselect must be 'IMSE1', 'IMSE2', 'US1', or 'manual'"
        exit 198
    }

    // FR-01-09: bw and bwselect consistency
    if "`bwselect'" == "manual" & "`bw'" == "" {
        di as error "bw must be specified manually when bwselect = 'manual'"
        exit 198
    }
    if "`bwselect'" != "manual" & "`bw'" != "" {
        di as error "bw can be specified manually only when bwselect = 'manual'"
        exit 198
    }

    // Automatic bandwidth selection integrates over the zeval grid and is
    // undefined for a single-point evaluation set. Single-point zeval remains
    // valid with an explicitly supplied manual bandwidth.
    local num_zeval : word count `zeval'
    if "`bwselect'" != "manual" & `num_zeval' < 2 {
        di as error "automatic bandwidth selection requires at least two zeval points; use bwselect(manual) with bw() for a single-point zeval()"
        exit 198
    }

    // FR-01-10: bw values must be positive (if specified)
    if "`bw'" != "" {
        foreach val of numlist `bw' {
            if `val' <= 0 | `val' >= . {
                di as error "bw must contain only positive values"
                exit 198
            }
        }

        local n_bw_tokens : word count `bw'
        if `n_bw_tokens' > 1 {
            if "`gteval'" == "" {
                di as error ///
                    "bw must be a positive scalar or vector whose length equals to the number of gteval"
                exit 198
            }

            local n_gteval_tokens : word count `gteval'
            if mod(`n_gteval_tokens', 2) != 0 {
                di as error ///
                    "gteval() must contain an even number of values (g1 t1 g2 t2 ...)"
                exit 198
            }

            local n_gteval_pairs = `n_gteval_tokens' / 2
            if `n_bw_tokens' != `n_gteval_pairs' {
                di as error ///
                    "bw must be a positive scalar or vector whose length equals to the number of gteval"
                exit 198
            }
        }
    }

    // =========================================================================
    // Step 2: Variable type validation (FR-03)
    // =========================================================================

    // FR-03-01: depvar must be numeric
    confirm numeric variable `depvar'

    // FR-03-03: z must be numeric
    confirm numeric variable `z'

    // FR-03-05: group must be numeric
    confirm numeric variable `group'

    // FR-03-06: time must be numeric
    confirm numeric variable `time'

    // FR-03-02: Promote depvar to double if needed
    local depvar_type : type `depvar'
    if "`depvar_type'" != "double" {
        di as text "Note: `depvar' recast to double for numerical precision"
        quietly recast double `depvar'
    }

    // FR-03-04: Promote z to double if needed
    local z_type : type `z'
    if "`z_type'" != "double" {
        di as text "Note: `z' recast to double for numerical precision"
        quietly recast double `z'
    }

    // FR-03-05: Check group contains integer values
    tempvar g_intcheck
    quietly gen byte `g_intcheck' = (`group' == floor(`group')) if !missing(`group')
    quietly summarize `g_intcheck', meanonly
    if r(min) == 0 {
        di as error "group variable must contain integer values"
        exit 198
    }

    // FR-03-06: Check time contains integer values
    tempvar t_intcheck
    quietly gen byte `t_intcheck' = (`time' == floor(`time')) if !missing(`time')
    quietly summarize `t_intcheck', meanonly
    if r(min) == 0 {
        di as error "time variable must contain integer values"
        exit 198
    }

    // =========================================================================
    // Step 3: Missing value detection and exclusion (FR-03-07/08)
    // =========================================================================

    // Detect missing values in key variables
    tempvar has_missing
    quietly gen byte `has_missing' = ///
        missing(`depvar') | missing(`z') | ///
        missing(`group') | missing(`time') | missing(`id')

    // Propagate to entire id (exclude whole individual if any obs missing)
    tempvar id_has_missing
    quietly bysort `id' (`time'): egen byte `id_has_missing' = max(`has_missing')

    quietly count if `id_has_missing' == 1
    local nmiss = r(N)
    if `nmiss' > 0 {
        di as text "Note: `nmiss' observations dropped due to missing values"
        quietly drop if `id_has_missing' == 1
        quietly count
        if r(N) == 0 {
            error 2000
        }
    }

    // Also handle xformula variables if specified
    if "`xformula'" != "" {
        foreach xvar of local xformula {
            capture confirm numeric variable `xvar'
            if _rc {
                di as error "xformula variable `xvar' must be numeric"
                exit 198
            }
            tempvar xmiss
            quietly gen byte `xmiss' = missing(`xvar')
            tempvar xid_miss
            quietly bysort `id' (`time'): egen byte `xid_miss' = max(`xmiss')
            quietly count if `xid_miss' == 1
            local xnmiss = r(N)
            if `xnmiss' > 0 {
                di as text "Note: `xnmiss' observations dropped due to missing values in `xvar'"
                quietly drop if `xid_miss' == 1
                quietly count
                if r(N) == 0 {
                    error 2000
                }
            }
        }
    }

    // =========================================================================
    // Step 4: Panel structure validation (FR-02)
    // =========================================================================

    // FR-02-01: id + time must uniquely identify observations
    tempvar dup_count
    quietly duplicates tag `id' `time', gen(`dup_count')
    quietly summarize `dup_count', meanonly
    if r(max) > 0 {
        di as error "id and time do not uniquely identify observations"
        exit 198
    }

    // FR-02-02: Balanced panel check
    tempvar T_i
    quietly bysort `id': gen long `T_i' = _N
    quietly summarize `T_i', meanonly
    if r(min) != r(max) {
        di as error ///
            "Panel is not balanced: each id must have the same number of time periods"
        exit 198
    }

    // FR-02-02b: Every id must share the same sorted time support
    // A balanced panel requires identical time values for each within-id
    // position after sorting by time, not just equal observation counts.
    tempvar t_seq t_min t_max
    quietly bysort `id' (`time'): gen long `t_seq' = _n
    quietly bysort `t_seq': egen double `t_min' = min(`time')
    quietly bysort `t_seq': egen double `t_max' = max(`time')
    quietly count if `t_min' != `t_max'
    if r(N) > 0 {
        di as error ///
            "Panel is not balanced: each id must share the same set of time periods"
        exit 198
    }

    // FR-02-03: group must be time-invariant within id
    tempvar g_sd
    quietly bysort `id': egen double `g_sd' = sd(`group')
    quietly summarize `g_sd', meanonly
    // Use tolerance for floating-point sd() artifacts (sd can be ~1e-15 for
    // truly constant values due to double-precision arithmetic)
    if r(max) > 1e-8 {
        di as error "group variable is not time-invariant within id"
        exit 198
    }

    // FR-02-04: z must be time-invariant within id
    tempvar z_sd
    quietly bysort `id': egen double `z_sd' = sd(`z')
    quietly summarize `z_sd', meanonly
    if r(max) > 1e-8 {
        di as error "z variable is not time-invariant within id"
        exit 198
    }

    // FR-02-04x: xformula variables must also be unit-level pre-treatment
    // covariates. The paper defines X_i, not X_it, so every resolved
    // xformula() column must be time-invariant within id.
    if "`xformula'" != "" {
        foreach xvar of local xformula {
            tempvar x_sd
            quietly bysort `id': egen double `x_sd' = sd(`xvar')
            quietly summarize `x_sd', meanonly
            if r(max) > 1e-8 {
                di as error "xformula variable `xvar' is not time-invariant within id"
                exit 198
            }
        }
    }

    // FR-02-04a: the continuous-Z estimator requires enough distinct support
    // points across ids to identify the cubic local polynomial used for the
    // Stage-1 density derivative estimate. Fewer than four distinct values
    // globally cannot satisfy this necessary rank condition.
    tempvar z_id_tag z_support_tag
    quietly egen byte `z_id_tag' = tag(`id')
    quietly egen byte `z_support_tag' = tag(`z') if `z_id_tag' == 1
    quietly count if `z_support_tag' == 1
    local n_z_support = r(N)
    if `n_z_support' < 4 {
        di as error ///
            "z() must have at least four distinct support points across ids; the continuous-Z estimator is not identified with fewer support points"
        exit 198
    }

    // FR-02-04b: treated cohorts must start strictly after the first period
    // Under staggered adoption, D_1 = 0 almost surely. The package encodes
    // never-treated units as group==0, so every treated cohort must satisfy
    // group > first observed period.
    quietly summarize `time', meanonly
    local first_time = r(min)
    quietly count if `group' != 0 & `group' <= `first_time'
    if r(N) > 0 {
        di as error ///
            "group() must be 0 for never-treated units or a treatment time strictly after the first observed period"
        exit 198
    }

    // FR-02-05: nevertreated control group requires at least one group==0 unit
    if "`control_group'" == "nevertreated" {
        quietly count if `group' == 0
        if r(N) == 0 {
            di as error ///
                "control_group(nevertreated) requires at least one never-treated unit with group==0"
            exit 198
        }
    }

    // =========================================================================
    // Step 5: Sample size warnings (FR-04)
    // =========================================================================

    // FR-04-01: Total sample size warning
    tempvar id_tag
    quietly egen byte `id_tag' = tag(`id')
    quietly count if `id_tag' == 1
    local n_total = r(N)
    if `n_total' < 100 {
        di as text ///
            "Warning: Sample size (n=`n_total') may be too small for reliable nonparametric estimation"
    }

    // FR-04-02/03: Per-group sample size warning
    quietly levelsof `group', local(group_levels)
    foreach g of local group_levels {
        quietly count if `id_tag' == 1 & `group' == `g'
        local n_g = r(N)
        if `n_g' < 30 {
            if `g' == 0 {
                di as text "Warning: Never-treated group has only `n_g' observations"
            }
            else {
                di as text "Warning: Treatment group g=`g' has only `n_g' individuals (< 30)"
            }
        }
    }

    // =========================================================================
    // Step 6: Sort data (FR-05-01)
    // Sort by id and time for consistent panel layout
    // =========================================================================
    sort `id' `time'

    // =========================================================================
    // Return validated parameters as locals for caller
    // =========================================================================
    c_local _dh_depvar    `depvar'
    c_local _dh_id        `id'
    c_local _dh_time      `time'
    c_local _dh_group     `group'
    c_local _dh_z         `z'
    c_local _dh_zeval     `zeval'
    c_local _dh_xformula  `xformula'
    c_local _dh_xformula_display `xformula_display'
    c_local _dh_xformula_has_intercept `xformula_has_intercept'
    c_local _dh_porder    `porder'
    c_local _dh_kernel    `kernel'
    c_local _dh_control   `control_group'
    c_local _dh_anticip   `anticipation'
    c_local _dh_alp       `alp'
    c_local _dh_biters    `biters'
    c_local _dh_bstrap    `bstrap_flag'
    c_local _dh_uniform   `uniformall_flag'
    c_local _dh_pretrend  `pretrend_flag'
    c_local _dh_bwselect  `bwselect'
    c_local _dh_bw        `bw'
    c_local _dh_n         `n_total'

end
