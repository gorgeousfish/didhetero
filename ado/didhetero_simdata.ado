*! didhetero_simdata.ado
*! Version 0.1.0
*!
*! Generate simulation data following the DGP in Section 6 of
*! Imai, Qin, and Yanagi (2025).
*!
*! Syntax:
*!   didhetero_simdata, n(integer) tau(integer) [seed(integer) hc dimx(integer 1)
*!                      dgpy(integer 1) discrete clear]
*!
*! Required:
*!   n(integer)     - number of cross-sectional units
*!   tau(integer)   - number of time periods (must be > 2)
*!
*! Optional:
*!   seed(integer)  - random seed for reproducibility
*!   hc             - use heteroscedastic error terms (default: homoscedastic)
*!   dimx(integer)  - dimension of covariates X (default: 1)
*!   dgpy(integer)  - DGP for treated potential outcome:
*!                    1 = nonlinear sin(pi * Z) term, 2 = linear Z * G / t term
*!   discrete       - use discrete Z in {-1, 0, 1} (default: continuous Z ~ N(0,1))
*!   clear          - allow replacing data in memory
*!
*! Output variables:
*!   id     - Individual ID (long)
*!   period - Time period (long)
*!   Y      - Observed outcome (double)
*!   G      - Treatment group: 0=never-treated, 2..tau (long)
*!   Z      - Covariate (double)
*!   X1..   - Additional covariates when dimx > 1 (double)

program define didhetero_simdata, rclass
    version 16.0

    // =========================================================================
    // Step 1: Parse syntax
    // =========================================================================
    syntax, n(integer) tau(integer)          ///
        [seed(integer -1) hc dimx(integer 1) ///
         dgpy(integer 1) discrete clear]

    // =========================================================================
    // Step 2: Validate inputs
    // =========================================================================
    if `n' <= 0 {
        display as error "n() must be a positive integer"
        exit 198
    }
    if `tau' <= 2 {
        display as error "tau() must be greater than 2"
        exit 198
    }
    if `dgpy' != 1 & `dgpy' != 2 {
        display as error "dgpy() must be 1 or 2"
        exit 198
    }
    if `dimx' <= 0 {
        display as error "dimx() must be a positive integer"
        exit 198
    }
    if `seed' < -1 {
        display as error "seed() must be -1 or a nonnegative integer"
        exit 198
    }

    // =========================================================================
    // Step 3: Check data in memory
    // =========================================================================
    if "`clear'" == "" {
        // Protect any existing dataset structure, including zero-observation
        // datasets that still carry variables and would be overwritten.
        if c(k) > 0 | c(N) > 0 {
            display as error "data in memory would be lost; specify clear option"
            exit 4
        }
    }

    // =========================================================================
    // Step 4: Map options to numeric flags
    // =========================================================================
    local hc_flag = ("`hc'" != "")
    local continuous_flag = ("`discrete'" == "")

    // =========================================================================
    // Step 5: Ensure Mata backend is available
    // =========================================================================
    // Use the shared backend loader first (handles mlib index + source fallback)
    capture noisily _dh_ensure_backend
    if _rc {
        // Fall back to the simdata-specific loader
        quietly _didhetero_simdata_ensure_loaded
    }
    else {
        // Verify simdata function is available after shared loader
        capture quietly mata: mata describe _didhetero_simdata_to_stata()
        if _rc {
            quietly _didhetero_simdata_ensure_loaded
        }
    }

    // =========================================================================
    // Step 6: Set seed if specified
    // =========================================================================
    local seed_explicit = (`seed' != -1)
    if `seed_explicit' {
        set seed `seed'
    }

    // =========================================================================
    // Step 7: Call Mata core and write data to Stata
    // =========================================================================
    local nobs = `n' * `tau'
    tempfile simdata_tmp

    // Preserve the caller's dataset until replacement data are fully generated.
    preserve
    clear

    capture noisily mata: _didhetero_simdata_to_stata(`n', `tau', `hc_flag', `dimx', `dgpy', `continuous_flag')
    local bridge_rc = _rc
    if `bridge_rc' {
        restore
        exit `bridge_rc'
    }

    capture noisily save "`simdata_tmp'", replace
    local save_rc = _rc
    restore
    if `save_rc' {
        exit `save_rc'
    }

    use "`simdata_tmp'", clear

    // =========================================================================
    // Step 8: Sort and label
    // =========================================================================
    sort id period

    label variable id     "Individual ID"
    label variable period "Time period"
    label variable Y      "Observed outcome"
    label variable G      "Treatment group (0=never-treated)"
    label variable Z      "Covariate"

    if `dimx' > 1 {
        forvalues k = 1/`=`dimx'-1' {
            capture label variable X`k' "Additional covariate `k'"
        }
    }

    // =========================================================================
    // Step 9: Return values
    // =========================================================================
    return scalar n = `n'
    return scalar tau = `tau'
    return scalar dimx = `dimx'
    return scalar dgpy = `dgpy'
    return scalar hc = `hc_flag'
    return scalar continuous = `continuous_flag'
    return scalar seed_explicit = `seed_explicit'
    if `seed_explicit' {
        return scalar seed = `seed'
    }
    else {
        return scalar seed = .
    }

    quietly tab G
    local ngroups = r(r)
    return scalar n_groups = `ngroups'

    // Build group list string
    quietly levelsof G, local(glevels)
    return local groups "`glevels'"

    quietly count if G == 0
    local n_never = r(N) / `tau'
    local n_treated = `n' - `n_never'
    local pct_treated : display %5.1f 100 * `n_treated' / `n'
    local groups_note " (including never-treated)"
    if `n_never' == 0 {
        local groups_note " (no never-treated realized in sample)"
    }

    // =========================================================================
    // Step 10: Display summary
    // =========================================================================

    display as text ""
    display as text "{hline 50}"
    display as text "  didhetero_simdata: simulation data generated"
    display as text "{hline 50}"
    display as text "  Units (n):          " as result `n'
    display as text "  Periods (tau):      " as result `tau'
    display as text "  Observations:       " as result `nobs'
    display as text "  Groups:             " as result `ngroups' as text "`groups_note'"
    display as text "  Treated share:      " as result "`pct_treated'%"
    display as text "  DGP-Y:             " as result `dgpy'
    display as text "  Covariates (dimx):  " as result `dimx'
    if "`hc'" != "" {
        display as text "  Errors:             " as result "heteroscedastic"
    }
    else {
        display as text "  Errors:             " as result "homoscedastic"
    }
    if "`discrete'" != "" {
        display as text "  Z type:             " as result "discrete {-1, 0, 1}"
    }
    else {
        display as text "  Z type:             " as result "continuous N(0,1)"
    }
    display as text "{hline 50}"
    display as text ""
end

program define _didhetero_simdata_ensure_loaded
    version 16.0

    // Prefer already-loaded Mata code.
    capture quietly mata: mata which _didhetero_simdata_to_stata()
    if !_rc {
        exit
    }

    // If a packaged mlib is available on adopath, indexing should expose it.
    capture quietly mata: mata mlib index
    capture quietly mata: mata which _didhetero_simdata_to_stata()
    if !_rc {
        exit
    }

    // Fall back to the source .mata file so ado-only users do not need a
    // prebuilt library artifact or a manual indexing step.
    quietly findfile didhetero_simdata.ado
    local ado_file `"`r(fn)'"'
    local mata_file = subinstr(`"`ado_file'"', "didhetero_simdata.ado", "../mata/didhetero_simdata.mata", 1)

    capture confirm file `"`mata_file'"'
    if _rc {
        display as error "didhetero_simdata Mata backend not found"
        exit 601
    }

    capture noisily do `"`mata_file'"'
    if _rc {
        display as error "didhetero_simdata Mata backend failed to load"
        exit _rc
    }

    capture quietly mata: mata which _didhetero_simdata_to_stata()
    if _rc {
        display as error "didhetero_simdata Mata backend remains unavailable"
        exit 3499
    }
end
