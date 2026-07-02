*! catt_gt_predict.ado
*! Post-estimation predict for catt_gt
*! Outputs CATT(g,t,z) point estimates from e(results) (paper Eq.2.3, Lemma 2)
*!
*! Usage (after catt_gt):
*!   predict varname            - CATT point estimates (default)
*!   predict varname, se        - standard errors
*!   predict varname, ci1       - analytical CI (generates varname_lb, varname_ub)
*!   predict varname, ci2       - bootstrap CI  (generates varname_lb, varname_ub)
*!   predict varname, bw        - bandwidth used
*!   predict varname, zval      - z evaluation points
*!   predict varname, gval      - group values
*!   predict varname, tval      - time period values

program define catt_gt_predict
    version 16.0

    // ─── Guard: must follow catt_gt or didhetero ────────────────────────────────
    if "`e(cmd)'" != "catt_gt" & "`e(cmd)'" != "didhetero" {
        di as error "catt_gt_predict requires estimation results from {bf:catt_gt} or {bf:didhetero}"
        error 301
    }

    // ─── Parse syntax ─────────────────────────────────────────────────────────
    syntax newvarname [, SE CI1 CI2 BW Zval Gval Tval]

    local varname `varlist'

    // Count options (ci1/ci2 are exclusive; at most one option allowed)
    local n_opts = ("`se'" != "") + ("`ci1'" != "") + ("`ci2'" != "") ///
                 + ("`bw'" != "") + ("`zval'" != "") + ("`gval'" != "") ///
                 + ("`tval'" != "")
    if `n_opts' > 1 {
        di as error "only one option may be specified"
        exit 198
    }

    // ─── Extract results matrix ───────────────────────────────────────────────
    capture confirm matrix e(results)
    if _rc {
        di as error "estimation results matrix e(results) not found"
        error 301
    }

    tempname R
    matrix `R' = e(results)
    local nrows = rowsof(`R')

    // Check dataset has enough observations
    if `nrows' > _N {
        di as error "dataset has fewer observations than evaluation points"
        di as error "  evaluation points: `nrows', observations: " _N
        di as error "  consider expanding the dataset or using a frame"
        error 2001
    }

    // ─── Determine column(s) to extract ──────────────────────────────────────
    // e(results) columns: 1=g, 2=t, 3=z, 4=est, 5=se,
    //                     6=ci1_lower, 7=ci1_upper, 8=ci2_lower, 9=ci2_upper, 10=bw

    local gen_pair 0

    if "`se'" != "" {
        local col 5
        local desc "standard error"
    }
    else if "`ci1'" != "" {
        local col_lb 6
        local col_ub 7
        local gen_pair 1
        local desc "analytical confidence interval"
    }
    else if "`ci2'" != "" {
        // Verify bootstrap was run
        if e(bstrap) == 0 {
            di as error "bootstrap CI not available: estimation was run without bootstrap"
            di as error "  re-run catt_gt with bootstrap enabled to obtain ci2"
            error 198
        }
        local col_lb 8
        local col_ub 9
        local gen_pair 1
        local desc "bootstrap confidence interval"
    }
    else if "`bw'" != "" {
        local col 10
        local desc "bandwidth"
    }
    else if "`zval'" != "" {
        local col 3
        local desc "z evaluation point"
    }
    else if "`gval'" != "" {
        local col 1
        local desc "group value"
    }
    else if "`tval'" != "" {
        local col 2
        local desc "time period value"
    }
    else {
        // Default: point estimate
        local col 4
        local desc "CATT point estimate"
    }

    // ─── Generate variable(s) ─────────────────────────────────────────────────

    if `gen_pair' {
        // CI options produce two variables: varname_lb, varname_ub
        local vname_lb "`varname'_lb"
        local vname_ub "`varname'_ub"

        // Check that derived names do not already exist
        capture confirm variable `vname_lb'
        if !_rc {
            di as error "variable {bf:`vname_lb'} already exists"
            exit 110
        }
        capture confirm variable `vname_ub'
        if !_rc {
            di as error "variable {bf:`vname_ub'} already exists"
            exit 110
        }

        quietly generate double `vname_lb' = .
        quietly generate double `vname_ub' = .

        forvalues i = 1/`nrows' {
            quietly replace `vname_lb' = `R'[`i', `col_lb'] in `i'
            quietly replace `vname_ub' = `R'[`i', `col_ub'] in `i'
        }

        label variable `vname_lb' "catt_gt `desc' lower bound"
        label variable `vname_ub' "catt_gt `desc' upper bound"

        di as text ""
        di as text "Generated variables:"
        di as text "  {bf:`vname_lb'} — `desc' lower bound"
        di as text "  {bf:`vname_ub'} — `desc' upper bound"
        di as text "  (`nrows' evaluation points stored in obs 1–`nrows')"
    }
    else {
        // Single variable
        quietly generate double `varname' = .

        forvalues i = 1/`nrows' {
            quietly replace `varname' = `R'[`i', `col'] in `i'
        }

        label variable `varname' "catt_gt `desc'"

        di as text ""
        di as text "Generated variable:"
        di as text "  {bf:`varname'} — `desc'"
        di as text "  (`nrows' evaluation points stored in obs 1–`nrows')"
    }
end
