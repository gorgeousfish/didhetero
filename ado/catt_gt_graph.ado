*! catt_gt_graph.ado
*! Version 0.1.0
*! Visualize CATT or aggregated parameter estimates
*! Post-estimation command for didhetero / catt_gt / aggte_gt
*!
*! Implements Imai, Qin, and Yanagi (2025).
*!
*! Syntax:
*!   catt_gt_graph, [plot_type(string) save_path(string) graph_opt(string)]
*!
*! Options:
*!   plot_type  - "CATT" or "Aggregated" (auto-detected if omitted)
*!   save_path  - file path to save the last graph (.png/.pdf via export, .gph via save)
*!   graph_opt  - additional Stata graph options passed to twoway
*!
*! e() matrix sources:
*!   catt_gt  stores the 10-column CATT matrix in both e(results) and e(Estimate)
*!            columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
*!   aggte_gt stores 9-column matrix as e(Estimate)
*!            columns: eval, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw

program define catt_gt_graph
    version 16.0

    // =========================================================================
    // Step 1: Parse syntax
    // =========================================================================
    syntax, [plot_type(string) save_path(string) graph_opt(string)]

    if `"`plot_type'"' != "" {
        if !inlist(`"`plot_type'"', "CATT", "Aggregated") {
            di as error "catt_gt_graph: plot_type() must be " `"""' "CATT" `"""' " or " `"""' "Aggregated" `"""'
            di as error "  received: `plot_type'"
            exit 198
        }
    }

    // =========================================================================
    // Step 2: Locate e() results matrix
    // =========================================================================
    // Auto mode preserves the historical precedence (prefer e(Estimate)).
    // Explicit plot_type() requests instead select the matching result object.

    tempname _Est
    local _source_matrix ""
    local has_estimate = 0
    local has_results = 0
    local mode ""

    capture confirm matrix e(Estimate)
    if !_rc {
        local has_estimate = 1
    }
    capture confirm matrix e(results)
    if !_rc {
        local has_results = 1
    }

    if !`has_estimate' & !`has_results' {
        di as error "catt_gt_graph: no estimation results found."
        di as error "  Neither e(Estimate) nor e(results) exists."
        di as error "  Run didhetero, catt_gt, or aggte_gt first."
        exit 301
    }

    if `"`plot_type'"' == "CATT" {
        if `has_results' {
            matrix `_Est' = e(results)
            local _source_matrix "results"
        }
        else {
            matrix `_Est' = e(Estimate)
            local _source_matrix "Estimate"
        }

        local ncols = colsof(`_Est')
        if `ncols' != 10 & `has_estimate' & "`_source_matrix'" != "Estimate" {
            matrix `_Est' = e(Estimate)
            local _source_matrix "Estimate"
            local ncols = colsof(`_Est')
        }
        if `ncols' != 10 {
            di as error "catt_gt_graph: plot_type(CATT) requires a 10-column CATT matrix."
            di as error "  e(results) and e(Estimate) do not currently store CATT results."
            exit 198
        }
        local mode "CATT"
    }
    else if `"`plot_type'"' == "Aggregated" {
        if `has_estimate' {
            matrix `_Est' = e(Estimate)
            local _source_matrix "Estimate"
        }
        else {
            matrix `_Est' = e(results)
            local _source_matrix "results"
        }

        local ncols = colsof(`_Est')
        if `ncols' != 9 & `has_results' & "`_source_matrix'" != "results" {
            matrix `_Est' = e(results)
            local _source_matrix "results"
            local ncols = colsof(`_Est')
        }
        if `ncols' != 9 {
            di as error "catt_gt_graph: plot_type(Aggregated) requires a 9-column aggregated matrix."
            di as error "  e(Estimate) and e(results) do not currently store aggregated results."
            exit 198
        }
        local mode "Aggregated"
    }
    else {
        if `has_estimate' {
            matrix `_Est' = e(Estimate)
            local _source_matrix "Estimate"
        }
        else {
            matrix `_Est' = e(results)
            local _source_matrix "results"
        }

        local ncols = colsof(`_Est')

        // =========================================================================
        // Step 3: Determine column count for mode detection
        // =========================================================================
        //   10 columns => CATT     (from catt_gt / e(results))
        //    9 columns => Aggregated (from aggte_gt / e(Estimate))
        if `ncols' == 10 {
            local mode "CATT"
        }
        else if `ncols' == 9 {
            local mode "Aggregated"
        }
        else {
            di as error "catt_gt_graph: e(`_source_matrix') has `ncols' columns;"
            di as error "  expected 10 (CATT) or 9 (Aggregated)."
            exit 198
        }
    }

    di as text ""
    di as text "catt_gt_graph: mode = `mode' (source = e(`_source_matrix'), columns = `ncols')"
    di as text ""

    // =========================================================================
    // Step 6: Confidence band selection (panel-aware)
    // =========================================================================
    // Determine CI column indices once, then decide the band source panel by
    // panel. A global "any ci2 exists" rule is incorrect for aggregated
    // outputs because some eval-specific panels may still only have ci1.
    local nrows = rowsof(`_Est')
    if "`mode'" == "CATT" {
        local z_col = 3
        local est_col = 4
        local ci1_lower_col = 6
        local ci1_upper_col = 7
        local ci2_lower_col = 8
        local ci2_upper_col = 9
    }
    else {
        local z_col = 2
        local est_col = 3
        local ci1_lower_col = 5
        local ci1_upper_col = 6
        local ci2_lower_col = 7
        local ci2_upper_col = 8
    }

    // Count non-missing bootstrap rows globally for summary display only.
    local ci2_nonmiss = 0
    forvalues i = 1/`nrows' {
        if `_Est'[`i', `ci2_lower_col'] < . & `_Est'[`i', `ci2_upper_col'] < . {
            local ci2_nonmiss = `ci2_nonmiss' + 1
        }
    }

    if `ci2_nonmiss' == 0 {
        di as text "  Confidence bands: analytical (ci1)"
    }
    else if `ci2_nonmiss' == `nrows' {
        di as text "  Confidence bands: bootstrap (ci2)"
    }
    else {
        di as text "  Confidence bands: panel-specific auto (prefer ci2; fallback to ci1)"
    }

    // =========================================================================
    // Step 7: Pre-trends zero reference line
    // =========================================================================
    local pretrend_opt ""
    capture scalar _pt = e(pretrend)
    if !_rc {
        if e(pretrend) == 1 {
            local pretrend_opt `"yline(0, lcolor(red) lpattern(dash))"'
        }
    }

    // =========================================================================
    // Step 8: Plotting loop
    // =========================================================================
    if "`mode'" == "CATT" {
        // -----------------------------------------------------------------
        // CATT mode: one plot per (g, t) pair from e(gteval)
        // -----------------------------------------------------------------
        tempname _gteval _subset _row
        capture confirm matrix e(gteval)
        if _rc {
            di as error "catt_gt_graph: e(gteval) not found. Required for CATT mode."
            exit 301
        }
        matrix `_gteval' = e(gteval)
        local n_gt = rowsof(`_gteval')
        local graph_count = 0

        forvalues k = 1/`n_gt' {
            local g1 = `_gteval'[`k', 1]
            local t1 = `_gteval'[`k', 2]

            // Count matching rows
            local match_count = 0
            forvalues i = 1/`nrows' {
                if `_Est'[`i', 1] == `g1' & `_Est'[`i', 2] == `t1' {
                    local match_count = `match_count' + 1
                }
            }

            if `match_count' == 0 {
                di as text "  Warning: no rows for g=`g1', t=`t1'. Skipping."
                continue
            }

            // Prefer bootstrap only when the current panel has complete ci2.
            local panel_ci2_nonmiss = 0
            forvalues i = 1/`nrows' {
                if `_Est'[`i', 1] == `g1' & `_Est'[`i', 2] == `t1' {
                    if `_Est'[`i', `ci2_lower_col'] < . & `_Est'[`i', `ci2_upper_col'] < . {
                        local panel_ci2_nonmiss = `panel_ci2_nonmiss' + 1
                    }
                }
            }

            local panel_ci_lower_col = `ci1_lower_col'
            local panel_ci_upper_col = `ci1_upper_col'
            local panel_ci_label "analytical (ci1 fallback)"
            if `panel_ci2_nonmiss' == `match_count' {
                local panel_ci_lower_col = `ci2_lower_col'
                local panel_ci_upper_col = `ci2_upper_col'
                local panel_ci_label "bootstrap (ci2)"
            }

            di as text "  Group `g1'. Time `t1': confidence bands = `panel_ci_label'"

            // Build subset matrix (match_count × 4): z, est, ci_lower, ci_upper
            matrix `_subset' = J(`match_count', 4, .)
            local r = 0
            forvalues i = 1/`nrows' {
                if `_Est'[`i', 1] == `g1' & `_Est'[`i', 2] == `t1' {
                    local r = `r' + 1
                    matrix `_subset'[`r', 1] = `_Est'[`i', `z_col']
                    matrix `_subset'[`r', 2] = `_Est'[`i', `est_col']
                    matrix `_subset'[`r', 3] = `_Est'[`i', `panel_ci_lower_col']
                    matrix `_subset'[`r', 4] = `_Est'[`i', `panel_ci_upper_col']
                }
            }

            // Convert to dataset and plot
            preserve
            clear
            quietly svmat `_subset', names(col)
            rename c1 z
            rename c2 est
            rename c3 ci_lower
            rename c4 ci_upper

            // Build graph name with Stata-safe tokens for decimal and negative values
            local g1_clean : subinstr local g1 "." "_", all
            local g1_clean : subinstr local g1_clean "-" "m", all
            local t1_clean : subinstr local t1 "." "_", all
            local t1_clean : subinstr local t1_clean "-" "m", all
            local gname "g`g1_clean'_t`t1_clean'"

            capture noisily twoway                                           ///
                (rarea ci_lower ci_upper z, color(gs12%60) lwidth(none))     ///
                (line est z, lcolor(black) lwidth(vthick)),                  ///
                title("Group `g1'. Time `t1'.")                              ///
                xtitle("z") ytitle("")                                       ///
                plotregion(fcolor(white))                                    ///
                legend(off)                                                  ///
                name(`gname', replace)                                       ///
                `pretrend_opt'                                               ///
                `graph_opt'
            local graph_rc = _rc

            restore
            if `graph_rc' {
                exit `graph_rc'
            }
            local graph_count = `graph_count' + 1
        }

        di as text ""
        di as text "catt_gt_graph: `graph_count' CATT graph(s) generated."
    }
    else {
        // -----------------------------------------------------------------
        // Aggregated mode: one plot per eval point (or single for simple)
        // -----------------------------------------------------------------
        tempname _subset

        // Get aggregation type
        local agg_type = e(type)
        local graph_count = 0

        if "`agg_type'" == "simple" {
            local panel_ci2_nonmiss = 0
            forvalues i = 1/`nrows' {
                if `_Est'[`i', `ci2_lower_col'] < . & `_Est'[`i', `ci2_upper_col'] < . {
                    local panel_ci2_nonmiss = `panel_ci2_nonmiss' + 1
                }
            }
            local panel_ci_lower_col = `ci1_lower_col'
            local panel_ci_upper_col = `ci1_upper_col'
            local panel_ci_label "analytical (ci1 fallback)"
            if `panel_ci2_nonmiss' == `nrows' {
                local panel_ci_lower_col = `ci2_lower_col'
                local panel_ci_upper_col = `ci2_upper_col'
                local panel_ci_label "bootstrap (ci2)"
            }
            di as text "  Simple weighted CATT: confidence bands = `panel_ci_label'"

            // Simple: single plot using all rows
            local subset_rows = `nrows'
            matrix `_subset' = J(`subset_rows', 4, .)
            forvalues i = 1/`nrows' {
                matrix `_subset'[`i', 1] = `_Est'[`i', `z_col']
                matrix `_subset'[`i', 2] = `_Est'[`i', `est_col']
                matrix `_subset'[`i', 3] = `_Est'[`i', `panel_ci_lower_col']
                matrix `_subset'[`i', 4] = `_Est'[`i', `panel_ci_upper_col']
            }

            preserve
            clear
            quietly svmat `_subset', names(col)
            rename c1 z
            rename c2 est
            rename c3 ci_lower
            rename c4 ci_upper

            capture noisily twoway                                           ///
                (rarea ci_lower ci_upper z, color(gs12%60) lwidth(none))     ///
                (line est z, lcolor(black) lwidth(vthick)),                  ///
                title("Simple weighted CATT.")                               ///
                xtitle("z") ytitle("")                                       ///
                plotregion(fcolor(white))                                    ///
                legend(off)                                                  ///
                name(simple, replace)                                        ///
                `pretrend_opt'                                               ///
                `graph_opt'
            local graph_rc = _rc

            restore
            if `graph_rc' {
                exit `graph_rc'
            }
            local graph_count = 1
        }
        else {
            // dynamic/group/calendar: one plot per unique eval point
            // Get unique eval values from column 1
            local eval_list ""
            forvalues i = 1/`nrows' {
                local val = `_Est'[`i', 1]
                local is_new = 1
                foreach existing of local eval_list {
                    if `val' == `existing' {
                        local is_new = 0
                    }
                }
                if `is_new' {
                    local eval_list `eval_list' `val'
                }
            }

            foreach e1 of local eval_list {
                // Count matching rows
                local match_count = 0
                forvalues i = 1/`nrows' {
                    if `_Est'[`i', 1] == `e1' {
                        local match_count = `match_count' + 1
                    }
                }

                if `match_count' == 0 {
                    continue
                }

                local panel_ci2_nonmiss = 0
                forvalues i = 1/`nrows' {
                    if `_Est'[`i', 1] == `e1' {
                        if `_Est'[`i', `ci2_lower_col'] < . & `_Est'[`i', `ci2_upper_col'] < . {
                            local panel_ci2_nonmiss = `panel_ci2_nonmiss' + 1
                        }
                    }
                }

                local panel_ci_lower_col = `ci1_lower_col'
                local panel_ci_upper_col = `ci1_upper_col'
                local panel_ci_label "analytical (ci1 fallback)"
                if `panel_ci2_nonmiss' == `match_count' {
                    local panel_ci_lower_col = `ci2_lower_col'
                    local panel_ci_upper_col = `ci2_upper_col'
                    local panel_ci_label "bootstrap (ci2)"
                }

                di as text "  Evaluation point `e1': confidence bands = `panel_ci_label'"

                // Build subset matrix
                matrix `_subset' = J(`match_count', 4, .)
                local r = 0
                forvalues i = 1/`nrows' {
                    if `_Est'[`i', 1] == `e1' {
                        local r = `r' + 1
                        matrix `_subset'[`r', 1] = `_Est'[`i', `z_col']
                        matrix `_subset'[`r', 2] = `_Est'[`i', `est_col']
                        matrix `_subset'[`r', 3] = `_Est'[`i', `panel_ci_lower_col']
                        matrix `_subset'[`r', 4] = `_Est'[`i', `panel_ci_upper_col']
                    }
                }

                preserve
                clear
                quietly svmat `_subset', names(col)
                rename c1 z
                rename c2 est
                rename c3 ci_lower
                rename c4 ci_upper

                // Build graph name (handle decimals/negatives)
                local e1_clean : subinstr local e1 "." "_", all
                local e1_clean : subinstr local e1_clean "-" "m", all
                local gname "eval`e1_clean'"

                capture noisily twoway                                        ///
                    (rarea ci_lower ci_upper z, color(gs12%60) lwidth(none))  ///
                    (line est z, lcolor(black) lwidth(vthick)),               ///
                    title("Evaluation point `e1'.")                           ///
                    xtitle("z") ytitle("")                                    ///
                    plotregion(fcolor(white))                                 ///
                    legend(off)                                               ///
                    name(`gname', replace)                                    ///
                    `pretrend_opt'                                            ///
                    `graph_opt'
                local graph_rc = _rc

                restore
                if `graph_rc' {
                    exit `graph_rc'
                }
                local graph_count = `graph_count' + 1
            }
        }

        di as text ""
        di as text "catt_gt_graph: `graph_count' Aggregated graph(s) generated (type = `agg_type')."
    }

    // =========================================================================
    // Step 9: Save graph if save_path() specified
    // =========================================================================
    if `"`save_path'"' != "" {
        local save_path_lc = lower(`"`save_path'"')
        if regexm(`"`save_path_lc'"', "\.gph$") {
            quietly graph save `"`save_path'"', replace
        }
        else {
            quietly graph export `"`save_path'"', replace
        }
        di as text "  Graph saved to: `save_path'"
    }

end
