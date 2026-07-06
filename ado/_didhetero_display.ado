*! Display CATT estimation results as a standard Stata estimation table

program define _didhetero_display
    version 16.0

    // Check that results exist
    capture confirm matrix e(results)
    if _rc {
        di as error "No estimation results found"
        exit 301
    }

    tempname results gteval_mat bw_mat
    matrix `results' = e(results)
    matrix `gteval_mat' = e(gteval)
    matrix `bw_mat' = e(bw)

    local nrows      = rowsof(`results')
    local num_gteval = e(num_gteval)
    local num_zeval  = e(num_zeval)
    local has_bstrap = e(bstrap)

    // Significance level and uniform-inference flag
    local level = e(level)
    if "`level'" == "" | `level' == . {
        local level = 95
    }
    local levdisp = strofreal(`level', "%9.0g")
    local uniform = e(uniformall)
    if "`uniform'" == "" | `uniform' == . {
        local uniform = 0
    }

    // Confidence-interval column label
    if `uniform' == 1 & `num_gteval' > 1 {
        local _ci_lab "[`levdisp'% Uniform CI]"
    }
    else {
        local _ci_lab "[`levdisp'% CI]"
    }

    // === Table header (standard Stata estimation-table framing) ===
    di as text "{hline 13}{c +}{hline 62}"
    di as text %12s "(g,t)   z" as text " {c |}" ///
       as text %11s "CATT" %11s "Std. err." "   " %-22s "`_ci_lab'" %10s "BW"
    di as text "{hline 13}{c +}{hline 62}"

    local row = 1
    forvalues gt = 1/`num_gteval' {
        local g1 = `gteval_mat'[`gt', 1]
        local t1 = `gteval_mat'[`gt', 2]
        local bw1 = `bw_mat'[1, `gt']

        // Group label row (dimension only; statistics blank)
        local _gtlab "(`g1',`t1')"
        di as text %-12s "`_gtlab'" as text " {c |}"

        forvalues r = 1/`num_zeval' {
            local z_val  = `results'[`row', 3]
            local est    = `results'[`row', 4]
            local se     = `results'[`row', 5]
            local ci1_l  = `results'[`row', 6]
            local ci1_u  = `results'[`row', 7]

            // Choose CI: prefer bootstrap (ci2), fallback to analytical (ci1)
            if `has_bstrap' == 1 {
                local ci2_l = `results'[`row', 8]
                local ci2_u = `results'[`row', 9]
                local ci_lo = cond(`ci2_l' != ., `ci2_l', `ci1_l')
                local ci_hi = cond(`ci2_u' != ., `ci2_u', `ci1_u')
            }
            else {
                local ci_lo = `ci1_l'
                local ci_hi = `ci1_u'
            }

            // Build a fixed-width bracketed CI cell
            local _b1 : di %9.4f `ci_lo'
            local _b2 : di %9.4f `ci_hi'
            local ci_str "[`_b1', `_b2']"

            di as result %11.3f `z_val' as text "  {c |}" ///
               as result %11.4f `est' %11.4f `se' ///
               as text "   `ci_str'" ///
               as result %10.4f `bw1'

            local row = `row' + 1
        }
    }

    di as text "{hline 13}{c +}{hline 62}"

    // === Critical value footer ===
    capture confirm matrix e(c_hat)
    local has_ana = (!_rc)
    local has_boot = 0
    if `has_bstrap' == 1 {
        capture confirm matrix e(c_check)
        local has_boot = (!_rc)
    }

    if `has_ana' | `has_boot' {
        tempname c_hat_mat c_check_mat

        if `has_ana' {
            matrix `c_hat_mat' = e(c_hat)
        }
        if `has_boot' {
            matrix `c_check_mat' = e(c_check)
        }

        // Detect whether all critical values are identical across (g,t)
        // (missing values are treated as "same").
        local all_same_ana = 1
        local all_same_boot = 1
        if `has_ana' & `num_gteval' > 1 {
            local _cv_ref = `c_hat_mat'[1, 1]
            forvalues gt = 2/`num_gteval' {
                local _cv_cur = `c_hat_mat'[1, `gt']
                if (`_cv_ref' == . & `_cv_cur' != .) | (`_cv_ref' != . & `_cv_cur' == .) {
                    local all_same_ana = 0
                }
                else if `_cv_ref' != . & `_cv_cur' != . {
                    if abs(`_cv_cur' - `_cv_ref') > 1e-10 {
                        local all_same_ana = 0
                    }
                }
            }
        }
        if `has_boot' & `num_gteval' > 1 {
            local _cv_ref = `c_check_mat'[1, 1]
            forvalues gt = 2/`num_gteval' {
                local _cv_cur = `c_check_mat'[1, `gt']
                if (`_cv_ref' == . & `_cv_cur' != .) | (`_cv_ref' != . & `_cv_cur' == .) {
                    local all_same_boot = 0
                }
                else if `_cv_ref' != . & `_cv_cur' != . {
                    if abs(`_cv_cur' - `_cv_ref') > 1e-10 {
                        local all_same_boot = 0
                    }
                }
            }
        }

        // Suffix noting uniform coverage across (g,t) pairs
        local _uni_suffix ""
        if `uniform' == 1 & `num_gteval' > 1 {
            local _uni_suffix "  (uniform across all g,t)"
        }

        if `all_same_ana' & `all_same_boot' {
            // Single critical value line
            if `has_boot' {
                local _cv_boot = `c_check_mat'[1, 1]
                if `_cv_boot' != . {
                    di as text "Critical value: Bootstrap = " ///
                       as result %6.4f `_cv_boot' as text "`_uni_suffix'"
                }
            }
            if `has_ana' {
                local _cv_ana = `c_hat_mat'[1, 1]
                if `_cv_ana' != . {
                    di as text "Critical value: Analytical = " ///
                       as result %6.4f `_cv_ana' as text "`_uni_suffix'"
                }
            }
        }
        else {
            // Per-pair critical values
            di as text "Critical values (per (g,t) pair):"
            forvalues gt = 1/`num_gteval' {
                local g1 = `gteval_mat'[`gt', 1]
                local t1 = `gteval_mat'[`gt', 2]
                if `has_boot' & `has_ana' {
                    local _cv_a = `c_hat_mat'[1, `gt']
                    local _cv_b = `c_check_mat'[1, `gt']
                    di as text %-13s "  (`g1',`t1')" ///
                       as text "Boot=" as result %8.4f `_cv_b' ///
                       as text "  Analytical=" as result %8.4f `_cv_a'
                }
                else if `has_boot' {
                    local _cv_b = `c_check_mat'[1, `gt']
                    di as text %-13s "  (`g1',`t1')" ///
                       as text "Boot=" as result %8.4f `_cv_b'
                }
                else {
                    local _cv_a = `c_hat_mat'[1, `gt']
                    di as text %-13s "  (`g1',`t1')" ///
                       as text "Analytical=" as result %8.4f `_cv_a'
                }
            }
        }
    }

    di as text ""

end
