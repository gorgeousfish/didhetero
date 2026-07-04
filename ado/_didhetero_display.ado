*! Display CATT estimation results in formatted tables

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

    local nrows = rowsof(`results')
    local num_gteval = e(num_gteval)
    local num_zeval = e(num_zeval)
    local has_bstrap = e(bstrap)

    // === Task 2: Unified table with single header ===
    di as text ""
    di as text "{hline 72}"
    if `has_bstrap' == 1 {
        di as text "       z        est         se       [95% CI]              bw"
    }
    else {
        di as text "       z        est         se       [95% CI]              bw"
    }
    di as text "{hline 72}"

    local row = 1
    forvalues gt = 1/`num_gteval' {
        local g1 = `gteval_mat'[`gt', 1]
        local t1 = `gteval_mat'[`gt', 2]
        local bw1 = `bw_mat'[1, `gt']

        // Lightweight group separator
        di as text " g=`g1', t=`t1'"

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

            // Format CI as [lo, hi] or [   .  ,    .  ]
            if `ci_lo' == . | `ci_hi' == . {
                local ci_str "[    .   ,    .   ]"
            }
            else {
                local _ci_lo_s : di %7.3f `ci_lo'
                local _ci_hi_s : di %7.3f `ci_hi'
                local ci_str "[`_ci_lo_s',`_ci_hi_s']"
            }

            di as result %8.4f `z_val' %10.4f `est' %10.4f `se' ///
               as text "  `ci_str'" ///
               as result %10.4f `bw1'

            local row = `row' + 1
        }
    }

    di as text "{hline 72}"

    // === Task 4: Critical Values matrix display ===
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

        // Check if all CV values are the same (uniformall case)
        // Note: missing values are treated as "same" (all missing = uniform)
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

        if `all_same_ana' & `all_same_boot' {
            // All CVs identical: single-line display
            local _cv_parts ""
            if `has_ana' {
                local _cv_ana = `c_hat_mat'[1, 1]
                if `_cv_ana' != . {
                    local _cv_ana_s : di %7.4f `_cv_ana'
                    local _cv_parts "Analytical=`_cv_ana_s'"
                }
            }
            if `has_boot' {
                local _cv_boot = `c_check_mat'[1, 1]
                if `_cv_boot' != . {
                    local _cv_boot_s : di %7.4f `_cv_boot'
                    if "`_cv_parts'" != "" {
                        local _cv_parts "`_cv_parts' | Bootstrap=`_cv_boot_s'"
                    }
                    else {
                        local _cv_parts "Bootstrap=`_cv_boot_s'"
                    }
                }
            }
            if "`_cv_parts'" != "" {
                di as text "Critical values: `_cv_parts'"
            }
        }
        else {
            // Different CVs: compact matrix
            di as text "Critical values:"
            if `has_ana' & `has_boot' {
                di as text "  (g,t)          Analytical  Bootstrap"
                forvalues gt = 1/`num_gteval' {
                    local g1 = `gteval_mat'[`gt', 1]
                    local t1 = `gteval_mat'[`gt', 2]
                    local _cv_a = `c_hat_mat'[1, `gt']
                    local _cv_b = `c_check_mat'[1, `gt']
                    di as text "  (`g1',`t1')" as result %12.4f `_cv_a' %12.4f `_cv_b'
                }
            }
            else if `has_ana' {
                di as text "  (g,t)          Analytical"
                forvalues gt = 1/`num_gteval' {
                    local g1 = `gteval_mat'[`gt', 1]
                    local t1 = `gteval_mat'[`gt', 2]
                    local _cv_a = `c_hat_mat'[1, `gt']
                    di as text "  (`g1',`t1')" as result %12.4f `_cv_a'
                }
            }
            else {
                di as text "  (g,t)          Bootstrap"
                forvalues gt = 1/`num_gteval' {
                    local g1 = `gteval_mat'[`gt', 1]
                    local t1 = `gteval_mat'[`gt', 2]
                    local _cv_b = `c_check_mat'[1, `gt']
                    di as text "  (`g1',`t1')" as result %12.4f `_cv_b'
                }
            }
        }
    }

    di as text ""

    // Check if ALL confidence intervals are missing (both ci1 and ci2)
    // Only show note if user has NO usable CI at all
    tempname _tmp_res
    matrix `_tmp_res' = e(results)
    local _all_ci_miss = 1
    forvalues i = 1/`nrows' {
        // Check ci1 (cols 6-7) and ci2 (cols 8-9 if bootstrap)
        if `_tmp_res'[`i', 6] != . | `_tmp_res'[`i', 7] != . {
            local _all_ci_miss = 0
        }
        if `has_bstrap' == 1 {
            if `_tmp_res'[`i', 8] != . | `_tmp_res'[`i', 9] != . {
                local _all_ci_miss = 0
            }
        }
    }
    if `_all_ci_miss' {
        di as text "Note: All confidence intervals are missing."
        di as text "  For bootstrap uniform confidence bands, add: biters(999)"
        di as text "  For analytical UCB, use more evaluation points in zeval()."
        di as text ""
    }

end
