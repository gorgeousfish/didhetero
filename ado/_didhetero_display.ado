*! _didhetero_display.ado
*! Display CATT estimation results in a formatted table
*! Shared by didhetero.ado and catt_gt.ado

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

    // Display per (g,t) pair
    local row = 1
    forvalues gt = 1/`num_gteval' {
        local g1 = `gteval_mat'[`gt', 1]
        local t1 = `gteval_mat'[`gt', 2]
        local bw1 = `bw_mat'[1, `gt']

        di as text ""
        di as text "{hline 72}"
        di as text "  Group g = " as result "`g1'" ///
           as text ",  Period t = " as result "`t1'" ///
           as text ",  Bandwidth = " as result %9.6f `bw1'
        di as text "{hline 72}"

        if `has_bstrap' == 1 {
            di as text "       z" ///
               as text "       est" ///
               as text "        se" ///
               as text "   ci1_low" ///
               as text "    ci1_up" ///
               as text "   ci2_low" ///
               as text "    ci2_up"
            di as text "{hline 72}"
        }
        else {
            di as text "       z" ///
               as text "       est" ///
               as text "        se" ///
               as text "   ci1_low" ///
               as text "    ci1_up"
            di as text "{hline 52}"
        }

        forvalues r = 1/`num_zeval' {
            local z_val  = `results'[`row', 3]
            local est    = `results'[`row', 4]
            local se     = `results'[`row', 5]
            local ci1_l  = `results'[`row', 6]
            local ci1_u  = `results'[`row', 7]

            if `has_bstrap' == 1 {
                local ci2_l = `results'[`row', 8]
                local ci2_u = `results'[`row', 9]

                di as result %8.4f `z_val' ///
                   as result %10.4f `est' ///
                   as result %10.4f `se' ///
                   as result %10.4f `ci1_l' ///
                   as result %10.4f `ci1_u' ///
                   as result %10.4f `ci2_l' ///
                   as result %10.4f `ci2_u'
            }
            else {
                di as result %8.4f `z_val' ///
                   as result %10.4f `est' ///
                   as result %10.4f `se' ///
                   as result %10.4f `ci1_l' ///
                   as result %10.4f `ci1_u'
            }

            local row = `row' + 1
        }
    }

    // Display critical values
    di as text ""
    di as text "{hline 72}"
    di as text "Critical values:"

    capture confirm matrix e(c_hat)
    if !_rc {
        tempname c_hat_mat
        matrix `c_hat_mat' = e(c_hat)
        forvalues gt = 1/`num_gteval' {
            local g1 = `gteval_mat'[`gt', 1]
            local t1 = `gteval_mat'[`gt', 2]
            local cv = `c_hat_mat'[1, `gt']
            di as text "  Analytical (g=`g1', t=`t1'): " as result %9.6f `cv'
        }
    }

    if `has_bstrap' == 1 {
        capture confirm matrix e(c_check)
        if !_rc {
            tempname c_check_mat
            matrix `c_check_mat' = e(c_check)
            forvalues gt = 1/`num_gteval' {
                local g1 = `gteval_mat'[`gt', 1]
                local t1 = `gteval_mat'[`gt', 2]
                local cv = `c_check_mat'[1, `gt']
                di as text "  Bootstrap  (g=`g1', t=`t1'): " as result %9.6f `cv'
            }
        }
    }

    di as text "{hline 72}"
    di as text ""

end
