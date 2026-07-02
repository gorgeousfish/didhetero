*! catt_gt_estat.ado
*! Post-estimation diagnostics for catt_gt
*! Diagnostic display tools - NOT formal statistical tests
*!
*! Subcommands:
*!   estat overlap   - Display GPS distribution diagnostics (Assumption 2.1)
*!   estat pretrend  - Display pre-treatment CATT estimates and UCB coverage

program catt_gt_estat, rclass
    version 16.0

    // Verify that catt_gt or didhetero was the last estimation command
    if "`e(cmd)'" != "catt_gt" & "`e(cmd)'" != "didhetero" {
        di as error "estat requires that catt_gt or didhetero be the last estimation command"
        error 301
    }

    gettoken subcmd 0 : 0, parse(" ,")

    if "`subcmd'" == "overlap" {
        _catt_gt_estat_overlap `0'
    }
    else if "`subcmd'" == "pretrend" {
        _catt_gt_estat_pretrend `0'
    }
    else {
        di as error "estat `subcmd' is not recognized"
        di as error "Available subcommands after catt_gt:"
        di as error "    estat overlap   - GPS distribution diagnostics"
        di as error "    estat pretrend  - Pre-treatment period diagnostics"
        exit 198
    }
end

// =============================================================================
// estat overlap
// Displays descriptive statistics of GPS values to diagnose overlap
// (Assumption 2.1: P(G=g|X) bounded away from 0 and 1)
//
// This is a DIAGNOSTIC TOOL, not a formal statistical test.
// =============================================================================
program _catt_gt_estat_overlap, rclass
    version 16.0
    syntax [, THreshold(real 0.01)]

    // Validate threshold
    if `threshold' <= 0 | `threshold' >= 0.5 {
        di as error "threshold() must be between 0 and 0.5"
        exit 198
    }

    // Retrieve GPS matrix from e()
    capture confirm matrix e(dh_gps_mat)
    if _rc {
        di as error "GPS matrix not found in estimation results"
        di as error "This may occur if the estimation was run without storing GPS values"
        exit 301
    }

    tempname gps_full
    matrix `gps_full' = e(dh_gps_mat)

    local ncols = colsof(`gps_full')
    local nrows = rowsof(`gps_full')

    // Determine control group type from column count
    // nevertreated: (id, g, p_hat) = 3 cols
    // notyettreated: (id, g, t, p_hat) = 4 cols
    local control_group "`e(control_group)'"

    if `ncols' == 3 {
        local pcol = 3
        local gcol = 2
    }
    else if `ncols' == 4 {
        local pcol = 4
        local gcol = 2
    }
    else {
        di as error "Unexpected GPS matrix dimensions"
        exit 498
    }

    // =========================================================================
    // Header
    // =========================================================================
    di as text ""
    di as text "{hline 72}"
    di as text "  Overlap Diagnostics (GPS Distribution)"
    di as text "  NOTE: This is a diagnostic tool, not a formal test."
    di as text "        (Assumption 2.1 requires GPS bounded away from 0 and 1)"
    di as text "{hline 72}"
    di as text ""
    di as text "  Control group: `control_group'"
    di as text "  Boundary threshold: `threshold'"
    di as text ""

    // =========================================================================
    // Overall GPS descriptive statistics
    // =========================================================================
    // Extract p_hat column into a Stata variable for computation
    preserve
    quietly {
        clear
        set obs `nrows'
        gen double gps_val = .
        gen double group_val = .
        forvalues i = 1/`nrows' {
            replace gps_val = `gps_full'[`i', `pcol'] in `i'
            replace group_val = `gps_full'[`i', `gcol'] in `i'
        }
    }

    // Overall statistics
    quietly summarize gps_val, detail
    local gps_min = r(min)
    local gps_p5  = r(p5)
    local gps_p25 = r(p25)
    local gps_p50 = r(p50)
    local gps_p75 = r(p75)
    local gps_p95 = r(p95)
    local gps_max = r(max)
    local gps_mean = r(mean)
    local gps_n   = r(N)

    di as text "  Overall GPS Distribution (N = `gps_n'):"
    di as text "  {hline 50}"
    di as text "    Min       P5        P25       Median    P75       P95       Max"
    di as result "   " %8.5f `gps_min' ///
       "  " %8.5f `gps_p5' ///
       "  " %8.5f `gps_p25' ///
       "  " %8.5f `gps_p50' ///
       "  " %8.5f `gps_p75' ///
       "  " %8.5f `gps_p95' ///
       "  " %8.5f `gps_max'
    di as text ""

    // =========================================================================
    // Boundary observations
    // =========================================================================
    quietly count if gps_val < `threshold'
    local n_low = r(N)
    quietly count if gps_val > (1 - `threshold')
    local n_high = r(N)
    local n_extreme = `n_low' + `n_high'
    local pct_extreme = 100 * `n_extreme' / `gps_n'

    local upper_thresh = 1 - `threshold'
    di as text "  Observations near boundaries (GPS < `threshold' or GPS > `upper_thresh'):"
    di as text "    GPS < `threshold':          " as result `n_low'
    di as text "    GPS > `upper_thresh':          " as result `n_high'
    di as text "    Total extreme:        " as result `n_extreme' ///
       as text " (" as result %4.1f `pct_extreme' as text "% of observations)"
    di as text ""

    if `n_extreme' > 0 {
        if `pct_extreme' > 25 {
            di as error "  WARNING: >25% of GPS values are near boundaries."
            di as error "  This suggests potential overlap violation (Assumption 2.1)."
            di as error "  Results may be unreliable."
        }
        else if `pct_extreme' > 10 {
            di as text "  {res}Note:{txt} >10% of GPS values are near boundaries."
            di as text "  Consider examining covariate balance more carefully."
        }
    }
    else {
        di as text "  No observations near boundaries - overlap appears satisfactory."
    }
    di as text ""

    // =========================================================================
    // By-group GPS distribution
    // =========================================================================
    di as text "  GPS Distribution by Treatment Group:"
    di as text "  {hline 66}"
    di as text "    Group       N      Min       P25       Median    P75       Max"
    di as text "  {hline 66}"

    quietly levelsof group_val, local(groups)
    foreach g of local groups {
        quietly summarize gps_val if group_val == `g', detail
        local g_n    = r(N)
        local g_min  = r(min)
        local g_p25  = r(p25)
        local g_p50  = r(p50)
        local g_p75  = r(p75)
        local g_max  = r(max)
        di as text "   " %6.0f `g' as result ///
           "  " %6.0f `g_n' ///
           "  " %8.5f `g_min' ///
           "  " %8.5f `g_p25' ///
           "  " %8.5f `g_p50' ///
           "  " %8.5f `g_p75' ///
           "  " %8.5f `g_max'
    }
    di as text "  {hline 66}"

    restore

    // =========================================================================
    // Return values
    // =========================================================================
    return scalar gps_min = `gps_min'
    return scalar gps_p5  = `gps_p5'
    return scalar gps_p25 = `gps_p25'
    return scalar gps_p50 = `gps_p50'
    return scalar gps_p75 = `gps_p75'
    return scalar gps_p95 = `gps_p95'
    return scalar gps_max = `gps_max'
    return scalar gps_mean = `gps_mean'
    return scalar n_extreme = `n_extreme'
    return scalar pct_extreme = `pct_extreme'
    return scalar n_obs = `gps_n'
    return scalar threshold = `threshold'

    di as text ""
    di as text "{hline 72}"
    di as text "  Note: This is a diagnostic display, not a formal test."
    di as text "  No p-values or test statistics are produced."
    di as text "  Visual inspection of GPS distributions is recommended."
    di as text "{hline 72}"
    di as text ""
end

// =============================================================================
// estat pretrend
// Displays pre-treatment CATT estimates and whether UCBs cover zero
// Based on Section 5/Appendix of the paper (Eq. S5.1-S5.2):
//   CATT_{g,t}(z) = 0 for t <= g - 2 (pre-treatment periods)
//
// This is a DIAGNOSTIC REFERENCE, not a formal pre-test.
// The paper explicitly warns against using this as a formal test
// (Manuscript lines 2942-2944).
// =============================================================================
program _catt_gt_estat_pretrend, rclass
    version 16.0
    syntax

    // Check pretrend flag
    if e(pretrend) != 1 {
        di as error "Pre-trend diagnostics require catt_gt to be run with the pretrend option"
        di as error "Usage: catt_gt ..., pretrend"
        exit 198
    }

    // Retrieve results
    tempname results gteval_mat
    matrix `results' = e(results)
    matrix `gteval_mat' = e(gteval)

    local nrows = rowsof(`results')
    local num_gteval = e(num_gteval)
    local num_zeval = e(num_zeval)
    local anticip = e(anticip)
    local has_bstrap = e(bstrap)
    local alp_est = e(alp)

    // =========================================================================
    // Header
    // =========================================================================
    di as text ""
    di as text "{hline 72}"
    di as text "  Pre-trend Diagnostics (Pre-treatment CATT Estimates)"
    di as text "  NOTE: This is a diagnostic reference, NOT a formal pre-test."
    di as text "        The paper recommends jointly examining pre and post results"
    di as text "        rather than conditioning on a pre-test outcome."
    di as text "{hline 72}"
    di as text ""
    di as text "  Significance level: `alp_est'"
    di as text "  Anticipation: `anticip'"
    if `has_bstrap' == 1 {
        di as text "  Confidence bands: Uniform (bootstrap)"
    }
    else {
        di as text "  Confidence bands: Pointwise (analytical)"
    }
    di as text ""

    // =========================================================================
    // Display pre-treatment period results
    // Pre-treatment: t < g - anticipation (treatment has not yet occurred)
    // =========================================================================
    di as text "  Pre-treatment period estimates (t < g - anticipation):"
    di as text "  Under the null: CATT_{g,t}(z) = 0 for pre-treatment periods"
    di as text ""

    if `has_bstrap' == 1 {
        di as text "     g        t        z       est       se      UCB_low   UCB_up  Covers 0?"
        di as text "  {hline 72}"
    }
    else {
        di as text "     g        t        z       est       se      CI_low    CI_up   Covers 0?"
        di as text "  {hline 66}"
    }

    // Loop through results and identify pre-treatment rows
    local row = 1
    local n_pretreat = 0
    local n_reject = 0

    forvalues gt = 1/`num_gteval' {
        local g1 = `gteval_mat'[`gt', 1]
        local t1 = `gteval_mat'[`gt', 2]

        // Pre-treatment condition: t < g - anticipation
        local is_pre = (`t1' < `g1' - `anticip')

        forvalues r = 1/`num_zeval' {
            if `is_pre' {
                local z_val = `results'[`row', 3]
                local est   = `results'[`row', 4]
                local se    = `results'[`row', 5]

                // Use uniform CB (ci2) if bootstrap available, else pointwise (ci1)
                if `has_bstrap' == 1 {
                    local cb_l = `results'[`row', 8]
                    local cb_u = `results'[`row', 9]
                }
                else {
                    local cb_l = `results'[`row', 6]
                    local cb_u = `results'[`row', 7]
                }

                // Check if confidence band covers zero
                local covers_zero = (`cb_l' <= 0 & `cb_u' >= 0)
                if `covers_zero' {
                    local flag "  Yes"
                }
                else {
                    local flag "  No *"
                    local n_reject = `n_reject' + 1
                }

                di as result "  " %6.0f `g1' ///
                   "   " %6.0f `t1' ///
                   "   " %7.4f `z_val' ///
                   "  " %8.4f `est' ///
                   "  " %8.4f `se' ///
                   "  " %8.4f `cb_l' ///
                   "  " %8.4f `cb_u' ///
                   as text "`flag'"

                local n_pretreat = `n_pretreat' + 1
            }
            local row = `row' + 1
        }
    }

    if `n_pretreat' == 0 {
        di as text "  No pre-treatment periods found in estimation results."
        di as text ""
        di as text "  This typically means all evaluated (g,t) pairs satisfy t >= g - anticipation."
        di as text "  To include pre-treatment periods, re-run catt_gt with the pretrend option"
        di as text "  and ensure the (g,t) domain includes pre-treatment periods."
    }
    else {
        di as text ""
        di as text "  {hline 50}"
        di as text "  Total pre-treatment (g,t,z) cells: " as result `n_pretreat'
        di as text "  Cells where CB does not cover 0:   " as result `n_reject'
        if `n_reject' > 0 {
            di as text "  (* denotes confidence band excludes zero)"
        }
    }

    // =========================================================================
    // Return values
    // =========================================================================
    return scalar n_pretreat = `n_pretreat'
    return scalar n_reject = `n_reject'
    return scalar alpha = `alp_est'
    return scalar anticipation = `anticip'

    // =========================================================================
    // Footer disclaimer
    // =========================================================================
    di as text ""
    di as text "{hline 72}"
    di as text "  IMPORTANT: This display is for diagnostic reference only."
    di as text "  The paper does NOT recommend using pre-trend results as a formal"
    di as text "  pre-test for conditioning subsequent inference."
    di as text "  Instead, jointly report pre-treatment and post-treatment estimates"
    di as text "  to allow readers to assess the plausibility of parallel trends."
    di as text "{hline 72}"
    di as text ""
end
