// This example requires: ssc install didhetero
// ============================================================================
// example_advanced.do — Advanced options for the didhetero Stata package
//
// Demonstrates side-by-side comparisons of estimation options.
// Each section changes exactly one option to isolate its effect.
// Dataset: min_wage_cs.dta (Callaway & Sant'Anna 2021)
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare data
// ============================================================================

use "../data/min_wage_cs.dta", clear
rename (first_treat year countyreal) (G period id)
bysort id (period): gen Z = pov[1]

summarize lemp G period Z

// Common options
local zeval_pts "0.105 0.136 0.181"
local gtv "gteval(2004 2004 2004 2005 2006 2006 2006 2007 2007 2007)"

// ============================================================================
// Section A: Kernel function comparison
// ============================================================================
// Gaussian: smooth infinite-support weights (default)
// Epanechnikov: compact-support, MSE-optimal among second-order kernels

display _newline(2)
display as text "=============================================="
display as text " Section A: Kernel function comparison"
display as text "=============================================="

display _newline
display as text "--- Gaussian kernel ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_gau = e(results)

display _newline
display as text "--- Epanechnikov kernel ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("epa") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_epa = e(results)

display _newline
display as text "--- Comparison: Gaussian vs Epanechnikov ---"
display as text "  Gaussian estimate | Epa estimate | Gaussian bw | Epa bw"

local nrows = rowsof(est_gau)
forvalues i = 1/`nrows' {
    local g   = est_gau[`i', 1]
    local t   = est_gau[`i', 2]
    local z   = est_gau[`i', 3]
    local e1  = est_gau[`i', 4]
    local e2  = est_epa[`i', 4]
    local b1  = est_gau[`i', 10]
    local b2  = est_epa[`i', 10]
    display as text "  g=`g' t=`t' z=" %6.3f `z' ///
        "  est_gau=" %8.4f `e1' "  est_epa=" %8.4f `e2' ///
        "  bw_gau=" %6.4f `b1' "  bw_epa=" %6.4f `b2'
}


// ============================================================================
// Section B: Bandwidth selection method comparison
// ============================================================================
// IMSE1: rule-of-thumb (default)
// IMSE2: plug-in estimator (more accurate, heavier)
// US1:   uniform-smoothing (preferred for uniform confidence bands)

display _newline(2)
display as text "=============================================="
display as text " Section B: Bandwidth selection comparison"
display as text "=============================================="

display _newline
display as text "--- IMSE1 (rule-of-thumb) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_imse1 = e(results)

display _newline
display as text "--- IMSE2 (plug-in) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE2") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_imse2 = e(results)

display _newline
display as text "--- US1 (uniform smoothing) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("US1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_us1 = e(results)

display _newline
display as text "--- Comparison: Bandwidth selection methods ---"
display as text "  IMSE1 bw | IMSE2 bw | US1 bw"

local nrows = rowsof(est_imse1)
forvalues i = 1/`nrows' {
    local g   = est_imse1[`i', 1]
    local t   = est_imse1[`i', 2]
    local z   = est_imse1[`i', 3]
    local b1  = est_imse1[`i', 10]
    local b2  = est_imse2[`i', 10]
    local b3  = est_us1[`i', 10]
    display as text "  g=`g' t=`t' z=" %6.3f `z' ///
        "  bw_IMSE1=" %6.4f `b1' ///
        "  bw_IMSE2=" %6.4f `b2' ///
        "  bw_US1="   %6.4f `b3'
}


// ============================================================================
// Section C: Manual bandwidth
// ============================================================================
// Fixed bandwidth via bwselect("manual") + bw() for robustness checks.

display _newline(2)
display as text "=============================================="
display as text " Section C: Manual bandwidth"
display as text "=============================================="

display _newline
display as text "--- Manual bandwidth = 0.04 (narrow) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("manual") bw(0.04) ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_narrow = e(results)

display _newline
display as text "--- Manual bandwidth = 0.10 (moderate) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("manual") bw(0.10) ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_moderate = e(results)

display _newline
display as text "--- Comparison: narrow vs moderate bandwidth ---"

local nrows = rowsof(est_narrow)
forvalues i = 1/`nrows' {
    local g   = est_narrow[`i', 1]
    local t   = est_narrow[`i', 2]
    local z   = est_narrow[`i', 3]
    local e1  = est_narrow[`i', 4]
    local e2  = est_moderate[`i', 4]
    display as text "  g=`g' t=`t' z=" %6.3f `z' ///
        "  est_narrow=" %8.4f `e1' "  est_moderate=" %8.4f `e2'
}


// ============================================================================
// Section D: Control group comparison
// ============================================================================
// notyettreated: never-treated + not-yet-treated (larger pool, lower variance)
// nevertreated:  never-treated only (more conservative)

display _newline(2)
display as text "=============================================="
display as text " Section D: Control group comparison"
display as text "=============================================="

display _newline
display as text "--- Not-yet-treated control group ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_nyt = e(results)

display _newline
display as text "--- Never-treated control group ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("nevertreated")

matrix est_nt = e(results)

display _newline
display as text "--- Comparison: Not-yet-treated vs Never-treated ---"

local nrows = rowsof(est_nyt)
forvalues i = 1/`nrows' {
    local g    = est_nyt[`i', 1]
    local t    = est_nyt[`i', 2]
    local z    = est_nyt[`i', 3]
    local e1   = est_nyt[`i', 4]
    local e2   = est_nt[`i', 4]
    local diff = `e1' - `e2'
    display as text "  g=`g' t=`t' z=" %6.3f `z' ///
        "  est_nyt=" %8.4f `e1' "  est_nt=" %8.4f `e2' ///
        "  diff=" %8.4f `diff'
}


// ============================================================================
// Section E: Polynomial order comparison
// ============================================================================
// porder(1): local linear - lower variance, higher bias
// porder(2): local quadratic - lower bias, higher variance

display _newline(2)
display as text "=============================================="
display as text " Section E: Polynomial order comparison"
display as text "=============================================="

display _newline
display as text "--- Local linear (porder = 1) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_p1 = e(results)

display _newline
display as text "--- Local quadratic (porder = 2) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix est_p2 = e(results)

display _newline
display as text "--- Comparison: porder(1) vs porder(2) ---"

local nrows = rowsof(est_p1)
forvalues i = 1/`nrows' {
    local g   = est_p1[`i', 1]
    local t   = est_p1[`i', 2]
    local z   = est_p1[`i', 3]
    local e1  = est_p1[`i', 4]
    local e2  = est_p2[`i', 4]
    local b1  = est_p1[`i', 10]
    local b2  = est_p2[`i', 10]
    display as text "  g=`g' t=`t' z=" %6.3f `z' ///
        "  est_p1=" %8.4f `e1' "  est_p2=" %8.4f `e2' ///
        "  bw_p1=" %6.4f `b1' "  bw_p2=" %6.4f `b2'
}


// ============================================================================
// Section F: GPS Trimming (Assumption 3 enforcement)
// ============================================================================
// Prevent extreme propensity scores from inflating weights

display _newline(2)
display as text "=============================================="
display as text " Section F: GPS Trimming"
display as text "=============================================="

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) nobstrap nouniformall ///
    control_group("notyettreated") ///
    gpstrim(0.01 0.99)

matrix list e(gps_diagnostics), format(%9.4f)


// ============================================================================
// Section G: Robust Bias-Corrected Inference (Section 3.3)
// ============================================================================
// RBC: uses LLR bandwidth with LQR estimation for bias correction

display _newline(2)
display as text "=============================================="
display as text " Section G: RBC Inference"
display as text "=============================================="

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    nobstrap nouniformall ///
    control_group("notyettreated") rbc

matrix list e(results), format(%9.4f)


// ============================================================================
// Section H: Undersmoothing Adjustment (Assumption 4(iii))
// ============================================================================
// Automatically shrinks bandwidth to satisfy undersmoothing condition

display _newline(2)
display as text "=============================================="
display as text " Section H: Undersmoothing"
display as text "=============================================="

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) nobstrap nouniformall ///
    control_group("notyettreated") undersmooth

matrix list e(results), format(%9.4f)


// ============================================================================
// Section I: Performance Profiling
// ============================================================================
// Display computation time breakdown by stage

display _newline(2)
display as text "=============================================="
display as text " Section I: Performance Profiling"
display as text "=============================================="

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) nobstrap nouniformall ///
    control_group("notyettreated") profile


// ============================================================================
// Section J: Anticipation Effects
// ============================================================================
// Allow 1 period of anticipation before treatment

display _newline(2)
display as text "=============================================="
display as text " Section J: Anticipation Effects"
display as text "=============================================="

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) nobstrap nouniformall ///
    control_group("notyettreated") anticipation(1)

matrix list e(results), format(%9.4f)


// ============================================================================
// Summary
// ============================================================================
display _newline(2)
display as text "=============================================="
display as text " Advanced options example complete"
display as text "=============================================="
