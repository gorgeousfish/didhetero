// ============================================================================
// example_advanced.do
// Advanced options example for the didhetero Stata package
//
// Implements Imai, Qin, and Yanagi (2025)
// "Doubly Robust Uniform Confidence Bands for Group-Time Conditional Average
// Treatment Effects in Difference-in-Differences"
//
// Empirical application: minimum wage increases and teen employment
// Dataset: min_wage_cs.dta — county-level panel, 2001–2007
//          Callaway & Sant'Anna (2021) / Imai, Qin & Yanagi (2025)
//
// Demonstrates how different estimation options affect CATT estimates by
// comparing settings side-by-side on the minimum wage data. Each section
// changes exactly one option so the effect is isolated.
//
// NOTE: Every call below explicitly disables bootstrap and uniform bands so
// the example remains directly runnable. The focus is on how estimation
// choices move point estimates and bandwidths.
//
// Date:    2026-04-12
// Author:  didhetero development team
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare the minimum wage dataset
// ============================================================================

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

describe
summarize lemp G period Z

// Common options shared across all comparisons below
// Evaluation points: IQR of the poverty rate (z)
local zeval_pts "0.105 0.136 0.181"
local gtv        "gteval(2004 2004 2006 2006 2007 2007)"
local inference_opts "bstrap(false) uniformall(false)"


// ============================================================================
// Section A: Kernel function comparison
// ============================================================================
// The kernel function controls how nearby observations are weighted when
// estimating the local polynomial regression at each evaluation point.
//
// Gaussian ("gau"):
//   Smooth, infinitely differentiable weights that decay exponentially.
//   All observations receive positive weight. Generally produces smoother
//   estimates and is the default.
//
// Epanechnikov ("epa"):
//   Compact-support kernel that assigns zero weight to observations
//   beyond the bandwidth. Theoretically MSE-optimal among second-order
//   kernels but can produce less smooth estimates.

display _newline(2)
display as text "=============================================="
display as text " Section A: Kernel function comparison"
display as text "=============================================="

// --- A1: Gaussian kernel ---
display _newline
display as text "--- Gaussian kernel ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_gau = e(results)

// --- A2: Epanechnikov kernel ---
display _newline
display as text "--- Epanechnikov kernel ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("epa") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_epa = e(results)

// Compare point estimates (column 4) and bandwidths (column 10)
display _newline
display as text "--- Comparison: Gaussian vs Epanechnikov ---"
display as text "Rows: (g, t, z) combinations.  Cols: estimate, bandwidth."
display as text ""
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
// The bandwidth controls the bias-variance trade-off in local polynomial
// estimation. Larger bandwidths reduce variance but increase bias.
//
// IMSE1: Integrated MSE-optimal bandwidth, rule-of-thumb (default).
// IMSE2: Integrated MSE-optimal bandwidth, plug-in estimator. More accurate
//        in finite samples but computationally heavier.
// US1:   Uniform-smoothing bandwidth. Ensures uniform bias order across z,
//        preferred when constructing uniform confidence bands.

display _newline(2)
display as text "=============================================="
display as text " Section B: Bandwidth selection comparison"
display as text "=============================================="

// --- B1: IMSE1 ---
display _newline
display as text "--- IMSE1 (rule-of-thumb) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_imse1 = e(results)

// --- B2: IMSE2 ---
display _newline
display as text "--- IMSE2 (plug-in) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE2") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_imse2 = e(results)

// --- B3: US1 ---
display _newline
display as text "--- US1 (uniform smoothing) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("US1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_us1 = e(results)

// Compare bandwidths across methods
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
// Supply a fixed bandwidth via bwselect("manual") and bw() for robustness
// checks or to replicate a published specification.
// For the poverty rate z ∈ [0.019, 0.467], a bandwidth of 0.04 is narrow
// (high local resolution) while 0.10 is moderate.

display _newline(2)
display as text "=============================================="
display as text " Section C: Manual bandwidth"
display as text "=============================================="

display _newline
display as text "--- Manual bandwidth = 0.04 (narrow) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("manual") bw(0.04) ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_narrow = e(results)

display _newline
display as text "--- Manual bandwidth = 0.10 (moderate) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("manual") bw(0.10) ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_moderate = e(results)

// Verify bandwidths and compare estimates
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
// The control group determines the counterfactual for each treated cohort.
//
// Not-yet-treated ("notyettreated"):
//   Includes never-treated units and units not yet treated at period t.
//   Larger control pool (lower variance) but requires a broader parallel
//   trends assumption.
//
// Never-treated ("nevertreated"):
//   Restricts the control pool to counties that never received a minimum
//   wage increase. More conservative; parallel trends is a weaker assumption
//   but variance is higher. Applicable only when never-treated units exist.

display _newline(2)
display as text "=============================================="
display as text " Section D: Control group comparison"
display as text "=============================================="

// --- D1: Not-yet-treated ---
display _newline
display as text "--- Not-yet-treated control group ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_nyt = e(results)

// --- D2: Never-treated ---
display _newline
display as text "--- Never-treated control group ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("nevertreated")

matrix est_nt = e(results)

// Compare estimates
display _newline
display as text "--- Comparison: Not-yet-treated vs Never-treated ---"
display as text "  NYT estimate | NT estimate | Difference"

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
// The polynomial order (porder) controls the degree of the local polynomial.
//
// Local linear (porder = 1): Lower variance, higher bias for curved CATT.
// Local quadratic (porder = 2): Reduces bias at the cost of higher variance.
//   Preferred when the sample is large (as here: 2,284 counties) and the
//   CATT function may be nonlinear in the poverty rate.

display _newline(2)
display as text "=============================================="
display as text " Section E: Polynomial order comparison"
display as text "=============================================="

// --- E1: Local linear (porder = 1) ---
display _newline
display as text "--- Local linear (porder = 1) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(1) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_p1 = e(results)

// --- E2: Local quadratic (porder = 2) ---
display _newline
display as text "--- Local quadratic (porder = 2) ---"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(`zeval_pts') `gtv' ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

matrix est_p2 = e(results)

// Compare estimates and bandwidths
display _newline
display as text "--- Comparison: porder(1) vs porder(2) ---"
display as text "  p1 estimate | p2 estimate | p1 bw | p2 bw"

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
// Summary
// ============================================================================
display _newline(2)
display as text "=============================================="
display as text " Advanced options example complete"
display as text "=============================================="
display as text ""
display as text "Stored matrices for further analysis:"
display as text "  est_gau      — Gaussian kernel estimates"
display as text "  est_epa      — Epanechnikov kernel estimates"
display as text "  est_imse1    — IMSE1 bandwidth estimates"
display as text "  est_imse2    — IMSE2 bandwidth estimates"
display as text "  est_us1      — US1 bandwidth estimates"
display as text "  est_narrow   — Manual bandwidth (0.04) estimates"
display as text "  est_moderate — Manual bandwidth (0.10) estimates"
display as text "  est_nyt      — Not-yet-treated control group estimates"
display as text "  est_nt       — Never-treated control group estimates"
display as text "  est_p1       — Local linear (porder=1) estimates"
display as text "  est_p2       — Local quadratic (porder=2) estimates"
