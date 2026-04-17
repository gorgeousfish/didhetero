// ============================================================================
// example_basic.do
// Basic usage example for the didhetero Stata package
//
// Implements Imai, Qin, and Yanagi (2025)
// "Doubly Robust Uniform Confidence Bands for Group-Time Conditional Average
// Treatment Effects in Difference-in-Differences"
//
// Empirical application: minimum wage increases and teen employment
// Dataset: min_wage_cs.dta — county-level panel, 2001–2007
//          Callaway & Sant'Anna (2021) / Imai, Qin & Yanagi (2025)
//
// Research question:
//   Does the employment effect of minimum wage increases vary with
//   the county poverty rate (z)?
//
// Date:    2026-04-12
// Author:  didhetero development team
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare the minimum wage dataset
// ============================================================================
// min_wage_cs.dta is installed with the package.
// Variables:
//   lemp        — log teen employment (outcome)
//   first_treat — first year of minimum wage increase (0 = never-treated;
//                 2004, 2006, or 2007 = treated)
//   year        — calendar year (2001–2007)
//   countyreal  — county FIPS identifier
//   pov         — county poverty rate (continuous covariate z)

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

// Inspect the data
describe
summarize lemp G period Z

// Examine treatment cohort and time period structure
tab G
tab period

// ============================================================================
// Step 2: Estimate CATT — poverty rate as the effect modifier
// ============================================================================
// Evaluate the CATT at three representative poverty-rate values:
//   z = 0.105  (25th percentile)
//   z = 0.136  (50th percentile / median)
//   z = 0.181  (75th percentile)
//
// Restrict to the three instantaneous treatment effects (g=t) to keep
// computation manageable. Uses Gaussian kernel, IMSE-optimal bandwidth,
// and a second-order local polynomial.

local inference_opts "bstrap(false) uniformall(false)"

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

// Display the full results matrix
// Columns: g, t, z, estimate, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
matrix list e(results), format(%9.4f)

// ============================================================================
// Step 3: Sensitivity check — manual bandwidth
// ============================================================================
// Re-estimate with a fixed bandwidth to assess robustness.
// A bandwidth of 0.04 is narrow relative to the poverty rate support [0,0.47].

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    porder(2) kernel("gau") bwselect("manual") bw(0.04) ///
    `inference_opts' ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

// ============================================================================
// Step 4: Visualize CATT estimates
// ============================================================================
// Restore the IMSE-optimal estimates and plot.
// Each panel corresponds to one (g,t) pair; the x-axis is the poverty rate.

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated")

catt_gt_graph

// ============================================================================
// Step 5: Event-study aggregation and visualization
// ============================================================================
// Aggregate group-time CATT into event-time parameters.
// e=0: year of the minimum wage increase
// e=1: one year after the increase
// e=2: two years after the increase

// (a) Without bootstrap — fast analytical confidence bands
aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")

// (b) With bootstrap — preferred for publication (may take a few minutes)
aggte_gt, type("dynamic") eval(0 1 2) bstrap("true") biters(500) seed(42)
catt_gt_graph, plot_type("Aggregated")
