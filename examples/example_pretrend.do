// This example requires: ssc install didhetero
// ============================================================================
// example_pretrend.do — Pre-trends testing for the didhetero Stata package
//
// The pretrend flag estimates CATT for pre-treatment periods (t < g).
// Under conditional parallel trends, these should be zero for all z.
// Dataset: min_wage_cs.dta (Callaway & Sant'Anna 2021)
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare data
// ============================================================================
// Pre-treatment structure:
//   Group g=2004: years 2001, 2002, 2003 are pre-treatment
//   Group g=2006: years 2001-2005 are pre-treatment
//   Group g=2007: years 2001-2006 are pre-treatment

use "../data/min_wage_cs.dta", clear
rename (first_treat year countyreal) (G period id)
bysort id (period): gen Z = pov[1]

summarize lemp G period Z

// ============================================================================
// Step 2: Run pre-trends test
// ============================================================================
// pretrend estimates CATT for t < g (excluding baseline t = g-1).
// When pretrend is specified, gteval() must be omitted; catt_gt automatically
// selects all valid pre-treatment (g,t) pairs.

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated") ///
    pretrend

matrix list e(results), format(%9.4f)
display as text "Pre-trend flag stored in e(pretrend): " e(pretrend)

// ============================================================================
// Step 3: Visualize pre-trends
// ============================================================================
// catt_gt_graph auto-detects e(pretrend)==1 and adds red dashed y=0 line.
// If confidence bands contain zero across all z, parallel trends supported.

catt_gt_graph

// ============================================================================
// Interpretation notes
// ============================================================================
// - Bands containing y=0: no evidence against parallel trends for that (g,t)
// - Bands excluding y=0: suggests differential pre-trends at those z values
// - Failure to reject != assumption holding (may reflect low power)
// - If pre-trends violated: try alternative control_group, add xformula(),
//   or vary bandwidth/porder for sensitivity analysis

display _newline
display as text "Pre-trends testing example complete."
display as text "  e(results)  - full results matrix"
display as text "  e(pretrend) - pre-trend flag (1 = pre-trends test)"
