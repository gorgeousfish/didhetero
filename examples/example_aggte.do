// This example requires: ssc install didhetero
// ============================================================================
// example_aggte.do — Aggregation parameters for the didhetero Stata package
//
// Four aggregation types: dynamic, group, calendar, simple
// Dataset: min_wage_cs.dta (Callaway & Sant'Anna 2021)
// ============================================================================

clear all
set more off

// -------------------------------------------------------------------------- //
// 1. Load and prepare data
// -------------------------------------------------------------------------- //

use "../data/min_wage_cs.dta", clear
rename (first_treat year countyreal) (G period id)
bysort id (period): gen Z = pov[1]

// Evaluation points: IQR of the county poverty rate
local zv   "zeval(0.105 0.136 0.181)"
local opts "porder(2) kernel(gau) bwselect(IMSE1) control_group(notyettreated)"

// Estimate CATT for representative (g,t) pairs covering event times 0-3
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `opts' ///
    gteval(2004 2004 2004 2005 2004 2006 2004 2007 2006 2006 2006 2007 2007 2007) ///
    nobstrap nouniformall

// -------------------------------------------------------------------------- //
// 2. Dynamic aggregation (event study)
//    Aggregates CATT(g,t,z) by event time e = t - g.
// -------------------------------------------------------------------------- //

aggte_gt, type("dynamic") bstrap("false")

display "--- Dynamic aggregation (event study) ---"
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 3. Group aggregation
//    Aggregates CATT(g,t,z) by treatment cohort g.
// -------------------------------------------------------------------------- //

aggte_gt, type("group") bstrap("false")

display "--- Group aggregation ---"
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 4. Calendar aggregation
//    Aggregates CATT(g,t,z) by calendar year t.
// -------------------------------------------------------------------------- //

aggte_gt, type("calendar") bstrap("false")

display "--- Calendar aggregation ---"
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 5. Simple aggregation
//    Weighted average across all (g,t) post-treatment pairs.
// -------------------------------------------------------------------------- //

aggte_gt, type("simple") bstrap("false")

display "--- Simple aggregation ---"
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")

display as text "Aggregation example complete."
