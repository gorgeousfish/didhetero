********************************************************************************
* example_aggte.do
*
* Aggregation Parameters Example for the didhetero Stata Package
*
* Empirical application: minimum wage increases and teen employment
* Dataset: min_wage_cs.dta — county-level panel, 2001–2007
*          Callaway & Sant'Anna (2021) / Imai, Qin & Yanagi (2025)
*
* This script demonstrates the four aggregation types available in aggte_gt.
* All four aggregate the group-time CATT(g,t,z) — the effect of minimum wage
* increases on log teen employment conditional on the county poverty rate z —
* along different dimensions.
*
* Aggregation types:
*   1. Dynamic  — aggregates by event time e = t - g
*   2. Group    — aggregates by treatment cohort g (2004, 2006, 2007)
*   3. Calendar — aggregates by calendar year t
*   4. Simple   — weighted average across all (g,t) pairs
*
* Reference: Imai, Qin, and Yanagi (2025)
*            "Doubly Robust Uniform Confidence Bands for Group-Time Conditional
*            Average Treatment Effects in Difference-in-Differences"
********************************************************************************

clear all
set more off

// -------------------------------------------------------------------------- //
// 1. Load and prepare data
// -------------------------------------------------------------------------- //

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

// Evaluation points: IQR of the county poverty rate
// p25 = 0.105, p50 = 0.136, p75 = 0.181
local zv   "zeval(0.105 0.136 0.181)"
local gtv  "gteval(2004 2004 2006 2006 2007 2007)"
local opts "porder(2) kernel(gau) bwselect(IMSE1) control_group(notyettreated)"

// -------------------------------------------------------------------------- //
// 2. Dynamic aggregation (event study)
//
//    Aggregates CATT(g,t,z) by event time e = t - g.
//      e = 0 : year of initial minimum wage increase
//      e = 1 : one year after the increase
//      e = 2 : two years after the increase
//    This shows how the poverty-rate-heterogeneous employment effect
//    evolves over time since treatment onset.
// -------------------------------------------------------------------------- //

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `opts' bstrap biters(500) uniformall

aggte_gt, type("dynamic") eval(0 1 2) bstrap("true") biters(500)

display "--- Dynamic aggregation (event study) ---"
matrix list e(Estimate), format(%9.4f)
display "Aggregation type: " e(type)

catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 3. Group aggregation
//
//    Aggregates CATT(g,t,z) by treatment cohort g.
//      g = 2004 : counties whose minimum wage first increased in 2004
//      g = 2006 : counties whose minimum wage first increased in 2006
//      g = 2007 : counties whose minimum wage first increased in 2007
//    This reveals whether earlier or later adopters experience different
//    poverty-rate-heterogeneous employment effects.
// -------------------------------------------------------------------------- //

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `opts' bstrap biters(500) uniformall

aggte_gt, type("group") eval(2004 2006 2007) bstrap("true") biters(500)

display "--- Group aggregation ---"
matrix list e(Estimate), format(%9.4f)
display "Aggregation type: " e(type)

catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 4. Calendar aggregation
//
//    Aggregates CATT(g,t,z) by calendar year t.
//      t = 2004 : average employment effect across treated cohorts in 2004
//      t = 2006 : average employment effect across treated cohorts in 2006
//      t = 2007 : average employment effect across treated cohorts in 2007
//    This captures aggregate time-varying impacts of minimum wage policy.
// -------------------------------------------------------------------------- //

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `opts' bstrap biters(500) uniformall

aggte_gt, type("calendar") eval(2004 2006 2007) bstrap("true") biters(500)

display "--- Calendar aggregation ---"
matrix list e(Estimate), format(%9.4f)
display "Aggregation type: " e(type)

catt_gt_graph, plot_type("Aggregated")

// -------------------------------------------------------------------------- //
// 5. Simple aggregation
//
//    Computes a single weighted-average employment effect curve as a function
//    of the poverty rate z, averaging across all (g,t) post-treatment pairs.
//    This is the most concise overall summary of how minimum wage effects
//    vary with the poverty rate of a county.
// -------------------------------------------------------------------------- //

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `opts' bstrap biters(500) uniformall

aggte_gt, type("simple") bstrap("true") biters(500)

display "--- Simple aggregation ---"
matrix list e(Estimate), format(%9.4f)
display "Aggregation type: " e(type)

catt_gt_graph, plot_type("Aggregated")
