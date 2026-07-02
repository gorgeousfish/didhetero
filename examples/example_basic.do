// This example requires: ssc install didhetero
// ============================================================================
// example_basic.do — Core workflow for the didhetero Stata package
//
// Minimum wage increases and teen employment (Callaway & Sant'Anna 2021)
// Research question: Does the employment effect vary with county poverty (z)?
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare data
// ============================================================================

use "../data/min_wage_cs.dta", clear
rename (first_treat year countyreal) (G period id)

* Time-invariant Z from first-period poverty rate
bysort id (period): gen Z = pov[1]

summarize lemp G period Z
tab G

// ============================================================================
// Step 2: Estimate CATT — poverty rate as the effect modifier
// ============================================================================
// Evaluate at IQR of poverty rate: p25=0.105, p50=0.136, p75=0.181
// nobstrap/nouniformall for fast analytical inference

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2004 2005 2004 2006 2006 2006 2006 2007 2007 2007) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

// Results: g, t, z, estimate, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
matrix list e(results), format(%9.4f)

// ============================================================================
// Step 3: Visualize CATT estimates
// ============================================================================
// Each panel: one (g,t) pair; x-axis = poverty rate

catt_gt_graph

// ============================================================================
// Step 4: Simple aggregation — overall effect curve
// ============================================================================
// Weighted average across all (g,t) post-treatment pairs

aggte_gt, type("simple") bstrap("false")
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")

display as text "Basic example complete."
