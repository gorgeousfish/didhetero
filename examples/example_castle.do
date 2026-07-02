// This example requires: ssc install didhetero
// ============================================================================
// example_castle.do — Castle doctrine / stand-your-ground laws
// Dataset: castle_doctrine.dta — state-level panel, 2000–2010
//          Cheng & Hoekstra (2013, JHR)
//
// Research question:
//   Does the homicide increase following castle doctrine adoption
//   differ across states by log population (z)?
// ============================================================================

clear all
set more off

use "../data/castle_doctrine.dta", clear

summarize Y G period Z
tab G if period == 2000

* Evaluate near Z's 25th/50th/75th percentiles
local zpts "9 11 13"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(2006 2006 2006 2007 2007 2007 2007 2008 2007 2009 2008 2008 2008 2009) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

catt_gt_graph

aggte_gt, type("dynamic") eval(0 1 2) bstrap("false")
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")

display as text "Castle doctrine example complete."
