// ============================================================================
// example_castle.do
// Empirical application: castle doctrine / stand-your-ground laws
// Dataset: castle_doctrine.dta — state-level panel, 2000–2010
//          Cheng & Hoekstra (2013, JHR)
//
// Research question:
//   Does the increase in homicides following castle doctrine adoption
//   differ across U.S. states by log state population (z)?
// ============================================================================

clear all
set more off

// castle_doctrine.dta is installed with the package (net install didhetero).
use castle_doctrine.dta, clear

* Y: log homicide rate (outcome)
* G: year of castle doctrine adoption (0 = never treated in sample window)
* Z: log state population (continuous covariate)
summarize Y G period Z

tab G if period == 2000

* Evaluate near Z's 25th/50th/75th percentiles (9.0, 10.3, 13.0).
* Points far in the tails produce numerically unstable estimates
* when cohort sizes are small (g=2008 has only 4 states).
local zpts "9 11 13"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(2006 2006 2006 2007 2007 2007 2007 2008 2007 2009 2008 2008 2008 2009) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

catt_gt_graph

aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")
