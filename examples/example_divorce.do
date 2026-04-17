// ============================================================================
// example_divorce.do
// Empirical application: no-fault divorce reform and female suicide rates
// Dataset: divorce_sw.dta — state-level panel, 1964–1996
//          Stevenson & Wolfers (2006, QJE)
//
// Research question:
//   Does the impact of unilateral divorce laws on female suicide rates vary
//   with log state population (z)?
// ============================================================================

clear all
set more off

// divorce_sw.dta is installed with the package (net install didhetero).
use divorce_sw.dta, clear

* Y: female suicide rate (outcome)
* G: year of unilateral divorce law adoption (0 = never adopted in sample)
* Z: log state population (continuous covariate)
summarize Y G period Z

tab G if period == 1964

local zpts "7.5 7.8 8.1"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(1969 1969 1969 1970 1969 1971 1971 1971 1971 1972 1971 1973 1973 1973 1973 1974 1973 1975) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

catt_gt_graph

aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")
