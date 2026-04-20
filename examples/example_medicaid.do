// ============================================================================
// example_medicaid.do
// Empirical application: ACA Medicaid expansion and working-age mortality
// Dataset: medicaid_county.dta
//   County-level balanced panel, 2009-2019 (11 periods, 2,697 counties).
//   1,069 counties expanded Medicaid eligibility in 2014, 172 in 2015,
//   93 in 2016, 140 in 2019; 1,223 counties never adopted through 2019.
//   Running example in the companion Stata Journal article.
//
// Research question:
//   Does the county-level mortality response to Medicaid expansion vary
//   with pre-expansion baseline mortality Z (2011-2013 mean of the
//   working-age crude mortality rate per 100,000)?
// ============================================================================

clear all
set more off

// medicaid_county.dta is installed with the package (net install didhetero).
use medicaid_county.dta, clear

* Y: working-age (20-64) crude mortality rate per 100,000
* G: year of Medicaid expansion (0 = never adopted through 2019)
* Z: 2011-2013 baseline mortality mean (pre-expansion effect modifier)
* Z_poverty: 2013 county poverty rate (%), retained for robustness
summarize Y G period Z

tab G if period == 2013

// Evaluation points: 25th, 50th, 75th percentiles of Z rounded.
local zpts "340 425 525"

// (g,t) grid: impact + 1 + 2 years for 2014, 2015, 2016 cohorts.
catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(2014 2014 2014 2015 2014 2016 2015 2015 2015 2016 2015 2017 2016 2016 2016 2017 2016 2018) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

catt_gt_graph

aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")
