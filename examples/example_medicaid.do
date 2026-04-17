clear all
set more off

adopath ++ "/Users/cxy/Desktop/2026project/didhetero/ado"

use "/Users/cxy/Desktop/2026project/didhetero/data/medicaid_county.dta", clear

summarize Y G period Z

tab G if period==2013

local zpts "12 18 24"

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
