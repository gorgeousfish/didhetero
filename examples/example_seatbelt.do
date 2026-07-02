// This example requires: ssc install didhetero
// ============================================================================
// example_seatbelt.do — Mandatory seat belt laws and traffic fatalities
// Dataset: seatbelt.dta — state-level panel, 1983–1997
//          Cohen & Einav (2003, RESTAT)
//
// Research question:
//   Does the fatality-reducing effect of seat belt laws vary
//   with 1983 per capita state income (z)?
// ============================================================================

clear all
set more off

use "../data/seatbelt.dta", clear

summarize Y G period Z
tab G if period == 1983

local zpts "10000 12000 14500"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(1986 1986 1986 1987 1986 1988 1987 1987 1987 1988 1987 1989 1988 1988 1988 1989) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    nobstrap nouniformall ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)

catt_gt_graph

aggte_gt, type("dynamic") eval(0 1 2) bstrap("false")
matrix list e(Estimate), format(%9.4f)

catt_gt_graph, plot_type("Aggregated")

display as text "Seat belt example complete."
