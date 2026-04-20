********************************************************************************
* example_visualization.do
*
* Visualization examples for the didhetero Stata package
* Implements Imai, Qin, and Yanagi (2025)
* "Doubly Robust Uniform Confidence Bands for Group-Time Conditional Average
* Treatment Effects in Difference-in-Differences"
*
* Empirical application: minimum wage increases and teen employment
* Dataset: min_wage_cs.dta — county-level panel, 2001–2007
*          Callaway & Sant'Anna (2021) / Imai, Qin & Yanagi (2025)
*
* This file demonstrates:
*   1. CATT mode plots (group-time level estimates by poverty rate)
*   2. Aggregated mode plots (dynamic event-study aggregation)
*   3. Saving graphs to file (save_path option)
*   4. Custom graph options (titles, axis labels, schemes)
*   5. Pre-trends visualization with zero reference line
*
* Graph elements:
*   - Shaded area (gs12%60): confidence bands
*     (bootstrap uniform CB if available, else analytical)
*   - Black thick line: point estimates
*   - Red dashed line at y=0: appears for pre-trend results
********************************************************************************

clear all
set more off

// ============================================================================
// Load and prepare the minimum wage dataset
// ============================================================================
// min_wage_cs.dta: county-level panel, 2001–2007
// lemp = log teen employment (outcome), pov = poverty rate (z)

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

describe
summarize lemp G period Z

local catt_inference "bstrap(false) uniformall(false)"
local agg_inference  "bstrap(false) uniformall(false)"

// Shared estimation options (IQR of poverty rate; instantaneous effects only)
local zv  "zeval(0.105 0.136 0.181)"
local gtv "gteval(2004 2004 2006 2006 2007 2007)"
local est_opts "porder(2) kernel(gau) bwselect(IMSE1) control_group(notyettreated)"


// ============================================================================
// Section 1: CATT mode plots
// ============================================================================
// CATT plots display group-time level treatment effect estimates as a
// function of the poverty rate z. One panel is produced per (g,t) pair.
//
// The graph auto-detects CATT mode from the 10-column e(results) matrix.
//   - Shaded area: analytical confidence bands around the point estimates
//   - Black line:  CATT point estimates as a function of the poverty rate

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `est_opts' `catt_inference'

* Plot CATT results — mode is auto-detected from e(results)
catt_gt_graph

disp as txt "CATT plot displayed: employment effect heterogeneity by poverty rate"


// ============================================================================
// Section 2: Aggregated mode plots
// ============================================================================
// Aggregated plots show treatment effects summarized across cohorts,
// here as a dynamic (event-study) aggregation.
//
// aggte_gt produces an e(Estimate) matrix whose column count depends
// on the aggregation type (9 columns for dynamic/group/calendar, 8 for
// simple); catt_gt_graph auto-detects the matching Aggregated mode.
//   - Shaded area: confidence bands for the aggregated estimates
//   - Black line:  aggregated CATT at each event time as a function of z

aggte_gt, type("dynamic") eval(0 1 2) `agg_inference'

* plot_type("Aggregated") can also be specified explicitly
catt_gt_graph, plot_type("Aggregated")

disp as txt "Aggregated (dynamic) plot displayed: event-study estimates by poverty rate"


// ============================================================================
// Section 3: Save graph to file
// ============================================================================
// catt_gt_graph can write the current plot to a file via save_path().
// Paths ending in .gph use graph save; image formats (.png, .pdf, .eps,
// .svg) use graph export.

* Re-run catt_gt because aggte_gt overwrites e(cmd)
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `est_opts' `catt_inference'

* Create output directory (capture suppresses error if it already exists)
capture mkdir output

* Save the CATT plot to PNG
* Explicitly specify plot_type("CATT") because a prior aggte_gt call
* may have left e(Estimate) in memory.
catt_gt_graph, plot_type("CATT") save_path("output/minwage_catt.png")

disp as txt "Graph saved to: output/minwage_catt.png"


// ============================================================================
// Section 4: Custom graph options
// ============================================================================
// The graph_opt() option passes additional Stata twoway graph options,
// allowing full control over titles, subtitles, axis labels, and schemes.

catt_gt_graph, plot_type("CATT") ///
    graph_opt(title("Minimum Wage Effect by County Poverty Rate") ///
    subtitle("CATT estimates — Imai, Qin, and Yanagi (2025)") ///
    xtitle("County poverty rate (z)") ///
    scheme(s2color))

disp as txt "Custom-styled CATT plot displayed"
disp as txt "  graph_opt() accepts any valid Stata twoway options:"
disp as txt "    title(), subtitle(), xtitle(), ytitle(), scheme(), etc."


// ============================================================================
// Section 5: Pre-trends visualization
// ============================================================================
// When catt_gt is run with the pretrend flag, the graph command
// automatically adds a red dashed horizontal line at y=0.
// This helps assess whether employment trends were already heterogeneous
// across poverty levels before the minimum wage changes.
//
// Graph elements in pre-trend mode:
//   - Shaded area: confidence bands around pre-treatment estimates
//   - Black line:  pre-treatment CATT(g,t,z) for t < g
//   - Red dashed line at y=0: null hypothesis (parallel trends)

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    `zv' `gtv' `est_opts' `catt_inference' pretrend

* Zero reference line appears automatically when e(pretrend) == 1
catt_gt_graph

disp as txt "Pre-trends plot displayed with red dashed zero reference line"
disp as txt "  Confidence bands containing zero support conditional parallel trends"


// ============================================================================
// Summary
// ============================================================================
disp as txt ""
disp as txt "============================================================"
disp as txt " Visualization example complete"
disp as txt "============================================================"
disp as txt ""
disp as txt " CATT mode:       group-time estimates across covariate values"
disp as txt " Aggregated mode: event-study dynamic aggregation"
disp as txt " Save to file:    save_path() writes .gph and exports .png/.pdf/.eps"
disp as txt " Custom options:  graph_opt() for titles, schemes, etc."
disp as txt " Pre-trends:      automatic y=0 reference line"
disp as txt ""
