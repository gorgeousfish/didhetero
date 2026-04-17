// ============================================================================
// example_pretrend.do
// Pre-trends testing example for the didhetero Stata package
//
// Implements Imai, Qin, and Yanagi (2025)
// "Doubly Robust Uniform Confidence Bands for Group-Time Conditional Average
// Treatment Effects in Difference-in-Differences"
//
// Empirical application: minimum wage increases and teen employment
// Dataset: min_wage_cs.dta — county-level panel, 2001–2007
//          Callaway & Sant'Anna (2021) / Imai, Qin & Yanagi (2025)
//
// Demonstrates how to use the `pretrend` option to test the conditional
// parallel trends assumption using real pre-treatment years (2001–2003).
// When the `pretrend` flag is specified, catt_gt estimates the CATT for
// pre-treatment periods. Under the null hypothesis of parallel trends,
// these pre-treatment effects should be zero at every poverty-rate value z.
//
// Date:    2026-04-12
// Author:  didhetero development team
// ============================================================================

clear all
set more off

// ============================================================================
// Step 1: Load and prepare the minimum wage dataset
// ============================================================================
// Variables:
//   lemp        — log teen employment (outcome)
//   first_treat — first year of minimum wage increase (0 = never-treated;
//                 2004, 2006, or 2007 = treated)
//   year        — calendar year (2001–2007)
//   countyreal  — county FIPS identifier
//   pov         — county poverty rate (continuous covariate z)
//
// Pre-treatment structure:
//   Group g=2004: years 2001, 2002, 2003 are pre-treatment
//   Group g=2006: years 2001–2005 are pre-treatment
//   Group g=2007: years 2001–2006 are pre-treatment
//
// If counties' employment trends were already diverging before minimum wage
// changes, the DiD parallel trends assumption would be violated. The
// pre-trends test using years 2001–2003 provides a direct diagnostic.

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

describe
summarize lemp G period Z

local inference_opts "bstrap(false) uniformall(false)"


// ============================================================================
// Step 2: Run pre-trends test
// ============================================================================
// The `pretrend` flag changes which (g,t) pairs are estimated. Without it,
// catt_gt estimates CATT only for post-treatment periods (t >= g). With
// `pretrend`, the command also estimates CATT for pre-treatment periods,
// excluding only the long-difference baseline t = g - 1.
//
// Under the conditional parallel trends assumption, the CATT at pre-treatment
// periods should be zero for every value of z (the poverty rate). Non-zero
// pre-treatment estimates indicate that employment trends were already
// heterogeneous across poverty levels before minimum wage changes, which
// would cast doubt on the identifying assumption.
//
// We evaluate at the IQR of the poverty rate (p25=0.105, p50=0.136, p75=0.181)
// and restrict to instantaneous-effect pairs to keep runtime manageable.

catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    `inference_opts' ///
    control_group("notyettreated") ///
    pretrend

// Display the full results matrix
// Columns: g, t, z, estimate, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
// Rows with t < g are pre-treatment estimates; these should be near zero.
matrix list e(results), format(%9.4f)

// Confirm that the pre-trends flag was recorded
display as text "Pre-trend flag stored in e(pretrend): " e(pretrend)


// ============================================================================
// Step 3: Visualize pre-trends
// ============================================================================
// catt_gt_graph detects e(pretrend) == 1 and automatically adds a red dashed
// horizontal line at y = 0. This reference line makes it easy to visually
// assess whether the estimated pre-treatment effects are consistent with
// the null of zero — i.e. whether the confidence bands contain the zero line.

catt_gt_graph


// ============================================================================
// Step 4: Interpretation guide
// ============================================================================
//
// PURPOSE OF PRE-TRENDS TESTING
// -----------------------------
// The difference-in-differences estimator identifies causal effects under the
// assumption that treated and control groups would have followed parallel
// outcome paths in the absence of treatment. This is fundamentally untestable
// because we never observe the counterfactual. However, we CAN check whether
// the groups were trending similarly BEFORE treatment — a necessary (though
// not sufficient) condition for parallel trends to hold after treatment.
//
// The `pretrend` option estimates the CATT at pre-treatment periods. If the
// parallel trends assumption holds, these estimates should be statistically
// indistinguishable from zero at every evaluation point z.
//
//
// HOW TO READ THE GRAPH
// ---------------------
// The graph produced by catt_gt_graph shows the estimated CATT(g,t,z) as a
// function of z for each (g,t) pair, with pointwise and uniform confidence
// bands. The red dashed line at y = 0 is the null hypothesis.
//
//   - If the confidence bands contain the zero line across all z values,
//     there is no evidence against parallel trends for that (g,t) cell.
//
//   - If the confidence bands exclude zero at some z values, this suggests
//     differential pre-trends at those covariate values. The treatment
//     effect heterogeneity along z may be confounded by pre-existing
//     differences.
//
// Pay attention to BOTH the pointwise bands (inner) and the uniform bands
// (outer). The uniform bands account for multiple testing across z values
// and provide simultaneous coverage. A rejection based on uniform bands is
// stronger evidence against parallel trends than a pointwise rejection.
//
//
// STATISTICAL LIMITATIONS
// -----------------------
// Pre-trends testing has well-known limitations that practitioners should
// keep in mind:
//
//   1. Failure to reject is NOT the same as the assumption holding.
//      A non-significant pre-trend test may simply reflect low statistical
//      power — the test may fail to detect real violations because the
//      sample is too small or the pre-treatment periods are too few.
//
//   2. Power depends on the number of pre-treatment periods and sample size.
//      With only one or two pre-treatment periods, the test has limited
//      ability to detect violations. More pre-treatment data strengthens
//      the test.
//
//   3. Pre-trends in levels vs. trends.
//      Parallel trends is an assumption about TRENDS (changes over time),
//      not about LEVELS. Groups can have different outcome levels and still
//      satisfy parallel trends. The pre-trend test checks for differential
//      changes, which is the correct object.
//
//   4. Pre-treatment parallel trends do not guarantee post-treatment
//      parallel trends. Even if groups trended identically before treatment,
//      they might have diverged after treatment for reasons unrelated to the
//      treatment itself (e.g., a contemporaneous shock affecting one group).
//
//
// WHAT TO DO IF PRE-TRENDS ARE VIOLATED
// --------------------------------------
// If the pre-trends test reveals significant pre-treatment effects, consider
// the following strategies:
//
//   1. Alternative control groups.
//      Switch between "notyettreated" and "nevertreated" control groups.
//      The parallel trends assumption differs across these choices, and one
//      may be more plausible than the other.
//
//   2. Additional covariates.
//      Include more variables in xformula() to condition on observable
//      differences that may drive the pre-trends. This strengthens the
//      conditional parallel trends assumption.
//
//   3. Different evaluation points.
//      The violation may be localised to certain values of z. Examine
//      whether the pre-trend is concentrated at the tails of the z
//      distribution, where estimation is noisier, or whether it is
//      pervasive across the support.
//
//   4. Sensitivity analysis.
//      Vary the bandwidth (bwselect or manual bw) and polynomial order
//      (porder) to check whether the pre-trend finding is robust to
//      estimation choices or driven by a particular specification.
//
//   5. Honest inference.
//      If pre-trends cannot be eliminated, consider bounds-based approaches
//      that allow for some degree of non-parallel trends (e.g., Rambachan
//      & Roth, 2023). These methods provide valid inference under weaker
//      assumptions at the cost of wider confidence intervals.

display _newline(2)
display as text "=============================================="
display as text " Pre-trends testing example complete"
display as text "=============================================="
display as text ""
display as text "Key stored results:"
display as text "  e(results)  — full results matrix (g, t, z, est, se, CIs, bw)"
display as text "  e(pretrend) — pre-trend flag (1 = pre-trends test was run)"
