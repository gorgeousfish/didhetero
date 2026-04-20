# didhetero

**Heterogeneous Treatment Effects in Staggered Difference-in-Differences for Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)]()

![didhetero](image/image.png)

## Overview

`didhetero` implements the **doubly robust estimator** for group-time conditional average treatment effects on the treated (CATT) proposed by Imai, Qin, and Yanagi (2025, *Journal of Business & Economic Statistics*). It extends the standard staggered difference-in-differences framework of Callaway and Sant'Anna (2021) to allow treatment effect heterogeneity driven by a **continuous pre-treatment covariate**.

The estimator combines inverse probability weighting with outcome regression using local polynomial methods, producing nonparametric estimates of how treatment effects vary as a function of a continuous covariate *z*. Uniform confidence bands are constructed via analytical distributional approximation and multiplier bootstrap procedures.

**Features**:

- **Doubly robust estimator** combining inverse probability weighting with outcome regression via local polynomial smoothing
- **Uniform confidence bands** for CATT functions via analytical or multiplier bootstrap methods
- **Flexible aggregation** into event-study (dynamic), group, calendar, and simple summary parameters
- **Built-in visualization** for estimated CATT curves and aggregated parameters
- **Pre-trends testing** for the conditional parallel trends assumption
- **Bundled datasets** including the ACA Medicaid expansion analyzed as the running example in the companion Stata Journal article, the minimum wage panel of Callaway and Sant'Anna (2021) used as the running example in Imai, Qin, and Yanagi (2025), and three additional empirical applications (mandatory seat belt laws, no-fault divorce reform, castle doctrine)

## Key Concepts

### Treatment Effect Heterogeneity in DiD

Standard DiD methods estimate average treatment effects that are homogeneous across covariate values. In many applications, however, treatment effects may vary systematically with pre-treatment characteristics. For example, the effect of minimum wage increases on teen employment may differ across counties with different poverty rates.

The **group-time conditional average treatment effect on the treated (CATT)** captures this heterogeneity:

```
CATT(g, t, z) = E[Y_t(g) - Y_t(0) | Z = z, G = g]
```

where *g* is the treatment group (first treatment period), *t* is the calendar period, and *z* is a continuous pre-treatment covariate.

### Three-Step Estimation Procedure

The estimation follows a three-step procedure:

| Stage | Description | Method |
|-------|-------------|--------|
| **Stage 1** | Estimate the generalized propensity score (GPS) and outcome regression (OR) function | Parametric (logit/OLS) |
| **Stage 2** | Construct intermediate doubly robust influence function variables | Local polynomial regression for nuisance parameters |
| **Stage 3** | Estimate the conditional DR estimand to obtain CATT | Local polynomial regression |

The doubly robust property ensures that CATT is consistently estimated if **either** the GPS or the OR function is correctly specified, but not necessarily both.

### Uniform Inference

Two methods for constructing uniform confidence bands are available:

- **Analytical**: Based on distributional approximation results for suprema of empirical processes
- **Multiplier bootstrap**: Based on weighted/multiplier bootstrap using Mammen's (1993) weights

The uniform confidence bands provide simultaneous coverage over all evaluation points *z* and optionally over all (g,t) pairs, controlling for multiple testing.

### Aggregation Parameters

After estimating CATT(g,t,z), the results can be aggregated into interpretable summary parameters:

| Type | Description | Aggregation Dimension |
|------|-------------|-----------------------|
| **Dynamic** | Event-study aggregation by relative time since treatment onset | Event time *e* |
| **Group** | Aggregation by treatment cohort | Group *g* |
| **Calendar** | Aggregation by calendar period | Period *t* |
| **Simple** | Weighted average across all post-treatment (g,t) pairs | Single summary curve |

## Requirements

- Stata 16.0 or later
- No additional dependencies (Mata included in Stata)

## Installation

### From GitHub

```stata
net install didhetero, from("https://raw.githubusercontent.com/gorgeousfish/didhetero/main/") replace
```

This installs the commands, help files, and compiled Mata library (`ldidhetero.mlib`) into your personal Stata ado directory.

### From a local directory

```stata
net install didhetero, from("/path/to/didhetero-stata/") replace
```

### Download example scripts and bundled datasets

Bundled datasets (`medicaid_county.dta`, `min_wage_cs.dta`, `seatbelt.dta`, `castle_doctrine.dta`, `divorce_sw.dta`, `simulated_data.dta`) and example do-files are distributed as ancillary files. Fetch them into your current working directory with:

```stata
net get didhetero
```

Run `net get` from a dedicated analysis folder; the files land in the current working directory and are available via plain `use filename.dta, clear` afterward.

### Verify Installation

```stata
which didhetero
help didhetero
```

## Quick Start

```stata
* Load the bundled minimum wage dataset
use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Estimate CATT at three poverty-rate values for instantaneous effects
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.10 0.14 0.18) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false)

* Visualize CATT curves as a function of the poverty rate
catt_gt_graph

* Aggregate into event-study parameters
aggte_gt, type(dynamic) bstrap(false) uniformall(false)
catt_gt_graph
```

## Commands

| Command | Description |
|---------|-------------|
| `didhetero` | Unified entry point for CATT estimation |
| `catt_gt` | Estimate group-time CATT functions CATT(g,t,z) |
| `aggte_gt` | Aggregate CATT into summary parameters (dynamic, group, calendar, simple) |
| `catt_gt_graph` | Visualize CATT or aggregated estimates |
| `didhetero_simdata` | Generate simulation data for testing |

## Options

### didhetero / catt_gt Options

`didhetero` and `catt_gt` share the same estimation engine and options.

| Option | Description | Default |
|--------|-------------|---------|
| `id(varname)` | Panel individual identifier | Required |
| `time(varname)` | Time period variable | Required |
| `group(varname)` | Treatment group variable; 0 = never-treated | Required |
| `z(varname)` | Continuous pre-treatment covariate | Required |
| `zeval(numlist)` | Evaluation points for *z* | Required |
| `xformula(string)` | Covariate specification for outcome regression | X = [1, z] |
| `gteval(numlist)` | Specific (g,t) pairs: `g1 t1 g2 t2 ...` | All valid pairs |
| `porder(#)` | Local polynomial order; 1 or 2 | 2 |
| `kernel(string)` | Kernel function; `gau` or `epa` | gau |
| `bwselect(string)` | Bandwidth selection; `IMSE1`, `IMSE2`, `US1`, or `manual` | IMSE1 |
| `bw(numlist)` | Manual bandwidth (requires `bwselect(manual)`) | - |
| `alp(#)` | Significance level (0, 1) | 0.05 |
| `bstrap(true\|false)` | Bootstrap toggle | true |
| `biters(#)` | Bootstrap iterations | 1000 |
| `seed(#)` | Random-number seed | - |
| `uniformall(true\|false)` | Joint uniform bands over all (g,t,z) | true |
| `control_group(string)` | `notyettreated` or `nevertreated` | notyettreated |
| `anticipation(#)` | Anticipation periods | 0 |
| `pretrend` | Include pre-treatment periods for testing | - |

**Covariate specification (`xformula`)**:

When `xformula()` is omitted, the design matrix defaults to X = [1, z]. Supported syntax includes bare variables, `I(expr)`, `:` (interaction only), `*` (main effects plus interaction), and no-intercept tokens `0 +` / `- 1`.

```stata
* Default: X = [1, Z]
catt_gt Y, ... z(Z) zeval(...)

* Explicit: same as default
catt_gt Y, ... z(Z) zeval(...) xformula(Z)

* Add covariates
catt_gt Y, ... z(Z) zeval(...) xformula("Z + X1")

* Interaction terms
catt_gt Y, ... z(Z) zeval(...) xformula("Z*X1")
```

### aggte_gt Options

`aggte_gt` is a post-estimation command. It must be run after `catt_gt` or `didhetero`.

| Option | Description | Default |
|--------|-------------|---------|
| `type(string)` | Aggregation type; `dynamic`, `group`, `calendar`, or `simple` | dynamic |
| `eval(numlist)` | Evaluation points for the aggregation dimension | All available |
| `porder(#)` | Local polynomial order; 1 or 2 | 2 |
| `bwselect(string)` | Bandwidth selection method | IMSE1 |
| `bw(numlist)` | Manual bandwidth override | - |
| `bstrap(string)` | `true` or `false`; enable bootstrap | true |
| `biters(#)` | Bootstrap iterations | 1000 |
| `seed(#)` | Random-number seed for bootstrap | - |
| `uniformall(string)` | `true` or `false`; joint uniform bands | true |

### catt_gt_graph Options

`catt_gt_graph` is a post-estimation command. It must be run after `didhetero`, `catt_gt`, or `aggte_gt`.

| Option | Description | Default |
|--------|-------------|---------|
| `plot_type(string)` | Force plot mode; `CATT` or `Aggregated` | Auto-detect |
| `save_path(string)` | File path to save the graph (.png, .pdf, .gph) | - |
| `graph_opt(string)` | Additional Stata `twoway` graph options | - |

### didhetero_simdata Options

| Option | Description | Default |
|--------|-------------|---------|
| `n(#)` | Number of cross-sectional units | Required |
| `tau(#)` | Number of time periods; must be > 2 | Required |
| `seed(#)` | Random seed for reproducibility | - |
| `hc` | Use heteroscedastic error terms | Off |
| `dimx(#)` | Dimension of covariates X | 1 |
| `dgpy(#)` | DGP for treated potential outcome; 1 or 2 | 1 |
| `discrete` | Use discrete Z in {-1, 0, 1} | Off |
| `clear` | Allow replacing data in memory | - |

## Recommended Workflow

### Step 1: Load Data

```stata
* Load the bundled minimum wage dataset (installed with the package)
use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Required variables in your own data:
*   Y      — outcome variable
*   id     — panel individual identifier
*   period — time period variable
*   G      — treatment group (0 = never-treated, g = first treated in period g)
*   Z      — continuous pre-treatment covariate (time-invariant)
```

### Step 2: Test Pre-trends (Optional but Recommended)

```stata
* Estimate CATT including pre-treatment periods (2001–2003 serve as pre-treatment)
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.10 0.14 0.18) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false) pretrend

* Visualize — zero reference line appears automatically
catt_gt_graph
```

Under the parallel trends assumption, pre-treatment CATT estimates should be near zero for all *z*.

### Step 3: Estimate CATT

```stata
* Estimate CATT for instantaneous effects (g=2004, 2006, 2007 at t=g)
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.10 0.14 0.18) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false)

* View results
matrix list e(results)

* Visualize CATT curves — how the minimum wage effect varies with poverty rate
catt_gt_graph
```

### Step 4: Aggregate Treatment Effects

```stata
* Event-study aggregation
aggte_gt, type(dynamic) bstrap(false) uniformall(false)
catt_gt_graph

* Group aggregation (re-run catt_gt first since aggte_gt overwrites e())
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.10 0.14 0.18) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false)
aggte_gt, type(group) bstrap(false) uniformall(false)
catt_gt_graph

* Simple weighted average
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.10 0.14 0.18) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false)
aggte_gt, type(simple) bstrap(false) uniformall(false)
catt_gt_graph
```

## Datasets

### ACA Medicaid Expansion (`medicaid_county.dta`)

County-level balanced panel on working-age (20–64) crude mortality rates, 2009–2019 (11 periods, 2,697 U.S. counties; 29,667 observations). States expanded Medicaid eligibility under the ACA in a staggered manner from 2014 onward: 1,069 counties in 2014, 172 in 2015, 93 in 2016, 140 in 2019, and 1,223 never adopted through 2019. The **running example** in the companion Stata Journal article.

**Data structure (`medicaid_county.dta`):**

| Variable | Description |
|----------|-------------|
| `Y` | Working-age (20–64) crude mortality rate per 100,000 |
| `G` | Year of Medicaid expansion (0 = never adopted through 2019) |
| `period` | Calendar year (2009–2019) |
| `id` | County FIPS identifier |
| `Z` | 2011–2013 mean of `Y`: pre-expansion baseline mortality per 100,000 |
| `Z_poverty` | 2013 county poverty rate (%), retained for robustness checks |
| `state`, `stfips` | State name and FIPS code |

**Research question**: Does the county-level mortality response to Medicaid expansion vary with pre-expansion baseline mortality?

**Data preparation**: See the replication materials (`replication/code/00_prep_baseline_mortality.do`) for the construction of `Z` from the raw 2009–2019 county-level mortality series.

### Minimum Wage — Callaway and Sant'Anna (2021)

The bundled `min_wage_cs.dta` dataset reproduces the empirical application in Imai, Qin, and Yanagi (2025). The data contain county-level panel observations on teen employment in the United States.

**Data structure (`min_wage_cs.dta`):**

| Variable | Description |
|----------|-------------|
| `lemp` | Log teen employment (outcome *Y*) |
| `first_treat` | First minimum wage increase year (0 = never-treated; 2004, 2006, or 2007 = treated) |
| `year` | Calendar year (2001–2007) |
| `countyreal` | County FIPS identifier |
| `pov` | County poverty rate — continuous covariate *z* |
| `lmedinc` | Log median county income — alternative covariate |
| `lpop` | Log county population |
| `col` | College attendance rate |
| `white`, `black` | Racial composition shares |

**Research question**: Does the employment effect of minimum wage increases vary systematically with the pre-treatment poverty rate of a county?

**Note on computation**: The full dataset (2,284 counties) requires substantial local polynomial estimation. All examples below use `gteval()` to evaluate only selected (g,t) pairs to keep runtimes manageable.

```stata
* Standard data loading and variable renaming used in all examples below
use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)
```

### Mandatory Seat Belt Laws (`seatbelt.dta`)

State-level panel on traffic fatalities (1983–1997) from the `USSeatBelts` dataset in the AER R package. U.S. states adopted mandatory seat belt enforcement laws in a staggered manner from 1985 to 1989.

**Data structure (`seatbelt.dta`):**

| Variable | Description |
|----------|-------------|
| `Y` | Traffic fatality rate (fatalities per million vehicle miles driven) |
| `G` | First year of mandatory seat belt enforcement law (0 = late/never adopter used as control) |
| `period` | Calendar year (1983–1997) |
| `id` | State identifier |
| `Z` | 1983 per capita income (USD, continuous covariate) |

**Research question**: Does the fatality-reducing effect of mandatory seat belt laws vary with pre-treatment state income? (Hypothesis: lower-income states had lower baseline compliance, so the mandate has a larger effect.)

**Data preparation**: Run `stata-mcp-folder/_build_seatbelt_data2.R` to construct from the `AER` R package (naturally balanced, no observations dropped).

### No-Fault Divorce Reform (`divorce_sw.dta`)

State-level panel on female suicide rates (1964–1996) from Stevenson and Wolfers (2006). U.S. states adopted unilateral (no-fault) divorce laws in a staggered manner from the late 1960s onward.

**Data structure (`divorce_sw.dta`):**

| Variable | Description |
|----------|-------------|
| `Y` | Female suicide rate (outcome) |
| `G` | Year of unilateral divorce law adoption (0 = never adopted in sample) |
| `period` | Calendar year (1964–1996) |
| `id` | State identifier |
| `Z` | Log state population (continuous covariate) |

**Research question**: Does the impact of divorce law liberalization on female suicide vary with state population size?

**Data preparation**: Run `stata-mcp-folder/_build_divorce_data.do` to construct this dataset from the `bacondecomp` R package.

### Castle Doctrine Laws (`castle_doctrine.dta`)

State-level panel on homicide rates (2000–2010) from Cheng and Hoekstra (2013), studying the staggered rollout of "stand your ground" / castle doctrine laws across U.S. states.

**Data structure (`castle_doctrine.dta`):**

| Variable | Description |
|----------|-------------|
| `Y` | Log homicide rate (outcome) |
| `G` | Year of castle doctrine adoption (0 = never treated in sample window) |
| `period` | Calendar year (2000–2010) |
| `id` | State identifier |
| `Z` | Log state population (continuous covariate) |

**Research question**: Does the increase in homicides following castle doctrine adoption differ across states by population size?

**Data preparation**: Run `stata-mcp-folder/_build_empirical_data.do` to download and construct this dataset.

## Examples

### Example 1: Basic Workflow — CATT Estimation

This example estimates the group-time CATT as a function of the county poverty rate, reproducing the core empirical analysis of Imai, Qin, and Yanagi (2025).

```stata
clear all
set more off

* Load and prepare data
use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Examine treatment groups and time periods
tab G
tab period

* Estimate CATT at the 25th, 50th, and 75th percentiles of the poverty rate
* (p25 = 0.105, p50 = 0.136, p75 = 0.181)
* Restrict to instantaneous effects: (g=2004,t=2004), (g=2006,t=2006), (g=2007,t=2007)
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2006 2006 2007 2007) ///
    bstrap(false) uniformall(false)

* View the results matrix
* Columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
matrix list e(results), format(%9.4f)

* Visualize CATT curves — one panel per (g,t) pair
catt_gt_graph
```

### Example 2: Event-Study Aggregation

This example aggregates CATT(g,t,z) into event-study parameters that show how the poverty-rate-heterogeneous effect evolves from the period of treatment onset onward.

```stata
clear all
set more off

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Step 1: Estimate CATT for all required (g,t) pairs
* aggte_gt reads from stored catt_gt results, so gteval must include every
* (g,t) pair needed for the requested eval points:
*   e=0: (2004,2004), (2006,2006), (2007,2007)
*   e=1: (2004,2005), (2006,2007)
*   e=2: (2004,2006)   [only G=2004 has a third post-period in this panel]
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2004 2005 2004 2006 2004 2007 ///
           2006 2006 2006 2007 ///
           2007 2007) ///
    bstrap(false) uniformall(false)

* Step 2: Dynamic (event-study) aggregation
* e=0: period of initial minimum wage increase
* e=1: one year after the increase
* e=2: two years after the increase
aggte_gt, type(dynamic) eval(0 1 2) bstrap(false) uniformall(false)

* View aggregated results
* Columns: eval, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
matrix list e(Estimate), format(%9.4f)

* Plot: CATT at each event time as a function of the poverty rate
catt_gt_graph, plot_type("Aggregated")
```

### Example 3: All Four Aggregation Types

This example demonstrates all four aggregation types on the minimum wage data, revealing how poverty-rate heterogeneity manifests across event time, treatment cohort, calendar time, and a simple overall summary.

```stata
clear all
set more off

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Estimate CATT once with all (g,t) pairs needed by all four aggregation types.
* aggte_gt can be called multiple times after one catt_gt — it reads from e().
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2004 2005 2004 2006 2004 2007 ///
           2006 2006 2006 2007 ///
           2007 2007) ///
    bstrap(false) uniformall(false)

* --- Dynamic aggregation (event study) ---
* Aggregates by relative time e = t - g across all cohorts
aggte_gt, type(dynamic) eval(0 1 2) bstrap(false) uniformall(false)
catt_gt_graph, plot_type("Aggregated")

* --- Group aggregation ---
* Average effect for each treatment cohort (g = 2004, 2006, 2007)
aggte_gt, type(group) eval(2004 2006 2007) bstrap(false) uniformall(false)
catt_gt_graph, plot_type("Aggregated")

* --- Calendar aggregation ---
* Average effect at each calendar year across all treated cohorts
aggte_gt, type(calendar) eval(2004 2006 2007) bstrap(false) uniformall(false)
catt_gt_graph, plot_type("Aggregated")

* --- Simple aggregation ---
* Single weighted-average effect curve as a function of the poverty rate
aggte_gt, type(simple) bstrap(false) uniformall(false)
catt_gt_graph, plot_type("Aggregated")
```

### Example 4: Pre-trends Testing

This example tests the conditional parallel trends assumption using the pre-treatment years 2001–2003. For groups first treated in 2004, 2006, or 2007, all years before their treatment onset are pre-treatment periods.

```stata
clear all
set more off

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* The pretrend option automatically includes all pre-treatment periods t < g
* for each cohort. It cannot be combined with gteval(); omit gteval() here
* so catt_gt evaluates all valid (g,t) pairs including pre-treatment ones.
* Note: estimating all pairs may take a few minutes on this dataset.
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    bstrap(false) uniformall(false) pretrend

* View results — pre-treatment estimates should be near zero
* if the conditional parallel trends assumption holds
matrix list e(results), format(%9.4f)

* Plot with automatic red dashed zero-reference line
* Confidence bands containing zero support the parallel trends assumption
catt_gt_graph
```

**Interpreting pre-trends**: If counties' employment paths were already diverging before minimum wage changes, that would invalidate the DiD identifying assumption. The pre-treatment CATT estimates — and their confidence bands relative to the zero line — provide a diagnostic check.

### Example 5: Bootstrap Uniform Confidence Bands

This example produces multiplier-bootstrap uniform confidence bands (the `ci2` columns), which provide simultaneous coverage over all evaluation points *z* and all (g,t) pairs — the preferred inference method for the analysis in Imai, Qin, and Yanagi (2025).

```stata
clear all
set more off

use min_wage_cs.dta, clear
rename (first_treat year countyreal pov) (G period id Z)

* Estimate with multiplier bootstrap (500 iterations for speed; use 1000 for publication)
* Include all (g,t) pairs required for dynamic aggregation at e=0,1,2
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(0.105 0.136 0.181) ///
    gteval(2004 2004 2004 2005 2004 2006 2004 2007 ///
           2006 2006 2006 2007 ///
           2007 2007) ///
    seed(42) biters(500)

* ci2 columns now contain bootstrap uniform confidence bands
matrix list e(results), format(%9.4f)
catt_gt_graph

* Event-study aggregation with bootstrap uniform bands
aggte_gt, type(dynamic) eval(0 1 2) seed(42) biters(500)
catt_gt_graph, plot_type("Aggregated")
```

### Example 6: Alternative Covariate — Log Median Income

This example substitutes the county log median income (`lmedinc`) for the poverty rate as the continuous effect modifier *z*, and adds population size as an additional control in `xformula`. It illustrates how to use covariates beyond the default specification.

```stata
clear all
set more off

use min_wage_cs.dta, clear

* Use log median income as the heterogeneity-driving covariate z
* Evaluate at the 25th, 50th, and 75th percentiles of lmedinc
* (p25 ≈ 3.33, p50 ≈ 3.46, p75 ≈ 3.59)
rename (first_treat year countyreal lmedinc) (G period id Z)

* Summarize the covariate distribution to choose evaluation points
summarize Z, detail

* Estimate CATT: does the employment effect of minimum wages vary
* more in low-income or high-income counties?
catt_gt lemp, group(G) time(period) id(id) z(Z) ///
    zeval(3.33 3.46 3.59) ///
    gteval(2004 2004 2004 2005 2004 2006 2004 2007 ///
           2006 2006 2006 2007 ///
           2007 2007) ///
    xformula("Z + lpop") ///
    bstrap(false) uniformall(false)

matrix list e(results), format(%9.4f)
catt_gt_graph

* Dynamic aggregation
aggte_gt, type(dynamic) eval(0 1 2) bstrap(false) uniformall(false)
catt_gt_graph, plot_type("Aggregated")
```

### Example 7: Empirical Application — Mandatory Seat Belt Laws

This example applies `didhetero` to the Cohen & Einav (2003) seat belt data to examine whether the fatality-reducing effect of mandatory enforcement laws varies with state income. States with lower pre-treatment income had lower baseline seat belt usage, so mandatory laws may have a larger effect there. Requires `seatbelt.dta` (naturally balanced, no observations dropped; see [data preparation](#mandatory-seat-belt-laws-seatbeltdta)).

```stata
clear all
set more off

use seatbelt.dta, clear

* Y: fatality rate per million miles, Z: 1983 per capita income (USD)
* G=0: states adopting in 1990 (used as control group)
summarize Y G period Z
tab G if period == 1983

* Evaluate at three income levels: low ($10k), median ($12k), high ($14.5k)
local zpts "10000 12000 14500"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(1986 1986 1986 1987 1986 1988 1987 1987 1987 1988 1987 1989 1988 1988 1988 1989) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)
catt_gt_graph

* Event-study aggregation
aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")
```

### Example 8: Empirical Application — No-Fault Divorce Reform

This example applies `didhetero` to Stevenson and Wolfers' (2006) divorce-reform panel to examine whether the impact of unilateral divorce laws on female suicide varies with state population. Requires `divorce_sw.dta` (see [data preparation](#no-fault-divorce-reform-divorce_swdta)).

```stata
clear all
set more off

use divorce_sw.dta, clear

* Y: female suicide rate, Z: log state population
summarize Y G period Z
tab G if period == 1964

* Evaluate at three log-population levels (7.5, 7.8, 8.1)
local zpts "7.5 7.8 8.1"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(1969 1969 1969 1970 1969 1971 1971 1971 1971 1972 1971 1973 1973 1973 1973 1974 1973 1975) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

* Note: Y is female suicide rate (~2e-4); use high-precision format to avoid 0.0000 display
matrix list e(results), format(%12.8f)
catt_gt_graph

* Event-study aggregation
aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%12.8f)
catt_gt_graph, plot_type("Aggregated")
```

### Example 9: Empirical Application — Castle Doctrine Laws

This example applies `didhetero` to Cheng and Hoekstra's (2013) data to test whether the increase in homicides following "stand your ground" law adoption differs across U.S. states by population size. Requires `castle_doctrine.dta` (see [data preparation](#castle-doctrine-laws-castle_doctrinedta)).

```stata
clear all
set more off

use castle_doctrine.dta, clear

* Y: log homicide rate, Z: log state population
summarize Y G period Z
tab G if period == 2000

* Evaluate near the 25th/50th/75th percentiles of log state population
* (p25 = 9.02, p50 = 10.31, p75 = 13.02); values outside the 5th–95th
* percentile range are ill-advised given the small n=50 and tiny cohorts.
local zpts "9 11 13"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(2006 2006 2006 2007 2007 2007 2007 2008 2007 2009 2008 2008 2008 2009) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)
catt_gt_graph

* Event-study aggregation
aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")
```

### Example 10: Empirical Application — ACA Medicaid Expansion

This example applies `didhetero` to the ACA Medicaid expansion panel used as the running example in the companion Stata Journal article. The research question is whether the county-level mortality response to Medicaid expansion varies with pre-expansion baseline mortality. Requires `medicaid_county.dta` (see [data preparation](#aca-medicaid-expansion-medicaid_countydta)).

```stata
clear all
set more off

use medicaid_county.dta, clear

* Y: working-age mortality rate, Z: 2011-2013 baseline mortality
summarize Y G period Z
tab G if period == 2013

* Evaluation points: 25th, 50th, 75th percentiles of Z rounded
local zpts "340 425 525"

catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(`zpts') ///
    gteval(2014 2014 2014 2015 2014 2016 2015 2015 2015 2016 2015 2017 2016 2016 2016 2017 2016 2018) ///
    porder(2) kernel("gau") bwselect("IMSE1") ///
    bstrap(false) uniformall(false) ///
    control_group("notyettreated")

matrix list e(results), format(%9.4f)
catt_gt_graph

* Event-study aggregation
aggte_gt, type("dynamic") eval(0 1 2) bstrap(false) uniformall(false)
matrix list e(Estimate), format(%9.4f)
catt_gt_graph, plot_type("Aggregated")
```

## Stored Results

### catt_gt / didhetero

**Scalars:**

| Result | Description |
|--------|-------------|
| `e(N)` | Number of cross-sectional units |
| `e(T)` | Number of time periods |
| `e(gbar)` | Upper bound for treatment-time support; `.` encodes +Inf |
| `e(num_gteval)` | Number of (g,t) pairs evaluated |
| `e(num_zeval)` | Number of z evaluation points |
| `e(porder)` | Polynomial order |
| `e(anticipation)` | Anticipation periods |
| `e(alp)` | Significance level |
| `e(bstrap)` | 1 if bootstrap was performed |
| `e(biters)` | Number of bootstrap iterations |
| `e(uniformall)` | 1 if joint uniform bands requested |
| `e(pretrend)` | 1 if pre-trends testing was requested |

**Macros:**

| Result | Description |
|--------|-------------|
| `e(cmd)` | Command name |
| `e(depvar)` | Outcome variable name |
| `e(idvar)` | Panel id variable |
| `e(timevar)` | Time variable |
| `e(groupvar)` | Group variable |
| `e(zvar)` | Covariate variable |
| `e(kernel)` | Kernel function used |
| `e(control_group)` | Control group specification |
| `e(bwselect)` | Bandwidth selection method |

**Matrices:**

| Result | Description |
|--------|-------------|
| `e(results)` | Main results (g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw) |
| `e(Estimate)` | Alias of `e(results)` |
| `e(gteval)` | Evaluated (g,t) pairs (K x 2) |
| `e(zeval)` | Evaluation points for z |
| `e(bw)` | Bandwidth per (g,t) pair |
| `e(c_hat)` | Analytical critical values |
| `e(c_check)` | Bootstrap critical values (if bootstrap) |
| `e(catt_est)` | CATT point estimates |
| `e(catt_se)` | CATT standard errors |

### aggte_gt

**Scalars:**

| Result | Description |
|--------|-------------|
| `e(porder)` | Aggregation polynomial order |
| `e(bstrap)` | 1 if bootstrap enabled |
| `e(uniformall)` | 1 if joint uniform bands effective |
| `e(alp)` | Significance level |
| `e(biters)` | Effective bootstrap iterations |

**Macros:**

| Result | Description |
|--------|-------------|
| `e(cmd)` | `aggte_gt` |
| `e(type)` | Aggregation type |
| `e(depvar)` | Dependent variable (inherited) |
| `e(kernel)` | Kernel function used |
| `e(bwselect)` | Bandwidth selection method |

**Matrices:**

| Result | Description |
|--------|-------------|
| `e(Estimate)` | Main results. 9 columns `(eval, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw)` for `type(dynamic)`/`type(group)`/`type(calendar)`; 8 columns `(z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw)` for `type(simple)` (no `eval` dimension) |
| `e(aggte_est)` | Aggregated point estimates |
| `e(aggte_se)` | Aggregated standard errors |
| `e(aggte_bw)` | Bandwidths used |
| `e(aggte_eval)` | Evaluation points |

## Methodology

### Doubly Robust Estimand

The conditional DR estimand based on the not-yet-treated group is:

```
DR(g,t)(z) = E[ (G_g / E[G_g | Z=z] - R(g,t) / E[R(g,t) | Z=z]) * (Y_t - Y_{g-1} - m(g,t)(X)) | Z=z ]
```

where:
- **G_g** is the group indicator
- **R(g,t)** is the inverse probability weight based on the generalized propensity score
- **m(g,t)(X)** is the outcome regression function
- The DR property ensures consistency if either the GPS or OR is correctly specified

Under the identification conditions (no anticipation, conditional parallel trends, overlap), CATT(g,t,z) = DR(g,t)(z).

### Three-Step Estimation

**Stage 1 — Parametric estimation:**
Estimate the GPS p(g,t)(X) and OR function m(g,t)(X) via logit and OLS, following Callaway and Sant'Anna (2021).

**Stage 2 — Intermediate variables:**
For each unit *i*, compute the doubly robust influence function:

```
Â(i,g,t) = (G(i,g) / μ̂_G(z) - R̂(i,g,t) / μ̂_R(z)) * (Y(i,t) - Y(i,g-1) - m̂(i,g,t))
```

where μ̂_G(z) and μ̂_R(z) are local polynomial regression estimates of E[G_g | Z=z] and E[R(g,t) | Z=z].

**Stage 3 — CATT estimation:**
Obtain the final CATT estimate via local polynomial regression of Â(i,g,t) on Z:

```
CATT_hat(g,t,z) = DR_hat(g,t)(z) = μ̂_A(z)
```

A common bandwidth *h* and kernel function *K* are used across all three nonparametric regressions (Stages 2–3), as required for uniform inference.

### Bandwidth Selection

| Method | Description |
|--------|-------------|
| **IMSE1** | Integrated MSE-optimal bandwidth, rule-of-thumb (default) |
| **IMSE2** | Integrated MSE-optimal bandwidth, plug-in estimator |
| **US1** | Uniform-smoothing bandwidth for uniform inference |
| **manual** | User-specified bandwidth via `bw()` |

For uniform inference, the bandwidth does not depend on (g, t, z). When `uniformall(true)`, the bandwidth is taken as the minimum across all (g,t) pairs to ensure simultaneous validity.

### Uniform Confidence Bands

The (1 − α) uniform confidence band for CATT(g,t,z) is:

```
[ CATT_hat(g,t,z) ± c * SE_hat(g,t,z) ]
```

where *c* is a critical value chosen to ensure:

```
P( CATT(g,t,z) ∈ CB for all (g,t,z) ∈ A ) ≥ 1 − α
```

Two critical values are provided:
- **c_hat** (analytical): Based on distributional approximation for suprema of Gaussian processes
- **c_check** (bootstrap): Based on multiplier bootstrap with Mammen's two-point distribution

The bootstrap critical value is generally preferred in finite samples.

### Aggregation

The aggregation parameters follow Section 5 of Imai, Qin, and Yanagi (2025). Each type computes a weighted average of CATT(g,t,z) over specific dimensions:

- **Dynamic**: θ_es(e, z) = Σ_g w(g,e) · CATT(g, g+e, z)
- **Group**: θ_group(g, z) = Σ_t w(g,t) · CATT(g, t, z)
- **Calendar**: θ_cal(t, z) = Σ_g w(g,t) · CATT(g, t, z)
- **Simple**: θ_simple(z) = Σ_{g,t} w(g,t) · CATT(g, t, z)

Aggregation uses its own bandwidth selection and local polynomial regression step.

## References

Imai, S., Qin, L., & Yanagi, T. (2025). Doubly robust uniform confidence bands for group-time conditional average treatment effects in difference-in-differences. *Journal of Business & Economic Statistics*, 1-13.

Callaway, B. & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.

Cheng, C. & Hoekstra, M. (2013). Does strengthening self-defense law deter crime or escalate violence? Evidence from expansions to castle doctrine. *Journal of Human Resources*, 48(3), 821-854.

Cohen, A., & Einav, L. (2003). The effects of mandatory seat belt laws on driving behavior and traffic fatalities. *Review of Economics and Statistics*, 85(4), 828-843.

Stevenson, B. & Wolfers, J. (2006). Bargaining in the shadow of the law: Divorce laws and family distress. *Quarterly Journal of Economics*, 121(1), 267-288.

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Shunsuke Imai**, Kyoto University
- **Lei Qin**, University of Niigata Prefecture
- **Takahide Yanagi**, Kyoto University

## License

AGPL-3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite both the methodology paper and the Stata implementation:

**APA Format:**

> Cai, X., & Xu, W. (2025). *didhetero: Stata module for heterogeneous treatment effects in staggered difference-in-differences* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/didhetero
>
> Imai, S., Qin, L., & Yanagi, T. (2025). Doubly robust uniform confidence bands for group-time conditional average treatment effects in difference-in-differences. *Journal of Business & Economic Statistics*, 1-13.

**BibTeX:**

```bibtex
@software{didhetero2025stata,
  title={didhetero: Stata module for heterogeneous treatment effects in staggered difference-in-differences},
  author={Xuanyu Cai and Wenli Xu},
  year={2025},
  version={0.1.0},
  url={https://github.com/gorgeousfish/didhetero}
}

@article{imai2025doubly,
  title={Doubly Robust Uniform Confidence Bands for Group-Time Conditional Average Treatment Effects in Difference-in-Differences},
  author={Imai, Shunsuke and Qin, Lei and Yanagi, Takahide},
  journal={Journal of Business \& Economic Statistics},
  pages={1--13},
  year={2025}
}
```
