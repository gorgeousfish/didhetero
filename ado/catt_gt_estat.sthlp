{smcl}
{* *! version 1.0.0  2026-06-30}{...}
{vieweralsosee "catt_gt" "help catt_gt"}{...}
{vieweralsosee "catt_gt_graph" "help catt_gt_graph"}{...}
{viewerjumpto "Syntax" "catt_gt_estat##syntax"}{...}
{viewerjumpto "Description" "catt_gt_estat##description"}{...}
{viewerjumpto "Options" "catt_gt_estat##options"}{...}
{viewerjumpto "Examples" "catt_gt_estat##examples"}{...}
{viewerjumpto "Stored results" "catt_gt_estat##results"}{...}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col:{bf:estat} {hline 2}}Post-estimation diagnostics for {cmd:catt_gt}{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Overlap diagnostics (GPS distribution)

{p 8 16 2}
{cmd:estat overlap} [{cmd:,} {opt th:reshold(#)}]


{pstd}
Pre-trend diagnostics (pre-treatment CATT estimates)

{p 8 16 2}
{cmd:estat pretrend}


{pstd}
These commands are available after {cmd:catt_gt}; see {help catt_gt}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:estat overlap} displays descriptive statistics of the Generalized Propensity
Score (GPS) values estimated during {cmd:catt_gt}. This diagnostic helps assess
whether Assumption 2.1 (overlap) is plausible: the GPS should be bounded away
from 0 and 1. Observations with extreme GPS values may indicate overlap
violations that could compromise estimation reliability.

{pstd}
{bf:Important:} {cmd:estat overlap} is a diagnostic display tool, not a formal
statistical test. No p-values or test statistics are produced. The paper does
not provide a formal overlap test statistic.

{pstd}
{cmd:estat pretrend} displays the CATT point estimates and confidence bands for
pre-treatment periods (t < g - anticipation). Under the identifying assumptions,
CATT_{g,t}(z) = 0 for pre-treatment periods. The display flags which (g,t,z)
cells have confidence bands that do not cover zero.

{pstd}
{bf:Important:} {cmd:estat pretrend} is a diagnostic reference only, not a
formal pre-test. The paper explicitly warns against conditioning subsequent
inference on a pre-test outcome (Section 5). Users should jointly report
pre-treatment and post-treatment estimates for transparency.


{marker options}{...}
{title:Options for estat overlap}

{phang}
{opt threshold(#)} specifies the boundary threshold for identifying extreme GPS
values. GPS values below {it:#} or above (1 - {it:#}) are flagged as near the
boundary. Default is {cmd:threshold(0.01)}. Must be between 0 and 0.5.


{marker examples}{...}
{title:Examples}

{pstd}Setup: run catt_gt estimation{p_end}
{phang2}{cmd:. catt_gt y, id(id) time(t) group(g) z(x) zeval(0.3 0.5 0.7)}{p_end}

{pstd}Overlap diagnostics with default threshold{p_end}
{phang2}{cmd:. estat overlap}{p_end}

{pstd}Overlap diagnostics with custom threshold{p_end}
{phang2}{cmd:. estat overlap, threshold(0.05)}{p_end}

{pstd}Pre-trend diagnostics (requires pretrend option){p_end}
{phang2}{cmd:. catt_gt y, id(id) time(t) group(g) z(x) zeval(0.3 0.5 0.7) pretrend}{p_end}
{phang2}{cmd:. estat pretrend}{p_end}

{pstd}Pre-trend diagnostics with bootstrap bands and custom level{p_end}
{phang2}{cmd:. catt_gt y, id(id) time(t) group(g) z(x) zeval(0.3 0.5 0.7) pretrend level(90) seed(42)}{p_end}
{phang2}{cmd:. estat pretrend}{p_end}
{phang2}{it:// Displays (g,t,z) cells where 90% CB does not cover zero}{p_end}
{phang2}{it:// A large number of flagged cells may suggest a violation of the}{p_end}
{phang2}{it:// parallel trends assumption at certain z values}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:estat overlap} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(gps_min)}}minimum GPS value{p_end}
{synopt:{cmd:r(gps_p5)}}5th percentile{p_end}
{synopt:{cmd:r(gps_p25)}}25th percentile{p_end}
{synopt:{cmd:r(gps_p50)}}median GPS value{p_end}
{synopt:{cmd:r(gps_p75)}}75th percentile{p_end}
{synopt:{cmd:r(gps_p95)}}95th percentile{p_end}
{synopt:{cmd:r(gps_max)}}maximum GPS value{p_end}
{synopt:{cmd:r(gps_mean)}}mean GPS value{p_end}
{synopt:{cmd:r(n_extreme)}}number of extreme observations{p_end}
{synopt:{cmd:r(pct_extreme)}}percentage of extreme observations{p_end}
{synopt:{cmd:r(n_obs)}}total number of GPS observations{p_end}
{synopt:{cmd:r(threshold)}}boundary threshold used{p_end}

{pstd}
{cmd:estat pretrend} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n_pretreat)}}number of pre-treatment (g,t,z) cells{p_end}
{synopt:{cmd:r(n_reject)}}number of cells where CB does not cover 0{p_end}
{synopt:{cmd:r(alpha)}}significance level used{p_end}
{synopt:{cmd:r(anticipation)}}anticipation parameter{p_end}


{title:Theoretical basis}

{pstd}
{bf:estat overlap:} Assumption 2.1 in the paper requires that the generalized
propensity score P(G=g|X) is bounded away from 0 and 1. This ensures that
the inverse probability weights used in the doubly-robust estimator remain
well-behaved. When GPS values approach 0 or 1, the corresponding weights
become extreme and estimation becomes unreliable.

{pstd}
{bf:estat pretrend:} Section 5 and the Supplementary Appendix (Eq. S5.1-S5.2)
define the testable implication that CATT_{g,t}(z) = 0 for pre-treatment
periods (t < g - anticipation). However, the paper explicitly states that
this should NOT be used as a formal pre-test for conditioning inference
(see discussion in Section 5). The recommended approach is to jointly display
pre-treatment and post-treatment results.


{title:Author}

{pstd}
Xinyu Chen
{p_end}
