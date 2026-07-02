{smcl}
{* *! version 1.0.0  2026-06-30}{...}
{vieweralsosee "catt_gt" "help catt_gt"}{...}
{vieweralsosee "catt_gt_estat" "help catt_gt_estat"}{...}
{vieweralsosee "catt_gt_graph" "help catt_gt_graph"}{...}
{viewerjumpto "Syntax" "catt_gt_predict##syntax"}{...}
{viewerjumpto "Description" "catt_gt_predict##description"}{...}
{viewerjumpto "Options" "catt_gt_predict##options"}{...}
{viewerjumpto "Examples" "catt_gt_predict##examples"}{...}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col:{bf:predict} {hline 2}}Post-estimation predictions for {cmd:catt_gt}{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:predict} {newvar} [{cmd:,} {it:statistic}]


{synoptset 12 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt se}}standard errors{p_end}
{synopt:{opt ci1}}analytical confidence interval (generates {it:newvar}{cmd:_lb} and {it:newvar}{cmd:_ub}){p_end}
{synopt:{opt ci2}}bootstrap confidence interval (generates {it:newvar}{cmd:_lb} and {it:newvar}{cmd:_ub}){p_end}
{synopt:{opt bw}}bandwidth used{p_end}
{synopt:{opt zval}}z evaluation points{p_end}
{synopt:{opt gval}}group values{p_end}
{synopt:{opt tval}}time period values{p_end}
{synoptline}

{pstd}
These statistics are available after {cmd:catt_gt} or {cmd:didhetero};
see {helpb catt_gt} and {helpb didhetero}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} after {cmd:catt_gt} (or {cmd:didhetero}) extracts components of
the estimation results matrix {cmd:e(results)} into new variables in the
current dataset. By default (no option), {cmd:predict} generates the CATT
point estimates (paper Eq. 2.3, Lemma 2).

{pstd}
Each generated variable contains one value per evaluation point in
observations 1 through {it:K}, where {it:K} = {cmd:e(num_gteval)} *
{cmd:e(num_zeval)} is the total number of rows in {cmd:e(results)}. Remaining
observations are set to missing.

{pstd}
Only one statistic may be specified at a time. The {opt ci1} and {opt ci2}
options each generate two variables ({it:newvar}{cmd:_lb} and
{it:newvar}{cmd:_ub}) for the lower and upper bounds of the confidence
interval.


{marker options}{...}
{title:Options}

{phang}
{opt se} stores the standard error of the CATT estimate at each evaluation
point (column 5 of {cmd:e(results)}).

{phang}
{opt ci1} stores the analytical (distributional approximation) confidence
interval bounds. Generates two variables: {it:newvar}{cmd:_lb} (lower bound)
and {it:newvar}{cmd:_ub} (upper bound). These correspond to columns 6 and 7
of {cmd:e(results)}.

{phang}
{opt ci2} stores the multiplier bootstrap confidence interval bounds. Generates
two variables: {it:newvar}{cmd:_lb} (lower bound) and {it:newvar}{cmd:_ub}
(upper bound). These correspond to columns 8 and 9 of {cmd:e(results)}.
Requires that the preceding estimation was run with bootstrap enabled
({cmd:bstrap(true)}, the default); otherwise an error is issued.

{phang}
{opt bw} stores the bandwidth used for local polynomial smoothing at each
evaluation point (column 10 of {cmd:e(results)}).

{phang}
{opt zval} stores the z evaluation point values (column 3 of
{cmd:e(results)}). Useful for merging results with external datasets.

{phang}
{opt gval} stores the group value for each evaluation point (column 1 of
{cmd:e(results)}).

{phang}
{opt tval} stores the time period value for each evaluation point (column 2
of {cmd:e(results)}).


{marker examples}{...}
{title:Examples}

{pstd}Setup: run estimation{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 0 0.8)}{p_end}

{pstd}Extract CATT point estimates (default){p_end}
{phang2}{cmd:. predict catt_hat}{p_end}

{pstd}Extract standard errors{p_end}
{phang2}{cmd:. predict catt_se, se}{p_end}

{pstd}Extract bootstrap confidence intervals{p_end}
{phang2}{cmd:. predict catt_ci, ci2}{p_end}
{phang2}{cmd:. list catt_ci_lb catt_ci_ub in 1/5}{p_end}

{pstd}Extract bandwidth values{p_end}
{phang2}{cmd:. predict h, bw}{p_end}

{pstd}Extract evaluation grid coordinates{p_end}
{phang2}{cmd:. predict z_pts, zval}{p_end}
{phang2}{cmd:. predict g_pts, gval}{p_end}
{phang2}{cmd:. predict t_pts, tval}{p_end}

{pstd}Combine for a custom results frame{p_end}
{phang2}{cmd:. predict my_est}{p_end}
{phang2}{cmd:. predict my_se, se}{p_end}
{phang2}{cmd:. predict my_z, zval}{p_end}
{phang2}{cmd:. predict my_g, gval}{p_end}
{phang2}{cmd:. list my_g my_z my_est my_se in 1/10}{p_end}


{title:Author}

{pstd}
Xinyu Chen
{p_end}
