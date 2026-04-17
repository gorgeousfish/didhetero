{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "didhetero_simdata##syntax"}{...}
{viewerjumpto "Description" "didhetero_simdata##description"}{...}
{viewerjumpto "Options" "didhetero_simdata##options"}{...}
{viewerjumpto "Output variables" "didhetero_simdata##variables"}{...}
{viewerjumpto "Examples" "didhetero_simdata##examples"}{...}
{viewerjumpto "Stored results" "didhetero_simdata##stored"}{...}
{viewerjumpto "Authors" "didhetero_simdata##authors"}{...}
{viewerjumpto "References" "didhetero_simdata##references"}{...}
{viewerjumpto "Also see" "didhetero_simdata##alsosee"}{...}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col:{cmd:didhetero_simdata} {hline 2}}Generate simulation data for heterogeneous DID estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:didhetero_simdata}{cmd:,}
{opt n(#)}
{opt tau(#)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt n(#)}}number of cross-sectional units{p_end}
{synopt:{opt tau(#)}}number of time periods; must be > 2{p_end}

{syntab:Optional}
{synopt:{opt seed(#)}}random seed; {cmd:-1} keeps the current RNG state; other values must be nonneg. integers{p_end}
{synopt:{opt hc}}use heteroscedastic error terms{p_end}
{synopt:{opt dimx(#)}}dimension of covariates X; default is {cmd:1}{p_end}
{synopt:{opt dgpy(#)}}DGP for treated potential outcome; 1 or 2; default is {cmd:1}{p_end}
{synopt:{opt discrete}}use discrete Z in {c -(}-1, 0, 1{c )-} instead of continuous{p_end}
{synopt:{opt clear}}allow replacing data in memory{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:didhetero_simdata} generates a balanced panel dataset following the
data-generating process (DGP) described in Section 6 of Imai,
Qin, and Yanagi (2025). The default continuous-{cmd:Z} DGP can be used for
testing and demonstrating {helpb catt_gt} and {helpb didhetero} on the current
main estimation path.

{pstd}
The DGP creates a staggered adoption design with {it:tau} time periods.
Treatment groups are assigned so that group {it:g} (for g = 2, ..., tau)
first receives treatment in period {it:g}, and group 0 is never treated.
Treatment effects are heterogeneous in the continuous covariate Z.

{pstd}
When {opt dgpy(1)} is specified (default), the treated potential outcome
uses the nonlinear term {cmd:(G / t) * sin(pi * Z)}. When {opt dgpy(2)} is
specified, the treated potential outcome uses the linear interaction term
{cmd:Z * G / t}.

{pstd}
This command implements the data-generating process described in Section 6 of
Imai, Qin, and Yanagi (2025).


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt n(#)} specifies the number of cross-sectional units (individuals). Must
be a positive integer.

{phang}
{opt tau(#)} specifies the number of time periods. Must be greater than 2.
The total number of observations generated is {it:n} x {it:tau}.

{dlgtab:Optional}

{phang}
{opt seed(#)} sets the random number seed for reproducibility. The special
value {cmd:-1} leaves the current Stata RNG state unchanged. Any other value
must be a nonnegative integer and is applied only after the Mata backend has
been confirmed available. After a successful run, {cmd:r(seed_explicit)}
records whether {opt seed()} was provided and {cmd:r(seed)} stores the
explicit seed value when one was used.

{phang}
{opt hc} requests heteroscedastic error terms. By default, errors are
homoscedastic.

{phang}
{opt dimx(#)} specifies the dimension of the covariate vector X. When
{cmd:dimx(1)}, only Z is generated. When {cmd:dimx(k)} with k > 1, additional
covariates X1, X2, ..., X(k-1) are generated. Default is {cmd:1}.

{phang}
{opt dgpy(#)} selects the DGP for the treated potential outcome. Must be 1 or
2. Default is {cmd:1}.

{phang}
{opt discrete} generates Z as a discrete variable taking values in
{c -(}-1, 0, 1{c )-} instead of the default continuous Z ~ N(0,1). The main
estimation commands {helpb catt_gt} / {helpb didhetero} require a continuous
{cmd:z()}; this option is intended for data-inspection purposes only and should
not be used as an estimation example.

{phang}
{opt clear} permits replacing data currently in memory. If data exists in
memory and {opt clear} is not specified, the command exits with an error.


{marker variables}{...}
{title:Output variables}

{synoptset 12 tabbed}{...}
{p2col:Variable}Description{p_end}
{synoptline}
{synopt:{cmd:id}}individual identifier{p_end}
{synopt:{cmd:period}}time period{p_end}
{synopt:{cmd:Y}}observed outcome{p_end}
{synopt:{cmd:G}}treatment group (0 = never-treated, 2..tau = first treated in period g){p_end}
{synopt:{cmd:Z}}covariate (continuous or discrete depending on options){p_end}
{synopt:{cmd:X1}, {cmd:X2}, ...}additional covariates when {cmd:dimx} > 1{p_end}
{synoptline}


{marker examples}{...}
{title:Examples}

{pstd}Example 1: Basic usage{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. describe}{p_end}
{phang2}{cmd:. summarize}{p_end}

{pstd}Example 2: Heteroscedastic errors with DGP-Y=2{p_end}
{phang2}{cmd:. didhetero_simdata, n(1000) tau(6) seed(99) hc dgpy(2) clear}{p_end}
{phang2}{cmd:. summarize Y, detail}{p_end}

{pstd}Example 3: Multiple covariates{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) dimx(3) clear}{p_end}
{phang2}{cmd:. describe}{p_end}


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:didhetero_simdata} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of cross-sectional units{p_end}
{synopt:{cmd:r(tau)}}number of time periods{p_end}
{synopt:{cmd:r(dimx)}}covariate dimension{p_end}
{synopt:{cmd:r(dgpy)}}DGP type{p_end}
{synopt:{cmd:r(hc)}}1 if heteroscedastic, 0 otherwise{p_end}
{synopt:{cmd:r(continuous)}}1 if continuous Z, 0 if discrete{p_end}
{synopt:{cmd:r(seed_explicit)}}1 if {opt seed()} was specified, 0 otherwise{p_end}
{synopt:{cmd:r(seed)}}explicit seed value; missing when the caller RNG state was used{p_end}
{synopt:{cmd:r(n_groups)}}number of groups (including never-treated){p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(groups)}}list of group values{p_end}


{marker references}{...}
{title:References}

{phang}
Callaway, B. and P. H. C. Sant'Anna. 2021.
"Difference-in-Differences with Multiple Time Periods."
{it:Journal of Econometrics} 225(2): 200-230.
{p_end}

{phang}
Imai, S., L. Qin, and T. Yanagi. 2025.
"Doubly Robust Uniform Confidence Bands for Group-Time Conditional
Average Treatment Effects in Difference-in-Differences."
{it:Journal of Business & Economic Statistics}, 1-13.
{p_end}


{marker authors}{...}
{title:Authors}

{pstd}
Xuanyu Cai, City University of Macau{break}
Wenli Xu, City University of Macau{p_end}


{marker alsosee}{...}
{title:Also see}

{psee}
Help: {helpb didhetero}, {helpb catt_gt}, {helpb aggte_gt}, {helpb catt_gt_graph}
{p_end}
