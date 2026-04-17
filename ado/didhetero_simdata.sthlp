{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "didhetero_simdata##syntax"}{...}
{viewerjumpto "Description" "didhetero_simdata##description"}{...}
{viewerjumpto "Options" "didhetero_simdata##options"}{...}
{viewerjumpto "Output variables" "didhetero_simdata##variables"}{...}
{viewerjumpto "Remarks" "didhetero_simdata##remarks"}{...}
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
{synopt:{opt n(#)}}number of cross-sectional units; must be a positive integer{p_end}
{synopt:{opt tau(#)}}number of time periods; must be an integer > 2{p_end}

{syntab:DGP}
{synopt:{opt dgpy(#)}}DGP for the treated potential outcome; {cmd:1} or {cmd:2}; default is {cmd:1}{p_end}
{synopt:{opt hc}}heteroscedastic error terms; default is homoscedastic{p_end}
{synopt:{opt dimx(#)}}dimension of covariate vector X; default is {cmd:1}{p_end}
{synopt:{opt discrete}}generate discrete Z in {c -(}-1, 0, 1{c )-} instead of continuous N(0, 1){p_end}

{syntab:RNG and data}
{synopt:{opt seed(#)}}RNG seed; {cmd:-1} keeps the current RNG state; nonnegative integers otherwise{p_end}
{synopt:{opt clear}}replace data currently in memory{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:didhetero_simdata} generates a balanced panel dataset following the
data-generating process (DGP) described in Section 6 of Imai, Qin, and
Yanagi (2025). The resulting dataset is intended for testing and
demonstrating {helpb catt_gt} and {helpb didhetero} under known data
generation.

{pstd}
The DGP creates a staggered adoption design over {it:tau} time periods:
group {it:g} (for g = 2, ..., tau) first receives treatment in period
{it:g} and group 0 is never treated. Treatment effects vary with a
continuous covariate {cmd:Z}. The treated potential outcome is

{phang2}
{opt dgpy(1)}: nonlinear term {cmd:(G / t) * sin(pi * Z)} (default){p_end}
{phang2}
{opt dgpy(2)}: linear interaction term {cmd:Z * G / t}{p_end}

{pstd}
All random draws come from the current Stata RNG. Pass {opt seed()} for
reproducible output. When additional covariates are requested via
{opt dimx()}, the extra columns are i.i.d. standard normal and are
written to variables {cmd:X1}, {cmd:X2}, ..., {cmd:X(dimx - 1)} in the
output dataset.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt n(#)} specifies the number of cross-sectional units (individuals).
Must be a positive integer. The total number of observations generated is
{it:n} x {it:tau}.

{phang}
{opt tau(#)} specifies the number of time periods. Must be an integer
strictly greater than 2.

{dlgtab:DGP}

{phang}
{opt dgpy(#)} selects the DGP for the treated potential outcome.
{cmd:dgpy(1)} uses the nonlinear term {cmd:(G / t) * sin(pi * Z)};
{cmd:dgpy(2)} uses the linear interaction {cmd:Z * G / t}. Default is
{cmd:1}.

{phang}
{opt hc} generates heteroscedastic error terms where the variance
of the idiosyncratic shock depends on {cmd:Z} and the treatment group.
Without this option, errors are homoscedastic standard normal.

{phang}
{opt dimx(#)} specifies the dimension of the covariate vector X. When
{cmd:dimx(1)}, only {cmd:Z} is generated (default). With {cmd:dimx(k)}
and k > 1, additional i.i.d. standard normal covariates
{cmd:X1}, ..., {cmd:X(k - 1)} are appended.

{phang}
{opt discrete} generates a discrete covariate Z taking values in
{c -(}-1, 0, 1{c )-} with equal probability instead of the default
continuous Z ~ N(0, 1). Because {helpb catt_gt} / {helpb didhetero}
require a continuous {cmd:z()} variable, the discrete option is intended
for data-inspection and diagnostic use only; do not feed the resulting
dataset to the main estimation commands.

{dlgtab:RNG and data}

{phang}
{opt seed(#)} sets the Stata RNG seed for reproducibility. The default
{cmd:seed(-1)} leaves the current RNG state unchanged. Any other value
must be a nonnegative integer and is applied only after the Mata backend
has been confirmed to load. After a successful run,
{cmd:r(seed_explicit)} records whether {opt seed()} was provided and
{cmd:r(seed)} stores the explicit value when one was used.

{phang}
{opt clear} replaces the dataset currently in memory. When variables or
observations already exist in memory and {opt clear} is not specified,
the command exits with error code 4 ({cmd:no; data in memory would be lost}).
An empty in-memory dataset still counts as existing data.


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


{marker remarks}{...}
{title:Remarks}

{pstd}
The generated dataset is sorted by {cmd:id} and {cmd:period} on return.
All columns have descriptive variable labels so that {cmd:describe}
displays the panel structure immediately after the command completes.
When {opt seed()} is provided, the RNG seed is set {it:after} the Mata
backend has been confirmed to load; this keeps the generated seed stream
independent of any RNG draws consumed while loading the backend.

{pstd}
Although the {cmd:discrete} option is available for diagnostic
purposes, the main estimation commands {helpb catt_gt} and
{helpb didhetero} require a continuous {cmd:z()} variable: the
underlying nonparametric estimator relies on kernel density estimates
of Z that are not defined for a discrete support. Use
{opt discrete} only to inspect the group-assignment mechanism or to
debug third-party code.


{marker examples}{...}
{title:Examples}

{pstd}Example 1: Basic usage with continuous Z.{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. describe}{p_end}
{phang2}{cmd:. summarize}{p_end}

{pstd}Example 2: Heteroscedastic errors with {cmd:dgpy(2)}.{p_end}
{phang2}{cmd:. didhetero_simdata, n(1000) tau(6) seed(99) hc dgpy(2) clear}{p_end}
{phang2}{cmd:. summarize Y, detail}{p_end}

{pstd}Example 3: Multiple pre-treatment covariates.{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) dimx(3) clear}{p_end}
{phang2}{cmd:. describe}{p_end}
{phang2}{cmd:. catt_gt Y, id(id) time(period) group(G) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.5 0 0.5) xformula("Z + X1 + X2") bstrap(false)}{p_end}

{pstd}Example 4: Inspect the discrete-Z DGP (not for main estimation).{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(1) discrete clear}{p_end}
{phang2}{cmd:. tab Z}{p_end}
{phang2}{cmd:. tab G Z}{p_end}


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
