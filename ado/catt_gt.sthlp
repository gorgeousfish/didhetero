{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "catt_gt##syntax"}{...}
{viewerjumpto "Description" "catt_gt##description"}{...}
{viewerjumpto "Options" "catt_gt##options"}{...}
{viewerjumpto "Remarks" "catt_gt##remarks"}{...}
{viewerjumpto "Examples" "catt_gt##examples"}{...}
{viewerjumpto "Stored results" "catt_gt##stored"}{...}
{viewerjumpto "References" "catt_gt##references"}{...}
{viewerjumpto "Authors" "catt_gt##authors"}{...}
{viewerjumpto "Also see" "catt_gt##alsosee"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{cmd:catt_gt} {hline 2}}Estimate group-time conditional average treatment effects on the treated{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:catt_gt} {depvar} {ifin}{cmd:,}
{opt id(varname)}
{opt time(varname)}
{opt group(varname)}
{opt z(varname)}
{opt zeval(numlist)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth id(varname)}}panel individual identifier{p_end}
{synopt:{opth time(varname)}}time period variable{p_end}
{synopt:{opth group(varname)}}treatment group variable; 0 = never-treated{p_end}
{synopt:{opth z(varname)}}continuous pre-treatment covariate{p_end}
{synopt:{opt zeval(numlist)}}evaluation points for {it:z}{p_end}

{syntab:Covariates}
{synopt:{opt xformula(string)}}complete covariate specification for outcome regression{p_end}

{syntab:Estimation}
{synopt:{opt porder(#)}}local polynomial order; {cmd:1} or {cmd:2}; default is {cmd:2}{p_end}
{synopt:{opt kernel(string)}}kernel function; {cmd:gau} or {cmd:epa}; default is {cmd:gau}{p_end}
{synopt:{opt bwselect(string)}}bandwidth selection method; default is {cmd:IMSE1}{p_end}
{synopt:{opt bw(numlist)}}manual bandwidth; requires {cmd:bwselect(manual)}{p_end}

{syntab:Inference}
{synopt:{opt alp(#)}}significance level; default is {cmd:0.05}{p_end}
{synopt:{opt bstrap(true|false)}}multiplier bootstrap toggle; default is {cmd:true}{p_end}
{synopt:{opt biters(#)}}bootstrap iterations; default is {cmd:1000}{p_end}
{synopt:{opt seed(#)}}RNG seed for bootstrap reproducibility; default is {cmd:-1} (current RNG state){p_end}
{synopt:{opt uniformall(true|false)}}joint uniform-band toggle; default is {cmd:true}{p_end}

{syntab:Panel}
{synopt:{opt gteval(numlist)}}(g,t) evaluation pairs: {it:g1 t1 g2 t2 ...}; default is all identifiable pairs{p_end}
{synopt:{opt control_group(string)}}{cmd:notyettreated} or {cmd:nevertreated}; default is {cmd:notyettreated}{p_end}
{synopt:{opt anticipation(#)}}anticipation periods; default is {cmd:0}{p_end}
{synopt:{opt pretrend}}include pre-treatment periods for testing{p_end}

{syntab:Legacy flags}
{synopt:{opt bstrap}}alias for {cmd:bstrap(true)}{p_end}
{synopt:{opt nobstrap}}alias for {cmd:bstrap(false)}{p_end}
{synopt:{opt nouniformall}}alias for {cmd:uniformall(false)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:catt_gt} estimates group-time conditional average treatment effects on
the treated, CATT(g,t,z), as a function of a continuous pre-treatment
covariate {it:z}. The estimator is doubly robust: it combines inverse
probability weighting with outcome regression through local polynomial
smoothing. Results include pointwise and uniform confidence bands over
{it:z} within each (g,t) pair, and (optionally) joint uniform bands over
(g,t,z).

{pstd}
The command implements Section 4 of Imai, Qin, and Yanagi (2025). It shares
the same estimation engine as {helpb didhetero}; the two commands are
interchangeable for estimation and differ only in legacy flag naming. After
estimation, use {helpb aggte_gt} to aggregate into summary parameters and
{helpb catt_gt_graph} to visualize the estimated CATT functions.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth id(varname)} specifies the variable that identifies individual panel
units. Each unit must appear in every time period (balanced panel required).

{phang}
{opth time(varname)} specifies the time period variable. Must contain integer
values.

{phang}
{opth group(varname)} specifies the treatment group variable, defined as the
first period in which a unit receives treatment. Units with {cmd:group = 0} are
never-treated. Must be time-invariant within each {opt id()}.

{phang}
{opth z(varname)} specifies the continuous pre-treatment covariate along which
treatment effect heterogeneity is estimated. Must be time-invariant within each
{opt id()}.

{phang}
{opt zeval(numlist)} specifies the evaluation points at which the CATT function
is estimated. These should lie within the support of {it:z}. The command sorts
{opt zeval()} into ascending order before estimation, and both {cmd:e(zeval)}
and the rows of {cmd:e(results)} are returned in that ascending-order {it:z}
layout rather than the raw input order.

{dlgtab:Covariates}

{phang}
{opt xformula(string)} specifies the {it:complete} covariate set for the
first-stage outcome regression and the generalized propensity score.
When {opt xformula()} is omitted, the design matrix defaults to
{it:X} = [1, {it:z}]. Once specified, {opt xformula()} fully replaces the
default set. A variable list or formula with an intercept (e.g.,
{cmd:xformula(Z)}, {cmd:xformula("~ Z + X1")}) produces
{it:X} = [1, columns given in xformula]. An explicit no-intercept formula
(e.g., {cmd:xformula("~ 0 + Z")} or {cmd:xformula("~ Z - 1")}) suppresses
the constant column.
{opt xformula()} must explicitly include {it:z} (e.g., {cmd:xformula(Z)});
otherwise the command exits with an error.
To add extra covariates, ensure they exist in the dataset and list them
together with {cmd:Z} in {opt xformula()}.
As with {opt z()}, every covariate in {opt xformula()} must be time-invariant
within each {opt id()} (pre-treatment covariates {it:X_i} in the paper);
time-varying variables cause an error.
Supported syntax includes bare variables, {cmd:I(expr)}, {cmd::} (interaction
only), {cmd:*} (main effects plus interaction), and the no-intercept tokens
{cmd:0 +} / {cmd:- 1}. General term removal beyond {cmd:- 1} is not supported.
For example, {cmd:xformula("Z*X1")} is equivalent to
{cmd:xformula("Z + X1 + Z:X1")}; {cmd:*} expands into main effects and
interactions. To include only the interaction term itself, write
{cmd:xformula("Z + Z:X1")}; {cmd::} adds the interaction without
automatically appending the right-hand-side main effect.
When {opt xformula()} is omitted the default {it:X} = [1, {it:z}] is a
convenience entry point. To match the full specification in the paper,
explicitly write {cmd:xformula(Z)} or a more general formula.

{dlgtab:Estimation}

{phang}
{opt porder(#)} sets the order of the local polynomial regression. Must be
{cmd:1} (local linear) or {cmd:2} (local quadratic). Default is {cmd:2}.

{phang}
{opt kernel(string)} specifies the kernel function for local polynomial
smoothing. Options are {cmd:gau} (Gaussian) and {cmd:epa} (Epanechnikov).
Default is {cmd:gau}.

{phang}
{opt bwselect(string)} specifies the bandwidth selection method. Options
are {cmd:IMSE1}, {cmd:IMSE2}, {cmd:US1}, and {cmd:manual}. Default is
{cmd:IMSE1}. When {cmd:manual} is specified, {opt bw()} must also be
provided. {cmd:IMSE1}/{cmd:IMSE2} are IMSE-optimal bandwidths for local
linear and local quadratic regressions, respectively; {cmd:US1} is a
rule-of-thumb undersmoothing bandwidth for the local linear fit and is
recommended when uniform inference is the primary goal.

{phang}
{opt bw(numlist)} specifies manual bandwidth value(s) and is valid only
when {cmd:bwselect(manual)}. All values must be positive. A scalar
{opt bw()} is always accepted. A vector {opt bw()} is accepted only when
{opt gteval()} is supplied, and its length must equal the number of
specified (g,t) pairs. When {opt gteval()} contains a single (g,t) pair,
the effective inference domain degenerates to uniform inference over
{it:z} alone.

{dlgtab:Inference}

{phang}
{opt alp(#)} sets the significance level for confidence intervals. Must
be strictly between 0 and 1. Default is {cmd:0.05}.

{phang}
{opt bstrap(true|false)} toggles the multiplier bootstrap for uniform
inference. Bootstrap is enabled by default; omit the option or pass
{cmd:bstrap(true)} to keep it on, and {cmd:bstrap(false)} to disable it.

{phang}
{opt biters(#)} sets the number of bootstrap iterations. Default is
{cmd:1000}. Must be a positive integer when the bootstrap is enabled.

{phang}
{opt seed(#)} sets the Stata RNG seed immediately before the bootstrap
weights are drawn. The default {cmd:seed(-1)} is a sentinel that leaves
the current RNG state unchanged; the effective seed is recorded in
{cmd:e(seed_request)} and {cmd:e(seed)}.

{phang}
{opt uniformall(true|false)} toggles joint uniform inference over
(g,t,z). The default {cmd:uniformall(true)} delivers a single critical
value that covers every evaluated (g,t) pair and every {it:z} grid point
simultaneously. {cmd:uniformall(false)} narrows coverage to bands that
are uniform over {it:z} within each (g,t) pair.

{dlgtab:Panel}

{phang}
{opt gteval(numlist)} specifies the (g,t) pairs to evaluate, in the
format {it:g1 t1 g2 t2 ...}. When omitted, {cmd:catt_gt} evaluates every
identifiable (g,t) pair implied by the sample support,
{opt control_group()}, and {opt anticipation()}.

{phang}
{opt control_group(string)} selects the comparison group.
{cmd:notyettreated} uses units not yet treated by period {it:t} as
controls. {cmd:nevertreated} restricts controls to never-treated units
(group = 0) and requires at least one such unit in the sample. Default
is {cmd:notyettreated}.

{phang}
{opt anticipation(#)} specifies the number of periods before treatment
in which units may anticipate the intervention. The excluded long-
difference baseline is period t = g - anticipation - 1; this reduces to
event_time = -1 when {cmd:anticipation(0)}. Default is {cmd:0}.

{phang}
{opt pretrend} includes pre-treatment periods in the estimation so that
CATT estimates can be plotted against a zero reference line for visual
pre-trends diagnostics. Only applicable when {opt gteval()} is omitted,
because an explicit {opt gteval()} already fixes the evaluation pairs.

{dlgtab:Legacy flags}

{phang}
{opt bstrap}, {opt nobstrap}, and {opt nouniformall} are legacy bare-flag
aliases for {cmd:bstrap(true)}, {cmd:bstrap(false)}, and
{cmd:uniformall(false)} respectively. They are retained for backward
compatibility with earlier versions of the package. Do not combine a
legacy flag with the corresponding {cmd:true|false} form.


{marker remarks}{...}
{title:Remarks}

{pstd}
Identification follows Callaway and Sant'Anna (2021) under a conditional
parallel trends assumption; heterogeneity in the continuous covariate
{it:z} is identified nonparametrically by the three-stage construction of
Imai, Qin, and Yanagi (2025). The first stage fits the generalized
propensity score and outcome regression parametrically; the second and
third stages estimate the doubly robust influence function and the CATT
surface by local polynomial regression. The covariate vector X used in
the first stage is controlled by {opt xformula()}; pass {cmd:xformula(Z)}
to match the default specification {it:X} = [1, {it:z}] explicitly.

{pstd}
{cmd:catt_gt} and {helpb didhetero} share the same estimation engine.
Both expose the modern controls {opt gteval()}, {opt seed()}, and
{opt bstrap(true|false)}. {cmd:catt_gt} keeps the legacy {opt nobstrap}
alias, while {helpb didhetero} keeps the symmetric {opt nobootstrap}
alias. Either command can feed {helpb aggte_gt} for aggregation and
{helpb catt_gt_graph} for visualization.

{pstd}
For persistence with {cmd:estimates store} and {cmd:estimates restore},
{cmd:catt_gt} posts a storable {cmd:e(b)} shell and a conformable zero
{cmd:e(V)}. Generic covariance-based postestimation commands such as
{cmd:lincom}, {cmd:test}, and {cmd:testparm} are therefore not a supported
inferential workflow. Use {cmd:e(results)}, {cmd:e(catt_se)}, and the
reported confidence bands for inference.


{marker examples}{...}
{title:Examples}

{pstd}Setup: generate simulation data.{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}

{pstd}Example 1: Deterministic path (analytical bands only).{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false)}{p_end}
{phang2}{cmd:. matrix list e(results)}{p_end}

{pstd}Example 2: Local linear estimator with a manual bandwidth.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.5 0 0.5) gteval(2 2) xformula(Z) porder(1) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.45) bstrap(false) uniformall(false)}{p_end}

{pstd}Example 3: Interaction syntax in {opt xformula()}.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.5 0 0.5) gteval(2 2 2 3 3 3) ///}{p_end}
{phang2}{cmd:        xformula("Z*X1") bwselect(manual) bw(0.45 0.45 0.45) ///}{p_end}
{phang2}{cmd:        bstrap(false) uniformall(false)}{p_end}
{phang2}{it:// xformula("Z*X1") expands to xformula("Z + X1 + Z:X1")}{p_end}
{phang2}{it:// Use xformula("Z + Z:X1") for the interaction only}{p_end}

{pstd}Example 4: Estimate and plot CATT curves.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 5: Pre-trends testing with zero reference line.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) pretrend bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 6: Selected (g,t) pairs with a bandwidth vector.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) gteval(2 2 3 3) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.30 0.45) bstrap(false) uniformall(false)}{p_end}

{pstd}Example 7: Reproducible bootstrap with an explicit seed.{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 0 0.8) biters(500) seed(42)}{p_end}


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:catt_gt} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of cross-sectional units{p_end}
{synopt:{cmd:e(T)}}number of time periods{p_end}
{synopt:{cmd:e(gbar)}}upper bound for treatment-time support; {cmd:.} encodes {cmd:+Inf} when never-treated units exist{p_end}
{synopt:{cmd:e(gbar_isinf)}}1 if {cmd:e(gbar)} encodes {cmd:+Inf}, 0 otherwise{p_end}
{synopt:{cmd:e(num_gteval)}}number of (g,t) pairs evaluated{p_end}
{synopt:{cmd:e(num_zeval)}}number of z evaluation points{p_end}
{synopt:{cmd:e(porder)}}polynomial order{p_end}
{synopt:{cmd:e(anticipation)}}anticipation periods{p_end}
{synopt:{cmd:e(anticip)}}anticipation periods (legacy alias){p_end}
{synopt:{cmd:e(alp)}}significance level{p_end}
{synopt:{cmd:e(bstrap)}}1 if bootstrap was performed, 0 otherwise{p_end}
{synopt:{cmd:e(biters)}}effective bootstrap iterations; 0 if bootstrap is disabled{p_end}
{synopt:{cmd:e(seed_request)}}requested bootstrap seed; {cmd:-1} means use the current RNG state if bootstrap runs{p_end}
{synopt:{cmd:e(seed)}}effective command-level bootstrap seed; missing if bootstrap is disabled or no command-level seed was applied{p_end}
{synopt:{cmd:e(uniformall)}}1 if joint uniform bands requested, 0 otherwise{p_end}
{synopt:{cmd:e(pretrend)}}1 if pre-trends testing was requested, 0 otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:catt_gt}{p_end}
{synopt:{cmd:e(depvar)}}outcome variable name{p_end}
{synopt:{cmd:e(idvar)}}panel id variable name{p_end}
{synopt:{cmd:e(timevar)}}time variable name{p_end}
{synopt:{cmd:e(groupvar)}}group variable name{p_end}
{synopt:{cmd:e(zvar)}}covariate variable name{p_end}
{synopt:{cmd:e(control_group)}}control group specification{p_end}
{synopt:{cmd:e(control)}}control group specification (legacy alias){p_end}
{synopt:{cmd:e(kernel)}}kernel function used{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection method used{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}storable estimate row used by {cmd:estimates store}; inference uses {cmd:e(results)} and {cmd:e(catt_se)} instead{p_end}
{synopt:{cmd:e(V)}}conformable zero matrix posted for persistence; no covariance inference is supported{p_end}
{synopt:{cmd:e(results)}}main results matrix with columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw; rows ordered by ascending {it:z} within each (g,t){p_end}
{synopt:{cmd:e(Estimate)}}alias of {cmd:e(results)}{p_end}
{synopt:{cmd:e(Estimate_b)}}vectorized point estimates used by the storable {cmd:e(b)} shell{p_end}
{synopt:{cmd:e(gteval)}}evaluated (g,t) pairs (K x 2){p_end}
{synopt:{cmd:e(zeval)}}evaluation points for {it:z}, sorted ascending{p_end}
{synopt:{cmd:e(bw)}}bandwidth per (g,t) pair{p_end}
{synopt:{cmd:e(c_hat)}}analytical joint critical values per (g,t) pair{p_end}
{synopt:{cmd:e(c_check)}}multiplier bootstrap critical values per (g,t) pair; only when {cmd:e(bstrap)} is 1{p_end}
{synopt:{cmd:e(catt_est)}}CATT point estimates{p_end}
{synopt:{cmd:e(catt_se)}}CATT standard errors{p_end}
{synopt:{cmd:e(kd0_Z)}}kernel density estimates at {it:z}{p_end}
{synopt:{cmd:e(kd1_Z)}}kernel density derivative estimates at {it:z}{p_end}
{synopt:{cmd:e(Z_supp)}}support of {it:Z}{p_end}

{pstd}
The remaining matrices capture the panel snapshot used by {helpb aggte_gt}
for two-pass aggregation and are not part of the user-facing inference
output:

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(B_g_t)}}flattened influence function matrix B_{i,g,t}(z_r){p_end}
{synopt:{cmd:e(G_g)}}group indicator matrix (n x K); G_{i,g} = 1 if unit {it:i} belongs to group {it:g}{p_end}
{synopt:{cmd:e(mu_G_g)}}conditional mean E[G_g | Z = z_r] at the eval grid (R x K){p_end}
{synopt:{cmd:e(Z)}}row vector of unit-level {it:z} values{p_end}
{synopt:{cmd:e(dh_Y_wide)}}wide-format outcome matrix (n x T){p_end}
{synopt:{cmd:e(dh_G_unit)}}row vector of unit-level treatment group assignments{p_end}
{synopt:{cmd:e(dh_t_vals)}}row vector of distinct time-period values{p_end}
{synopt:{cmd:e(dh_gps_mat)}}stacked generalized propensity score fits by (g, unit){p_end}
{synopt:{cmd:e(dh_or_mat)}}stacked outcome regression fits by (g, unit){p_end}


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
Help: {helpb didhetero}, {helpb aggte_gt}, {helpb catt_gt_graph}, {helpb didhetero_simdata}
{p_end}
