{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "catt_gt##syntax"}{...}
{viewerjumpto "Description" "catt_gt##description"}{...}
{viewerjumpto "Options" "catt_gt##options"}{...}
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
{synopt:{opt gteval(numlist)}}(g,t) evaluation pairs: {cmd:g1 t1 g2 t2 ...}{p_end}
{synopt:{opt porder(#)}}local polynomial order; 1 or 2; default is {cmd:2}{p_end}
{synopt:{opt kernel(string)}}kernel function; {cmd:gau} or {cmd:epa}; default is {cmd:gau}{p_end}
{synopt:{opt bwselect(string)}}bandwidth selection method; default is {cmd:IMSE1}{p_end}
{synopt:{opt bw(numlist)}}manual bandwidth; scalar or vector matching {cmd:gteval()} rows{p_end}

{syntab:Inference}
{synopt:{opt alp(#)}}significance level; default is {cmd:0.05}{p_end}
{synopt:{opt bstrap(true|false)}}explicit bootstrap toggle; default is {cmd:true}{p_end}
{synopt:{opt nobstrap}}legacy bootstrap-off flag{p_end}
{synopt:{opt biters(#)}}bootstrap iterations; default is {cmd:1000}{p_end}
{synopt:{opt seed(#)}}random-number seed for bootstrap reproducibility; default is {cmd:-1} (use current RNG state){p_end}
{synopt:{opt uniformall(true|false)}}explicit uniform-band toggle; default is {cmd:true}{p_end}
{synopt:{opt nouniformall}}legacy flag alias for {cmd:uniformall(false)}{p_end}

{syntab:Control group}
{synopt:{opt control_group(string)}}comparison group; {cmd:notyettreated} or {cmd:nevertreated}; default is {cmd:notyettreated}{p_end}
{synopt:{opt anticipation(#)}}number of anticipation periods; default is {cmd:0}{p_end}

{syntab:Pre-trends}
{synopt:{opt pretrend}}include pre-treatment periods for testing{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:catt_gt} estimates group-time conditional average treatment effects on the
treated, CATT(g,t,z), as a function of a continuous pre-treatment covariate
{it:z}. The estimator is doubly robust, combining inverse probability weighting
with outcome regression using local polynomial smoothing.

{pstd}
The method is designed for staggered difference-in-differences settings where
treatment effects may vary with a continuous pre-treatment characteristic. It
produces nonparametric estimates of the CATT function at user-specified
evaluation points {opt zeval()}, along with pointwise and uniform confidence
bands.

{pstd}
This command implements the estimation procedure described in Section 4 of
Imai, Qin, and Yanagi (2025).

{pstd}
After estimation, use {helpb aggte_gt} to aggregate results into summary
parameters and {helpb catt_gt_graph} to visualize the estimated CATT functions.


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
{opt gteval(numlist)} explicitly specifies the (g,t) pairs to evaluate, in
the format {cmd:g1 t1 g2 t2 ...}. When omitted, {cmd:catt_gt} automatically
constructs all identifiable (g,t) pairs based on the sample support,
{opt control_group()}, and {opt anticipation()}.

{phang}
{opt porder(#)} sets the order of the local polynomial regression. Must be 1
(local linear) or 2 (local quadratic). Default is {cmd:2}.

{phang}
{opt kernel(string)} specifies the kernel function for local polynomial
smoothing. Options are {cmd:gau} (Gaussian) and {cmd:epa} (Epanechnikov).
Default is {cmd:gau}.

{phang}
{opt bwselect(string)} specifies the bandwidth selection method. Options are
{cmd:IMSE1}, {cmd:IMSE2}, {cmd:US1}, and {cmd:manual}. Default is {cmd:IMSE1}.
When {cmd:manual} is specified, the {opt bw()} option must also be provided.

{phang}
{opt bw(numlist)} specifies manual bandwidth value(s). Only valid when
{cmd:bwselect(manual)} is specified. When {opt gteval()} is omitted,
{opt bw()} must be a single positive scalar. When {opt gteval()} is
explicitly provided, {opt bw()} may be either a single positive scalar or a
positive vector whose length equals the number of rows of {cmd:e(gteval)}.
If {opt gteval()} contains only one (g,t) pair, the effective inference
domain degenerates to uniform inference over {it:z} alone.

{dlgtab:Inference}

{phang}
{opt alp(#)} sets the significance level for confidence intervals. Must be
strictly between 0 and 1. Default is {cmd:0.05}.

{phang}
{opt bstrap(true|false)} explicitly controls the multiplier bootstrap for
inference. Bootstrap is enabled by default, so omitting the option keeps
bootstrap on, {cmd:bstrap(true)} turns it on explicitly, and
{cmd:bstrap(false)} turns it off. The legacy {opt nobstrap} flag remains
available for backward compatibility.

{phang}
{opt biters(#)} sets the number of bootstrap iterations. Default is {cmd:1000}.
Must be a positive integer.

{phang}
{opt seed(#)} sets the random-number seed immediately before bootstrap
resampling. The default {cmd:seed(-1)} leaves the current Stata RNG state
unchanged and therefore uses the seed already active in the session.

{phang}
{opt uniformall(true|false)} explicitly controls whether inference is joint
over all evaluation points and (g,t) pairs. Omitting the option keeps the
default {cmd:uniformall(true)} behavior. The legacy flag {opt nouniformall}
remains available as an alias for {cmd:uniformall(false)}.

{phang}
{opt nouniformall} is a legacy flag alias for {cmd:uniformall(false)}. Use
either form when you want confidence bands that are only uniform over {it:z}
within each (g,t) pair instead of a single joint critical value over (g,t,z).

{dlgtab:Control group}

{phang}
{opt control_group(string)} specifies the comparison group. {cmd:notyettreated}
uses units not yet treated by period {it:t} as controls. {cmd:nevertreated}
uses only never-treated units (group = 0). Default is {cmd:notyettreated}.

{phang}
{opt anticipation(#)} specifies the number of periods before treatment in which
anticipation effects may occur. Default is {cmd:0}.

{dlgtab:Pre-trends}

{phang}
{opt pretrend} includes pre-treatment periods in the estimation for
pre-trends testing. When specified, CATT estimates for pre-treatment periods
are computed and can be visualized with {helpb catt_gt_graph}. This option is
only applicable when {opt gteval()} is omitted, because explicit
{opt gteval()} already fixes the evaluation pairs.
More generally, the excluded long-difference baseline period is
t = g - anticipation - 1; this reduces to event_time = -1 only when
{cmd:anticipation(0)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup: generate simulation data{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}

{pstd}Example 1: Basic estimation on the safe deterministic path{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false)}{p_end}
{phang2}{cmd:. matrix list e(results)}{p_end}
{phang2}{it:// bstrap(false) runs the deterministic analytical-CI-only path}{p_end}

{pstd}Example 1b: Deterministic inference with explicit boolean toggles{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.5 0 0.5) ///}{p_end}
{phang2}{cmd:        gteval(2 2) xformula(Z) porder(1) kernel(gau) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.45) bstrap(false) uniformall(false)}{p_end}

{pstd}Example 2: Explicitly pass the default covariate set via {opt xformula()}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) xformula(Z)}{p_end}
{phang2}{it:// Note: omitting xformula() gives X=[1,Z]; xformula(Z) also gives X=[1,Z]}{p_end}

{pstd}Example 2b: Interaction syntax in {opt xformula()}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.5 0 0.5) gteval(2 2 2 3 3 3) ///}{p_end}
{phang2}{cmd:        xformula("Z*X1") bwselect(manual) bw(0.45 0.45 0.45) ///}{p_end}
{phang2}{cmd:        bstrap(false) uniformall(false)}{p_end}
{phang2}{it:// Note: xformula("Z*X1") is equivalent to xformula("Z + X1 + Z:X1")} {p_end}
{phang2}{it:// To include only the interaction, write xformula("Z + Z:X1")} {p_end}

{pstd}Example 3: Estimation and visualization without bootstrap{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 4: Pre-trends testing{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) pretrend}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 5: Estimate only selected {it:(g,t)} pairs{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) ///}{p_end}
{phang2}{cmd:        gteval(2 2 3 3) bwselect(manual) bw(0.30 0.45) uniformall(false)}{p_end}

{pstd}
For persistence, {cmd:catt_gt} posts a storable {cmd:e(b)} shell but does not
publish a covariance matrix. Generic covariance-based postestimation such as
{cmd:lincom}, {cmd:test}, and {cmd:testparm} is therefore not a supported
inferential workflow. Use {cmd:e(results)}, {cmd:e(catt_se)}, and the reported
confidence bands for inference.


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
{synopt:{cmd:e(results)}}main results matrix with columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw; rows ordered by ascending z within each (g,t){p_end}
{synopt:{cmd:e(Estimate)}}alias of {cmd:e(results)}{p_end}
{synopt:{cmd:e(Estimate_b)}}vectorized point estimates used by the storable eclass shell{p_end}
{synopt:{cmd:e(B_g_t)}}flattened influence function matrix B_{i,g,t}(z_r); used internally by {helpb aggte_gt}{p_end}
{synopt:{cmd:e(G_g)}}group indicator matrix (n x K); G_{i,g} = 1 if unit i belongs to group g{p_end}
{synopt:{cmd:e(mu_G_g)}}conditional mean of G at evaluation points (R x K); E[G_g | Z = z_r]{p_end}
{synopt:{cmd:e(gteval)}}evaluated (g,t) pairs{p_end}
{synopt:{cmd:e(zeval)}}evaluation points for z, sorted in ascending order{p_end}
{synopt:{cmd:e(bw)}}bandwidth per (g,t) pair{p_end}
{synopt:{cmd:e(c_hat)}}analytical joint critical values per (g,t) pair{p_end}
{synopt:{cmd:e(c_check)}}bootstrap critical values per (g,t) pair (if {opt bstrap}){p_end}
{synopt:{cmd:e(catt_est)}}CATT point estimates{p_end}
{synopt:{cmd:e(catt_se)}}CATT standard errors{p_end}
{synopt:{cmd:e(kd0_Z)}}kernel density estimates at z{p_end}
{synopt:{cmd:e(kd1_Z)}}kernel density derivative estimates at z{p_end}
{synopt:{cmd:e(Z_supp)}}support of Z{p_end}

{pstd}
The following additional matrices capture the panel snapshot used by
{helpb aggte_gt} for two-pass aggregation and are not part of the user-facing
inference output:

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(Z)}}row vector of unit-level {it:z} values used by the estimation{p_end}
{synopt:{cmd:e(dh_Y_wide)}}wide-format outcome matrix (n x T){p_end}
{synopt:{cmd:e(dh_G_unit)}}row vector of unit-level treatment group assignments{p_end}
{synopt:{cmd:e(dh_t_vals)}}row vector of distinct time-period values used by the estimation{p_end}
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
