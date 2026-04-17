{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "aggte_gt##syntax"}{...}
{viewerjumpto "Description" "aggte_gt##description"}{...}
{viewerjumpto "Options" "aggte_gt##options"}{...}
{viewerjumpto "Remarks" "aggte_gt##remarks"}{...}
{viewerjumpto "Examples" "aggte_gt##examples"}{...}
{viewerjumpto "Stored results" "aggte_gt##stored"}{...}
{viewerjumpto "References" "aggte_gt##references"}{...}
{viewerjumpto "Authors" "aggte_gt##authors"}{...}
{viewerjumpto "Also see" "aggte_gt##alsosee"}{...}
{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{cmd:aggte_gt} {hline 2}}Aggregate group-time CATT into summary treatment effect parameters{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:aggte_gt}{cmd:,}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Aggregation}
{synopt:{opt type(string)}}aggregation type; {cmd:dynamic}, {cmd:group}, {cmd:calendar}, or {cmd:simple}; default is {cmd:dynamic}{p_end}
{synopt:{opt eval(numlist)}}evaluation points along the aggregation dimension; not allowed with {cmd:type(simple)}{p_end}

{syntab:Estimation}
{synopt:{opt porder(#)}}local polynomial order; {cmd:1} or {cmd:2}; default is {cmd:2}{p_end}
{synopt:{opt bwselect(string)}}bandwidth selection method; default is {cmd:IMSE1}{p_end}
{synopt:{opt bw(numlist)}}manual bandwidth; requires {cmd:bwselect(manual)}{p_end}

{syntab:Inference}
{synopt:{opt bstrap(true|false)}}multiplier bootstrap toggle; default is {cmd:true}{p_end}
{synopt:{opt biters(#)}}bootstrap iterations; default is {cmd:1000}{p_end}
{synopt:{opt seed(#)}}RNG seed for bootstrap reproducibility; default is {cmd:-1} (current RNG state){p_end}
{synopt:{opt uniformall(true|false)}}joint uniform-band toggle; default is {cmd:true}{p_end}
{synoptline}

{pstd}
{cmd:aggte_gt} is a post-estimation command. It must be preceded by a
successful {helpb catt_gt} or {helpb didhetero} run; both commands produce
the {cmd:e()} matrices that {cmd:aggte_gt} consumes for aggregation.


{marker description}{...}
{title:Description}

{pstd}
{cmd:aggte_gt} aggregates the group-time CATT(g,t,z) estimates produced by
{helpb catt_gt} or {helpb didhetero} into interpretable summary parameters
indexed by the continuous covariate {it:z}. The command implements Section
5 of Imai, Qin, and Yanagi (2025). Four aggregation types are available:

{phang2}
{cmd:dynamic} {hline 1} event-study aggregation by relative time since
treatment onset (event time {it:e}). Here {it:e} is the elapsed period
count on the ordered panel support, not the raw difference between the
original time labels. Produces one CATT curve per event time.

{phang2}
{cmd:group} {hline 1} cohort aggregation by treatment group {it:g}.
Produces one CATT curve per cohort.

{phang2}
{cmd:calendar} {hline 1} calendar-time aggregation by period {it:t}.
Produces one CATT curve per calendar period.

{phang2}
{cmd:simple} {hline 1} a single weighted average across all
post-treatment (g,t) pairs. Produces one CATT curve.

{pstd}
Aggregation proceeds in two passes: an aggregation-specific bandwidth is
computed for each eval point, then the CATT surface is re-smoothed at
those bandwidths before the weighted average is taken. The kernel
function and significance level are inherited from the upstream
{helpb catt_gt}/{helpb didhetero} result to match the closed algorithm
path described in the paper. To change the kernel or the significance
level, set {opt kernel()} or {opt alp()} when running the upstream
estimator.


{marker options}{...}
{title:Options}

{dlgtab:Aggregation}

{phang}
{opt type(string)} specifies the aggregation type. Must be one of {cmd:dynamic},
{cmd:group}, {cmd:calendar}, or {cmd:simple}. Default is {cmd:dynamic}.

{phang}
{opt eval(numlist)} specifies the evaluation points along the aggregation
dimension. For {cmd:dynamic}, these are event times measured as elapsed period
counts on the ordered panel support (e.g., {cmd:eval(0 1 2)}), not raw
differences between time labels. For {cmd:group}, these are group values. For
{cmd:calendar}, these are calendar periods. If omitted, all available values
are used. {cmd:type(simple)} has no eval dimension, so {cmd:eval()} must be
omitted.

{dlgtab:Estimation}

{phang}
{opt porder(#)} sets the local polynomial order for the aggregation step.
Must be 1 or 2. Default is {cmd:2}.

{phang}
{opt bwselect(string)} specifies the bandwidth selection method for the
aggregation step. Options are {cmd:IMSE1}, {cmd:IMSE2}, {cmd:US1}, and
{cmd:manual}. Default is {cmd:IMSE1}.

{phang}
{opt bw(numlist)} specifies manual bandwidth value(s) and is valid only
when {cmd:bwselect(manual)}. All values must be positive. {opt bw()} must
be either a positive scalar or a positive vector whose length equals the
number of final {cmd:eval()} points used by {cmd:aggte_gt}. This rule
applies whether {cmd:eval()} is user-supplied or constructed automatically
from the upstream (g,t) support. When {opt uniformall(true)} is in
effect, {cmd:aggte_gt} further collapses the per-eval bandwidths to their
common minimum before the joint critical value is computed.

{dlgtab:Inference}

{phang}
{opt bstrap(true|false)} toggles the multiplier bootstrap for uniform
inference. Bootstrap is enabled by default.

{phang}
{opt biters(#)} sets the number of bootstrap iterations. Default is
{cmd:1000}. Must be a positive integer when the bootstrap is enabled.

{phang}
{opt seed(#)} sets the Stata RNG seed immediately before the multiplier
bootstrap weights are drawn. The default {cmd:seed(-1)} leaves the
current RNG state unchanged; the effective seed is recorded in
{cmd:e(aggte_seed_request)} and {cmd:e(aggte_seed)}.

{phang}
{opt uniformall(true|false)} toggles joint uniform inference over the
aggregation dimension and {it:z}. The default is {cmd:true}. When the
final aggregation uses only one {cmd:eval} point, the effective problem
degenerates to z-only uniform inference and {cmd:e(aggte_uniformall)} is
reported as {cmd:0} even if the user typed {cmd:uniformall(true)}.


{marker remarks}{...}
{title:Remarks}

{pstd}
The summary parameters produced by {cmd:aggte_gt} follow Section 5 of Imai,
Qin, and Yanagi (2025). Let {it:θ} denote the chosen summary parameter.
With weights {it:w(g,t)} determined by the aggregation type, the parameter is

{phang2}
{cmd:dynamic}: {it:θ_es(e, z) = Σ_g w(g,e) CATT(g, g+e, z)}{p_end}
{phang2}
{cmd:group}: {it:θ_grp(g, z) = Σ_t w(g,t) CATT(g, t, z)}{p_end}
{phang2}
{cmd:calendar}: {it:θ_cal(t, z) = Σ_g w(g,t) CATT(g, t, z)}{p_end}
{phang2}
{cmd:simple}: {it:θ_smp(z) = Σ_{g,t} w(g,t) CATT(g, t, z)}{p_end}

{pstd}
The command does not re-estimate the first-stage generalized propensity
score or outcome regression. It reuses the influence-function panel
snapshot stored in {cmd:e(B_g_t)} and {cmd:e(G_g)} by the upstream
{helpb catt_gt}/{helpb didhetero} call. Running {cmd:aggte_gt} multiple
times with different {opt type()} or {opt eval()} values therefore does
not require re-estimating the CATT surface.

{pstd}
After {cmd:aggte_gt}, the aggregated results live in {cmd:e(Estimate)}
(9 columns: {it:eval, z, est, se, ci1_lower, ci1_upper, ci2_lower,
ci2_upper, bw}), while the upstream CATT matrix remains in
{cmd:e(results)} (10 columns). Both surfaces are accessible to
{helpb catt_gt_graph}; pass {cmd:plot_type(Aggregated)} to plot the
summary parameter or {cmd:plot_type(CATT)} to plot the underlying CATT
object. For persistence, {cmd:aggte_gt} posts {cmd:e(b)} as the
aggregated point-estimate row and {cmd:e(V)} as a diagonal covariance
matrix formed from {cmd:e(aggte_se)^2}; the off-diagonal elements are
zero because the package's authoritative inference is the uniform
confidence bands in {cmd:e(Estimate)}.


{marker examples}{...}
{title:Examples}

{pstd}Setup: upstream CATT estimates on simulation data.{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.5 0 0.5) gteval(2 2 2 3 3 3) ///}{p_end}
{phang2}{cmd:        xformula(Z) porder(1) kernel(gau) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.45 0.45 0.45) ///}{p_end}
{phang2}{cmd:        bstrap(false) uniformall(false)}{p_end}

{pstd}Example 1: Event-study aggregation (default {cmd:type}).{p_end}
{phang2}{cmd:. aggte_gt, bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 2: Dynamic aggregation at selected event times with bootstrap.{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic) eval(0 1) biters(500) seed(42)}{p_end}

{pstd}Example 3: Cohort (group) aggregation.{p_end}
{phang2}{cmd:. aggte_gt, type(group) bstrap(false)}{p_end}
{phang2}{cmd:. matrix list e(Estimate)}{p_end}

{pstd}Example 4: Calendar-time aggregation.{p_end}
{phang2}{cmd:. aggte_gt, type(calendar) bstrap(false)}{p_end}

{pstd}Example 5: Simple weighted average across all post-treatment (g,t).{p_end}
{phang2}{cmd:. aggte_gt, type(simple) bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph, plot_type(Aggregated)}{p_end}

{pstd}Example 6: Manual bandwidth vector matching the eval dimension.{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic) eval(0 1) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.40 0.45) bstrap(false) uniformall(false)}{p_end}


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:aggte_gt} exposes three layers of stored results: the current
aggregated state (unprefixed and {cmd:aggte_*} prefixed, which carry the
same information), the upstream {helpb catt_gt}/{helpb didhetero} snapshot
({cmd:aggte_base_*}), and the panel snapshot preserved for re-running
{helpb catt_gt_graph} against the underlying CATT surface.

{pstd}
{cmd:aggte_gt} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars (current aggregation state)}{p_end}
{synopt:{cmd:e(N)}}number of cross-sectional units{p_end}
{synopt:{cmd:e(T)}}number of time periods{p_end}
{synopt:{cmd:e(gbar)}}upper bound for treatment-time support; {cmd:.} encodes {cmd:+Inf} when never-treated units exist{p_end}
{synopt:{cmd:e(gbar_isinf)}}1 if {cmd:e(gbar)} encodes {cmd:+Inf}, 0 otherwise{p_end}
{synopt:{cmd:e(num_gteval)}}number of (g,t) pairs evaluated by the upstream result{p_end}
{synopt:{cmd:e(num_zeval)}}number of z evaluation points in the upstream result{p_end}
{synopt:{cmd:e(porder)}}aggregation polynomial order{p_end}
{synopt:{cmd:e(alp)}}significance level (inherited from the upstream result){p_end}
{synopt:{cmd:e(bstrap)}}1 if bootstrap was performed during the current aggregation, 0 otherwise{p_end}
{synopt:{cmd:e(biters)}}effective bootstrap iterations; 0 if bootstrap is disabled{p_end}
{synopt:{cmd:e(seed_request)}}requested bootstrap seed; {cmd:-1} means use the current RNG state{p_end}
{synopt:{cmd:e(seed)}}effective command-level bootstrap seed; missing if bootstrap is disabled or no seed was applied{p_end}
{synopt:{cmd:e(uniformall)}}1 if joint uniform bands over {it:eval} and {it:z} are effective, 0 otherwise{p_end}
{synopt:{cmd:e(anticipation)}}anticipation periods inherited from the upstream result{p_end}
{synopt:{cmd:e(anticip)}}anticipation periods (legacy alias){p_end}
{synopt:{cmd:e(pretrend)}}1 if the upstream result includes pretrend periods, 0 otherwise{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars (explicit {cmd:aggte_*} alias)}{p_end}
{synopt:{cmd:e(aggte_porder)}}duplicate of {cmd:e(porder)}{p_end}
{synopt:{cmd:e(aggte_alp)}}duplicate of {cmd:e(alp)}{p_end}
{synopt:{cmd:e(aggte_bstrap)}}duplicate of {cmd:e(bstrap)}{p_end}
{synopt:{cmd:e(aggte_biters)}}duplicate of {cmd:e(biters)}{p_end}
{synopt:{cmd:e(aggte_seed_request)}}duplicate of {cmd:e(seed_request)}{p_end}
{synopt:{cmd:e(aggte_seed)}}duplicate of {cmd:e(seed)}{p_end}
{synopt:{cmd:e(aggte_uniformall)}}duplicate of {cmd:e(uniformall)}{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars (upstream snapshot)}{p_end}
{synopt:{cmd:e(aggte_base_porder)}}upstream {helpb catt_gt}/{helpb didhetero} polynomial order{p_end}
{synopt:{cmd:e(aggte_base_alp)}}upstream significance level{p_end}
{synopt:{cmd:e(aggte_base_bstrap)}}upstream bootstrap flag{p_end}
{synopt:{cmd:e(aggte_base_biters)}}upstream effective bootstrap iterations{p_end}
{synopt:{cmd:e(aggte_base_seed_request)}}upstream requested seed{p_end}
{synopt:{cmd:e(aggte_base_seed)}}upstream effective seed; missing if none was applied{p_end}
{synopt:{cmd:e(aggte_base_uniformall)}}upstream uniform-band flag{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:aggte_gt}{p_end}
{synopt:{cmd:e(type)}}aggregation type: {cmd:dynamic}, {cmd:group}, {cmd:calendar}, or {cmd:simple}{p_end}
{synopt:{cmd:e(aggte_type)}}duplicate of {cmd:e(type)}{p_end}
{synopt:{cmd:e(aggte_source_cmd)}}upstream command that produced the CATT object ({cmd:catt_gt} or {cmd:didhetero}){p_end}
{synopt:{cmd:e(kernel)}}kernel function used for the aggregation step{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection method used by {cmd:aggte_gt}{p_end}
{synopt:{cmd:e(aggte_kernel)}}duplicate of {cmd:e(kernel)}{p_end}
{synopt:{cmd:e(aggte_bwselect)}}duplicate of {cmd:e(bwselect)}{p_end}
{synopt:{cmd:e(aggte_base_kernel)}}upstream kernel function{p_end}
{synopt:{cmd:e(aggte_base_bwselect)}}upstream bandwidth selection method{p_end}
{synopt:{cmd:e(depvar)}}outcome variable inherited from the upstream result{p_end}
{synopt:{cmd:e(idvar)}}panel identifier inherited from the upstream result{p_end}
{synopt:{cmd:e(timevar)}}time variable inherited from the upstream result{p_end}
{synopt:{cmd:e(groupvar)}}group variable inherited from the upstream result{p_end}
{synopt:{cmd:e(zvar)}}continuous covariate variable inherited from the upstream result{p_end}
{synopt:{cmd:e(control_group)}}control-group rule inherited from the upstream result{p_end}
{synopt:{cmd:e(control)}}control-group rule (legacy alias){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices (aggregation output)}{p_end}
{synopt:{cmd:e(b)}}storable aggregated point-estimate row used by {cmd:estimates store}{p_end}
{synopt:{cmd:e(V)}}diagonal covariance matrix with {cmd:e(aggte_se)^2} entries; off-diagonals are zero{p_end}
{synopt:{cmd:e(Estimate)}}main results matrix (9 columns): {it:eval, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw}{p_end}
{synopt:{cmd:e(Estimate_b)}}vectorized aggregated point estimates used by {cmd:e(b)}{p_end}
{synopt:{cmd:e(aggte_est)}}aggregated point estimates{p_end}
{synopt:{cmd:e(aggte_se)}}aggregated standard errors{p_end}
{synopt:{cmd:e(aggte_ci1_lower)}}analytical CI lower bounds{p_end}
{synopt:{cmd:e(aggte_ci1_upper)}}analytical CI upper bounds{p_end}
{synopt:{cmd:e(aggte_ci2_lower)}}multiplier bootstrap CI lower bounds (missing when {cmd:e(aggte_bstrap)} is 0){p_end}
{synopt:{cmd:e(aggte_ci2_upper)}}multiplier bootstrap CI upper bounds (missing when {cmd:e(aggte_bstrap)} is 0){p_end}
{synopt:{cmd:e(aggte_bw)}}bandwidth per eval point used in the aggregation step{p_end}
{synopt:{cmd:e(aggte_eval)}}eval-dimension evaluation points{p_end}
{synopt:{cmd:e(aggte_zeval)}}z evaluation points used in the aggregation step{p_end}

{pstd}
The upstream CATT matrices produced by {helpb catt_gt}/{helpb didhetero}
are preserved in {cmd:e()} so that {helpb catt_gt_graph} can still plot
the original CATT surface after {cmd:aggte_gt} is run:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices (upstream snapshot)}{p_end}
{synopt:{cmd:e(results)}}upstream CATT results matrix (10 columns): {it:g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw}{p_end}
{synopt:{cmd:e(gteval)}}upstream (g,t) evaluation pairs{p_end}
{synopt:{cmd:e(zeval)}}upstream z evaluation points{p_end}
{synopt:{cmd:e(bw)}}upstream per-(g,t) bandwidth vector{p_end}
{synopt:{cmd:e(c_hat)}}upstream analytical critical values per (g,t) pair{p_end}
{synopt:{cmd:e(c_check)}}upstream bootstrap critical values (only when the upstream result was bootstrapped){p_end}
{synopt:{cmd:e(catt_est)}}upstream CATT point estimates{p_end}
{synopt:{cmd:e(catt_se)}}upstream CATT standard errors{p_end}


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
Help: {helpb didhetero}, {helpb catt_gt}, {helpb catt_gt_graph}, {helpb didhetero_simdata}
{p_end}
