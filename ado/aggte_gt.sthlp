{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "aggte_gt##syntax"}{...}
{viewerjumpto "Description" "aggte_gt##description"}{...}
{viewerjumpto "Options" "aggte_gt##options"}{...}
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
{synopt:{opt eval(numlist)}}evaluation points for the aggregation dimension; not allowed with {cmd:type(simple)}{p_end}

{syntab:Estimation}
{synopt:{opt porder(#)}}local polynomial order; 1 or 2; default is {cmd:2}{p_end}
{synopt:{opt bwselect(string)}}bandwidth selection method; default is {cmd:IMSE1}{p_end}
{synopt:{opt bw(numlist)}}manual bandwidth override; length must equal the final {cmd:eval()} count used by {cmd:aggte_gt}{p_end}

{syntab:Inference}
{synopt:{opt bstrap(string)}}{cmd:true} or {cmd:false}; enable bootstrap; default is {cmd:true}{p_end}
{synopt:{opt biters(#)}}bootstrap iterations; default is {cmd:1000}{p_end}
{synopt:{opt seed(#)}}random-number seed for bootstrap reproducibility; default is {cmd:-1} (use current RNG state){p_end}
{synopt:{opt uniformall(string)}}{cmd:true} or {cmd:false}; joint uniform bands; default is {cmd:true}{p_end}
{synoptline}

{pstd}
{cmd:aggte_gt} is a post-estimation command. It must be run after {helpb catt_gt}
or {helpb didhetero}. Both commands produce the required {cmd:e()} matrices that
{cmd:aggte_gt} consumes for aggregation.


{marker description}{...}
{title:Description}

{pstd}
{cmd:aggte_gt} aggregates the group-time CATT(g,t,z) estimates produced by
{helpb catt_gt} into interpretable summary parameters. Four aggregation types
are available:

{phang2}
{cmd:dynamic} — event-study aggregation. Summarizes treatment effects by
relative time since treatment onset (event time {it:e}). Here {it:e} is the
elapsed period count on the ordered panel support, not the raw difference
between the original time labels. Produces one CATT curve per event time.

{phang2}
{cmd:group} — group aggregation. Summarizes treatment effects by treatment
cohort (group {it:g}). Produces one CATT curve per group.

{phang2}
{cmd:calendar} — calendar-time aggregation. Summarizes treatment effects by
calendar period {it:t}. Produces one CATT curve per period.

{phang2}
{cmd:simple} — simple weighted average across all post-treatment (g,t) pairs.
Produces a single CATT curve.

{pstd}
The aggregation procedure involves a two-pass estimation: first, aggregation-
specific bandwidths are computed; then CATT estimates are re-evaluated at those
bandwidths before weighted averaging.

{pstd}
The aggregation step inherits the kernel function from the preceding
{helpb catt_gt} or {helpb didhetero} result, matching the closed algorithm path
described in Imai, Qin, and Yanagi (2025).

{pstd}
The aggregation step also inherits the significance level from the preceding
{helpb catt_gt} or {helpb didhetero} result via {cmd:e(alp)}. To change the
confidence level, set {opt alp()} when estimating the upstream CATT object,
then run {cmd:aggte_gt} on that result.


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
{opt bw(numlist)} specifies manual bandwidth(s), overriding automatic
selection. Values must be positive. When {opt bwselect(manual)} is used, the
vector whose length equals the number of final {cmd:eval()} points used by
{cmd:aggte_gt}. This rule also applies when {cmd:eval()} is omitted and the
command constructs multiple evaluation points automatically.

{p 8 12 2}When {opt uniformall(true)} is in effect, {cmd:aggte_gt} still
collapses the per-eval bandwidths to their common minimum for inference after
validating that the user-supplied {cmd:bw()} length matches the number of eval
points.{p_end}

{dlgtab:Inference}

{phang}
{opt bstrap(string)} enables or disables the multiplier bootstrap.
{cmd:true} (default) or {cmd:false}.

{phang}
{opt biters(#)} sets the number of bootstrap iterations. Default is {cmd:1000}.

{phang}
{opt seed(#)} sets the random-number seed immediately before multiplier
bootstrap weights are generated for aggregation inference. The default
{cmd:seed(-1)} leaves the current Stata RNG state unchanged.

{phang}
{opt uniformall(string)} enables or disables joint uniform confidence bands.
{cmd:true} (default) or {cmd:false}.

{p 8 12 2}If the final aggregation uses only one {cmd:eval} point, the
effective inference problem degenerates to z-only uniform inference. In that
case, {cmd:e(aggte_uniformall)} is reported as {cmd:0} even if the user typed
{cmd:uniformall(true)}.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.5 0 0.5) ///}{p_end}
{phang2}{cmd:        gteval(2 2 2 3 3 3) xformula(Z) porder(1) kernel(gau) ///}{p_end}
{phang2}{cmd:        bwselect(manual) bw(0.45 0.45 0.45) bstrap(false) uniformall(false)}{p_end}

{pstd}Example 1: Event-study aggregation (dynamic){p_end}
{phang2}{cmd:. aggte_gt, type(dynamic)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 1b: Reproducible aggregation bootstrap with an explicit seed{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic) eval(0 1) bstrap(true) biters(500) seed(42)}{p_end}

{pstd}Example 2: Group aggregation{p_end}
{phang2}{cmd:. aggte_gt, type(group)}{p_end}
{phang2}{cmd:. matrix list e(Estimate)}{p_end}

{pstd}Example 3: Simple weighted average{p_end}
{phang2}{cmd:. aggte_gt, type(simple)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:aggte_gt} stores aggregation-specific results in the unprefixed
{cmd:e()} namespace for the current post-estimation layer. When you need the
underlying {helpb catt_gt}/{helpb didhetero} state, use the documented
{cmd:aggte_base_*} snapshot entries:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(pretrend)}}1 if the upstream result includes pretrend periods, 0 otherwise{p_end}
{synopt:{cmd:e(porder)}}current aggregation polynomial order{p_end}
{synopt:{cmd:e(bstrap)}}1 if bootstrap is enabled for the current {cmd:aggte_gt} result, 0 otherwise{p_end}
{synopt:{cmd:e(uniformall)}}1 if joint uniform bands over multiple eval points are effective for the current {cmd:aggte_gt} result, 0 otherwise{p_end}
{synopt:{cmd:e(alp)}}significance level used by the current {cmd:aggte_gt} result{p_end}
{synopt:{cmd:e(biters)}}effective bootstrap iterations for the current {cmd:aggte_gt} result; 0 if bootstrap is disabled{p_end}
{synopt:{cmd:e(seed_request)}}requested bootstrap seed for the current {cmd:aggte_gt} result; {cmd:-1} means use the current RNG state if bootstrap runs{p_end}
{synopt:{cmd:e(seed)}}effective command-level bootstrap seed for the current {cmd:aggte_gt} result; missing if bootstrap is disabled or no command-level seed was applied{p_end}
{synopt:{cmd:e(aggte_porder)}}polynomial order used{p_end}
{synopt:{cmd:e(aggte_bstrap)}}1 if bootstrap is enabled, 0 otherwise{p_end}
{synopt:{cmd:e(aggte_uniformall)}}1 if joint uniform bands over multiple eval points are effective, 0 otherwise{p_end}
{synopt:{cmd:e(aggte_alp)}}significance level inherited from the upstream result{p_end}
{synopt:{cmd:e(aggte_biters)}}effective bootstrap iterations; 0 if bootstrap is disabled{p_end}
{synopt:{cmd:e(aggte_seed_request)}}requested bootstrap seed for {cmd:aggte_gt}; {cmd:-1} means use the current RNG state if bootstrap runs{p_end}
{synopt:{cmd:e(aggte_seed)}}effective command-level bootstrap seed for {cmd:aggte_gt}; missing if bootstrap is disabled or no command-level seed was applied{p_end}
{synopt:{cmd:e(aggte_base_seed_request)}}requested bootstrap seed inherited from the upstream {cmd:catt_gt}/{cmd:didhetero} result{p_end}
{synopt:{cmd:e(aggte_base_seed)}}effective command-level bootstrap seed inherited from the upstream result; missing if none was applied{p_end}
{synopt:{cmd:e(aggte_base_porder)}}upstream {helpb catt_gt}/{helpb didhetero} polynomial order snapshot{p_end}
{synopt:{cmd:e(aggte_base_bstrap)}}upstream bootstrap flag snapshot{p_end}
{synopt:{cmd:e(aggte_base_uniformall)}}upstream uniform-band flag snapshot{p_end}
{synopt:{cmd:e(aggte_base_alp)}}upstream significance level snapshot{p_end}
{synopt:{cmd:e(aggte_base_biters)}}upstream effective bootstrap iterations snapshot{p_end}
{synopt:{cmd:e(aggte_base_seed_request)}}upstream requested seed snapshot{p_end}
{synopt:{cmd:e(aggte_base_seed)}}upstream effective seed snapshot{p_end}
{synopt:{cmd:e(gbar_isinf)}}1 indicates the current {cmd:e(gbar)} encodes {cmd:+Inf}; 0 indicates a finite upper bound{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:aggte_gt}{p_end}
{synopt:{cmd:e(type)}}aggregation type{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name inherited from the upstream result{p_end}
{synopt:{cmd:e(idvar)}}panel identifier inherited from the upstream result{p_end}
{synopt:{cmd:e(timevar)}}time variable inherited from the upstream result{p_end}
{synopt:{cmd:e(groupvar)}}group variable inherited from the upstream result{p_end}
{synopt:{cmd:e(zvar)}}continuous pre-treatment covariate variable inherited from the upstream result{p_end}
{synopt:{cmd:e(kernel)}}kernel function used by the current {cmd:aggte_gt} result{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection method used by the current {cmd:aggte_gt} result{p_end}
{synopt:{cmd:e(control_group)}}control-group rule inherited from the upstream result{p_end}
{synopt:{cmd:e(control)}}control-group alias inherited from the upstream result{p_end}
{synopt:{cmd:e(aggte_type)}}aggregation type (alias){p_end}
{synopt:{cmd:e(aggte_source_cmd)}}upstream command that produced the aggregated result object{p_end}
{synopt:{cmd:e(aggte_kernel)}}kernel function used{p_end}
{synopt:{cmd:e(aggte_bwselect)}}bandwidth selection method used{p_end}
{synopt:{cmd:e(aggte_base_kernel)}}upstream kernel snapshot{p_end}
{synopt:{cmd:e(aggte_base_bwselect)}}upstream bandwidth-selection snapshot{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(Estimate)}}main results matrix (9 columns): eval, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw{p_end}
{synopt:{cmd:e(aggte_est)}}aggregated point estimates{p_end}
{synopt:{cmd:e(aggte_se)}}aggregated standard errors{p_end}
{synopt:{cmd:e(aggte_ci1_lower)}}analytical CI lower bounds{p_end}
{synopt:{cmd:e(aggte_ci1_upper)}}analytical CI upper bounds{p_end}
{synopt:{cmd:e(aggte_ci2_lower)}}bootstrap CI lower bounds{p_end}
{synopt:{cmd:e(aggte_ci2_upper)}}bootstrap CI upper bounds{p_end}
{synopt:{cmd:e(aggte_bw)}}bandwidths used{p_end}
{synopt:{cmd:e(aggte_eval)}}evaluation points{p_end}
{synopt:{cmd:e(aggte_zeval)}}z evaluation points{p_end}


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
