{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "catt_gt_graph##syntax"}{...}
{viewerjumpto "Description" "catt_gt_graph##description"}{...}
{viewerjumpto "Options" "catt_gt_graph##options"}{...}
{viewerjumpto "Examples" "catt_gt_graph##examples"}{...}
{viewerjumpto "Remarks" "catt_gt_graph##remarks"}{...}
{viewerjumpto "Stored results" "catt_gt_graph##stored"}{...}
{viewerjumpto "References" "catt_gt_graph##references"}{...}
{viewerjumpto "Authors" "catt_gt_graph##authors"}{...}
{viewerjumpto "Also see" "catt_gt_graph##alsosee"}{...}
{title:Title}

{p2colset 5 25 27 2}{...}
{p2col:{cmd:catt_gt_graph} {hline 2}}Visualize CATT or aggregated treatment effect estimates{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:catt_gt_graph}{cmd:,}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt plot_type(string)}}force plot mode; {cmd:CATT} or {cmd:Aggregated}{p_end}
{synopt:{opt save_path(string)}}file path to save the last graph{p_end}
{synopt:{opt graph_opt(string)}}additional Stata {cmd:twoway} graph options{p_end}
{synoptline}

{pstd}
{cmd:catt_gt_graph} is a post-estimation command. It must be run after
{helpb didhetero}, {helpb catt_gt}, or {helpb aggte_gt}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:catt_gt_graph} plots estimated CATT functions or aggregated
treatment-effect curves from {helpb didhetero}, {helpb catt_gt}, or
{helpb aggte_gt} results stored in {cmd:e()}. The command produces one
{helpb twoway} graph per (g,t) pair (or per eval point in aggregated mode)
showing the point estimate and its confidence band as functions of the
continuous covariate {it:z}.

{pstd}
The plot mode is determined automatically from the column count of the
preferred results matrix: 10 columns ({helpb catt_gt}/{helpb didhetero})
triggers {cmd:CATT} mode, while 9 columns ({helpb aggte_gt}) triggers
{cmd:Aggregated} mode. When {opt plot_type()} is supplied, it overrides
auto-detection: {cmd:plot_type(CATT)} pulls the 10-column matrix from
{cmd:e(results)} (useful after {cmd:aggte_gt}, which preserves the
upstream CATT object), while {cmd:plot_type(Aggregated)} pulls the
9-column matrix from {cmd:e(Estimate)}.

{pstd}
Confidence bands are selected panel by panel. Bootstrap uniform bands
(the {cmd:ci2} columns) are used whenever they are non-missing for every
point on the panel; otherwise the analytical bands ({cmd:ci1}) are
used. When pre-trends testing was requested via {opt pretrend} on the
upstream estimator, a red dashed zero reference line is added to every
plot.


{marker options}{...}
{title:Options}

{phang}
{opt plot_type(string)} overrides the auto-detected plot mode. Accepts
{cmd:CATT} or {cmd:Aggregated}; quoted forms such as {cmd:plot_type("CATT")}
are also accepted. {cmd:CATT} requires a 10-column matrix and reads from
{cmd:e(results)}; {cmd:Aggregated} requires a 9-column matrix and reads
from {cmd:e(Estimate)}. When omitted, the mode is inferred from the
preferred results matrix dimensions.

{phang}
{opt save_path(string)} specifies a file path for the last generated
graph. Paths ending in {cmd:.gph} are saved with {helpb graph save};
other suffixes such as {cmd:.png} and {cmd:.pdf} are exported with
{helpb graph export}. The full path (absolute or relative) is honored
verbatim.

{phang}
{opt graph_opt(string)} passes additional options to the underlying
{helpb twoway} command. Options are appended to the defaults listed in
the Remarks section. For example,
{cmd:graph_opt(xlabel(, format(%4.1f)))} sets a custom x-axis format.


{marker examples}{...}
{title:Examples}

{pstd}Setup: deterministic (analytical-CI-only) CATT base.{p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false)}{p_end}

{pstd}Example 1: Auto-detected CATT plot (one graph per (g,t) pair).{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 2: Force CATT mode directly after {helpb didhetero}.{p_end}
{phang2}{cmd:. didhetero Y, id(id) time(period) group(G) z(Z) ///}{p_end}
{phang2}{cmd:        zeval(-0.8 -0.4 0 0.4 0.8) gteval(2 2) xformula(Z) ///}{p_end}
{phang2}{cmd:        bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph, plot_type(CATT)}{p_end}

{pstd}Example 3: Plot aggregated estimates after {helpb aggte_gt}.{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic) bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph, plot_type(Aggregated)}{p_end}

{pstd}Example 4: Recover the original CATT plot after aggregation.{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic) bstrap(false)}{p_end}
{phang2}{cmd:. catt_gt_graph, plot_type(CATT)}{p_end}

{pstd}Example 5: Save the last graph to a file.{p_end}
{phang2}{cmd:. catt_gt_graph, save_path(my_catt_plot.png)}{p_end}

{pstd}Example 6: Custom twoway options.{p_end}
{phang2}{cmd:. catt_gt_graph, graph_opt(scheme(s2color) xlabel(, format(%4.1f)))}{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
Default graph style settings used by {cmd:catt_gt_graph}:

{p2colset 9 30 32 2}{...}
{p2col:Element}Default{p_end}
{p2line}
{p2col:CI band color}{cmd:gs12} at 60% opacity{p_end}
{p2col:Estimate line}black, very thick{p_end}
{p2col:Plot region}white background{p_end}
{p2col:Legend}off{p_end}
{p2col:Pre-trend line}red dashed at {it:y} = 0{p_end}
{p2line}

{pstd}
The graph name pattern is {cmd:g{it:G}_t{it:T}} in CATT mode and
{cmd:eval{it:E}} in Aggregated mode, with decimal points replaced by
underscores and negative signs replaced by the letter {cmd:m}. For
example, CATT at {it:g} = 2, {it:t} = 3 becomes graph
{cmd:g2_t3}; an aggregated plot at eval = -0.5 becomes {cmd:evalm0_5}.
These names are reusable with {helpb graph use} or {helpb graph combine}
after the command returns.

{pstd}
Confidence-band selection is panel-specific. For each (g,t) pair in
CATT mode or each eval point in Aggregated mode, {cmd:catt_gt_graph}
checks whether the {cmd:ci2} columns are non-missing across the entire
panel. If so, the bootstrap uniform band is plotted; otherwise the
analytical band is used. A summary line reports the choice for each
panel before plotting.

{pstd}
{cmd:catt_gt_graph} does not modify {cmd:e()} and can therefore be
called multiple times in a row with different {opt plot_type()} or
{opt graph_opt()} values on the same estimation result.


{marker stored}{...}
{title:Stored results}

{pstd}
{cmd:catt_gt_graph} is not an estimation command and does not post results
to {cmd:e()}. It reads from the existing {cmd:e()} contents produced by
{helpb catt_gt}, {helpb didhetero}, or {helpb aggte_gt} and produces
graphs only.


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
Help: {helpb didhetero}, {helpb catt_gt}, {helpb aggte_gt}, {helpb didhetero_simdata}
{p_end}
