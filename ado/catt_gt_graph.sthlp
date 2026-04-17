{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "catt_gt_graph##syntax"}{...}
{viewerjumpto "Description" "catt_gt_graph##description"}{...}
{viewerjumpto "Options" "catt_gt_graph##options"}{...}
{viewerjumpto "Examples" "catt_gt_graph##examples"}{...}
{viewerjumpto "Remarks" "catt_gt_graph##remarks"}{...}
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
{cmd:catt_gt_graph} plots estimated CATT functions or aggregated treatment
effect curves from {helpb didhetero}, {helpb catt_gt}, or {helpb aggte_gt}
results stored in {cmd:e()}.

{pstd}
The command auto-detects the plot mode based on the column count of the
results matrix:

{phang2}
10 columns (from {helpb catt_gt} or {helpb didhetero}) {hline 1} {cmd:CATT}
mode. One plot is generated per (g,t) pair, showing CATT(g,t,z) as a function
of z with confidence bands.

{phang2}
9 columns (from {helpb aggte_gt}) {hline 1} {cmd:Aggregated} mode. One plot
is generated per evaluation point (or a single plot for {cmd:simple} type),
showing the aggregated CATT as a function of z.

{pstd}
Confidence bands are selected automatically: bootstrap-based uniform bands
(ci2) are used when available; otherwise analytical bands (ci1) are used.

{pstd}
If {opt plot_type(string)} is specified, {cmd:catt_gt_graph} first selects the
matching result object before checking dimensions:
{cmd:plot_type(CATT)} looks for the CATT matrix in {cmd:e(results)}, while
{cmd:plot_type(Aggregated)} looks for the aggregated matrix in {cmd:e(Estimate)}.
This allows users to revisit the original CATT plots even after running
{helpb aggte_gt}.

{pstd}
When pre-trends testing was requested ({cmd:pretrend} option in {helpb catt_gt}),
a dashed red zero reference line is added to all plots.


{marker options}{...}
{title:Options}

{phang}
{opt plot_type(string)} overrides the auto-detected plot mode. Must be
{cmd:CATT} (requires 10-column matrix) or {cmd:Aggregated} (requires 9-column
matrix). If omitted, the mode is determined automatically from the preferred
results matrix dimensions. If specified, the command searches the corresponding
stored result first ({cmd:e(results)} for {cmd:CATT}; {cmd:e(Estimate)} for
{cmd:Aggregated}).

{phang}
{opt save_path(string)} specifies a file path to save the last generated graph.
Supported formats include {cmd:.png}, {cmd:.pdf}, and {cmd:.gph}. Paths ending
in {cmd:.gph} are saved with Stata's {cmd:graph save}; other supported suffixes
are exported with {cmd:graph export}.

{phang}
{opt graph_opt(string)} passes additional options to the underlying Stata
{cmd:twoway} graph command. These are appended to the default graph options.
For example: {cmd:graph_opt(xlabel(, format(%4.1f)))}.


{marker examples}{...}
{title:Examples}

{pstd}Setup (safe non-bootstrap base){p_end}
{phang2}{cmd:. didhetero_simdata, n(500) tau(4) seed(12345) clear}{p_end}
{phang2}{cmd:. catt_gt Y, group(G) time(period) id(id) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) bstrap(false) uniformall(true)}{p_end}
{phang2}{it:// The setup above uses the deterministic (analytical-CI-only) base}{p_end}

{pstd}Example 1: Plot CATT estimates (one graph per (g,t) pair){p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 2: Plot CATT estimates directly after didhetero{p_end}
{phang2}{cmd:. didhetero Y, id(id) time(period) group(G) z(Z) zeval(-0.8 -0.4 0 0.4 0.8) gteval(2 2) xformula(Z)}{p_end}
{phang2}{cmd:. catt_gt_graph, plot_type(CATT)}{p_end}

{pstd}Example 3: Plot aggregated estimates after aggte_gt{p_end}
{phang2}{cmd:. aggte_gt, type(dynamic)}{p_end}
{phang2}{cmd:. catt_gt_graph}{p_end}

{pstd}Example 4: Save graph to file{p_end}
{phang2}{cmd:. catt_gt_graph, save_path(my_catt_plot.png)}{p_end}

{pstd}Example 5: Custom graph options{p_end}
{phang2}{cmd:. catt_gt_graph, graph_opt(scheme(s2color) xlabel(, format(%4.1f)))}{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
Default graph style settings:

{p2colset 9 30 32 2}{...}
{p2col:Element}Default{p_end}
{p2line}
{p2col:CI band color}gs12 at 60% opacity{p_end}
{p2col:Estimate line}black, very thick{p_end}
{p2col:Plot region}white background{p_end}
{p2col:Legend}off{p_end}
{p2col:Pre-trend line}red dashed at y=0{p_end}
{p2line}

{pstd}
Graph names follow the pattern {cmd:g{it:G}_t{it:T}} in CATT mode and
{cmd:eval{it:E}} in Aggregated mode (with decimals replaced by underscores
and negative signs replaced by {cmd:m}).


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
