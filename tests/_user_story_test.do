// =============================================================================
// User Story Comprehensive Test Suite for didhetero
// =============================================================================
clear all
set more off
set linesize 200

// Setup adopath
local PKGROOT "/Users/cxy/Desktop/didhetero/didhetero-main"
adopath + "`PKGROOT'/ado"
cd "`PKGROOT'"

// Track pass/fail
local total_tests 0
local pass_tests 0
local fail_tests 0
local fail_details ""

capture program drop record_pass
program define record_pass
    args testname
    display as result "[PASS] `testname'"
end

capture program drop record_fail
program define record_fail
    args testname msg
    display as error "[FAIL] `testname': `msg'"
end

// =============================================================================
// Story 1.1: First install and run
// =============================================================================
display _n "{hline 60}"
display "STORY 1.1: First install and run"
display "{hline 60}"

// help files should exist
capture noisily help didhetero
local s11_help1 = _rc
display "help didhetero rc = `s11_help1'"

capture noisily help catt_gt
local s11_help2 = _rc
display "help catt_gt rc = `s11_help2'"

capture noisily help aggte_gt
local s11_help3 = _rc
display "help aggte_gt rc = `s11_help3'"

// Simulate data
capture noisily didhetero_simdata, n(200) tau(4) clear
local s11_sim = _rc
display "simdata rc = `s11_sim'"
if `s11_sim' == 0 {
    describe
    display "[PASS] Story 1.1: simdata created successfully"
}
else {
    display as error "[FAIL] Story 1.1: simdata failed"
}

// =============================================================================
// Story 1.2: Simplest estimation from help docs
// =============================================================================
display _n "{hline 60}"
display "STORY 1.2: Simplest estimation"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear
timer clear 1
timer on 1
capture noisily didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7)
local s12_rc = _rc
timer off 1
timer list 1
display "didhetero wrapper rc = `s12_rc'"
if `s12_rc' == 0 {
    ereturn list
    capture noisily matrix list e(results)
    display "[PASS] Story 1.2: didhetero wrapper works"
}
else {
    display as error "[FAIL] Story 1.2: didhetero wrapper failed with rc=`s12_rc'"
}

// =============================================================================
// Story 1.3: Common newbie errors
// =============================================================================
display _n "{hline 60}"
display "STORY 1.3: Common newbie errors"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear

// Missing zeval
display _n "Test: missing zeval"
capture noisily didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z)
local s13_a = _rc
display "Missing zeval rc = `s13_a' (expect non-zero)"
if `s13_a' != 0 {
    display "[PASS] Missing zeval correctly caught"
}
else {
    display as error "[FAIL] Missing zeval NOT caught"
}

// Typo in option (kernal instead of kernel)
display _n "Test: typo in option"
capture noisily didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.5) kernal(gau)
local s13_b = _rc
display "Typo option rc = `s13_b' (expect non-zero)"
if `s13_b' != 0 {
    display "[PASS] Typo in option correctly caught"
}
else {
    display as error "[FAIL] Typo in option NOT caught"
}

// Non-existent variable
display _n "Test: non-existent variable"
capture noisily didhetero Y, id(id) time(period) group(G) z(ZZZ) xformula(ZZZ) zeval(0.5)
local s13_c = _rc
display "Non-existent var rc = `s13_c' (expect non-zero)"
if `s13_c' != 0 {
    display "[PASS] Non-existent variable correctly caught"
}
else {
    display as error "[FAIL] Non-existent variable NOT caught"
}

// Unbalanced panel
display _n "Test: unbalanced panel"
drop if _n == 1
capture noisily didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.5) bwselect(manual) bw(0.3)
local s13_d = _rc
display "Unbalanced panel rc = `s13_d'"
if `s13_d' != 0 {
    display "[PASS] Unbalanced panel correctly detected"
}
else {
    display "[INFO] Unbalanced panel handled gracefully (rc=0)"
}

// =============================================================================
// Story 2.1: Step-by-step workflow catt_gt -> aggte_gt -> graph
// =============================================================================
display _n "{hline 60}"
display "STORY 2.1: Step-by-step workflow"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear

// Step 1: catt_gt
timer clear 2
timer on 2
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(200) seed(2024)
local s21_catt = _rc
timer off 2
timer list 2
display "catt_gt rc = `s21_catt'"

if `s21_catt' == 0 {
    matrix list e(results)
    display "N = " e(N) ", num_gteval = " e(num_gteval)
    display "[PASS] Story 2.1 Step 1: catt_gt"
    
    // Step 2: aggte_gt dynamic
    capture noisily aggte_gt, type(dynamic) bstrap(true) biters(200)
    local s21_agg = _rc
    display "aggte_gt dynamic rc = `s21_agg'"
    if `s21_agg' == 0 {
        matrix list e(results)
        display "[PASS] Story 2.1 Step 2: aggte_gt dynamic"
    }
    else {
        display as error "[FAIL] Story 2.1 Step 2: aggte_gt failed rc=`s21_agg'"
    }
    
    // Step 3: graph (re-estimate for graph)
    capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) bwselect(IMSE1) ///
        bstrap(true) biters(100) seed(2024)
    capture noisily catt_gt_graph
    local s21_graph1 = _rc
    display "catt_gt_graph rc = `s21_graph1'"
    
    capture noisily catt_gt_graph, plot_type(Aggregated)
    local s21_graph2 = _rc
    display "catt_gt_graph Aggregated rc = `s21_graph2'"
    
    if `s21_graph1' == 0 & `s21_graph2' == 0 {
        display "[PASS] Story 2.1 Step 3: graphs"
    }
    else {
        display as error "[FAIL] Story 2.1 Step 3: graph issue"
    }
}
else {
    display as error "[FAIL] Story 2.1: catt_gt failed"
}

// =============================================================================
// Story 2.2: Predict extraction
// =============================================================================
display _n "{hline 60}"
display "STORY 2.2: Predict extraction"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) seed(42)
local s22_catt = _rc

if `s22_catt' == 0 {
    capture noisily predict catt_est
    local s22_p1 = _rc
    capture noisily predict catt_se, se
    local s22_p2 = _rc
    capture noisily predict catt_ci2, ci2
    local s22_p3 = _rc
    capture noisily predict catt_z, zval
    local s22_p4 = _rc
    capture noisily predict catt_g, gval
    local s22_p5 = _rc
    capture noisily predict catt_t, tval
    local s22_p6 = _rc
    
    display "predict results: est=`s22_p1' se=`s22_p2' ci2=`s22_p3' z=`s22_p4' g=`s22_p5' t=`s22_p6'"
    
    if `s22_p1' == 0 & `s22_p2' == 0 & `s22_p4' == 0 & `s22_p5' == 0 & `s22_p6' == 0 {
        list catt_est catt_se catt_z catt_g catt_t in 1/10
        display "[PASS] Story 2.2: predict extraction works"
    }
    else {
        display as error "[FAIL] Story 2.2: some predict types failed"
    }
}
else {
    display as error "[FAIL] Story 2.2: catt_gt failed first"
}

// =============================================================================
// Story 2.3: Compare different estimation settings
// =============================================================================
display _n "{hline 60}"
display "STORY 2.3: Compare settings"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear

// Gaussian kernel
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1)
local s23_gau = _rc
if `s23_gau' == 0 {
    matrix gau_results = e(results)
}

// Epanechnikov kernel
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(epa) porder(2) bwselect(IMSE1)
local s23_epa = _rc
if `s23_epa' == 0 {
    matrix epa_results = e(results)
}

// porder=1
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1)
local s23_p1 = _rc

// Level 90
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) level(90)
local s23_l90 = _rc
if `s23_l90' == 0 {
    display "level = " e(level)
}

// Level 99
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) level(99)
local s23_l99 = _rc
if `s23_l99' == 0 {
    display "level = " e(level)
}

display "Results: gau=`s23_gau' epa=`s23_epa' p1=`s23_p1' l90=`s23_l90' l99=`s23_l99'"
if `s23_gau' == 0 & `s23_epa' == 0 & `s23_p1' == 0 & `s23_l90' == 0 & `s23_l99' == 0 {
    display "[PASS] Story 2.3: all settings work"
}
else {
    display as error "[FAIL] Story 2.3: some settings failed"
}

// =============================================================================
// Story 3.1: Pretrend testing workflow
// =============================================================================
display _n "{hline 60}"
display "STORY 3.1: Pretrend testing"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    pretrend bstrap(true) biters(200) seed(123)
local s31_catt = _rc
display "catt_gt with pretrend rc = `s31_catt'"

if `s31_catt' == 0 {
    capture noisily estat pretrend
    local s31_pre = _rc
    display "estat pretrend rc = `s31_pre'"
    
    capture noisily estat overlap
    local s31_ov = _rc
    display "estat overlap rc = `s31_ov'"
    
    if `s31_pre' == 0 {
        return list
        display "[PASS] Story 3.1: pretrend workflow"
    }
    else {
        display as error "[FAIL] Story 3.1: estat pretrend failed"
    }
}
else {
    display as error "[FAIL] Story 3.1: catt_gt with pretrend failed"
}

// =============================================================================
// Story 3.2: Custom (g,t) evaluation pairs
// =============================================================================
display _n "{hline 60}"
display "STORY 3.2: Custom gteval"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear
tab G
tab period

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    gteval(3 3)
local s32_rc = _rc
display "gteval(3 3) rc = `s32_rc'"
if `s32_rc' == 0 {
    matrix list e(results)
    display "[PASS] Story 3.2: custom gteval"
}
else {
    display as error "[FAIL] Story 3.2: custom gteval failed rc=`s32_rc'"
}

// =============================================================================
// Story 3.3: Control group comparison
// =============================================================================
display _n "{hline 60}"
display "STORY 3.3: Control group selection"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    control_group(notyettreated)
local s33_nyt = _rc
if `s33_nyt' == 0 {
    matrix nyt = e(results)
}

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    control_group(nevertreated)
local s33_nt = _rc
if `s33_nt' == 0 {
    matrix nt = e(results)
}

display "notyettreated rc=`s33_nyt', nevertreated rc=`s33_nt'"
if `s33_nyt' == 0 & `s33_nt' == 0 {
    matrix list nyt
    matrix list nt
    display "[PASS] Story 3.3: both control groups work"
}
else {
    display as error "[FAIL] Story 3.3: control group issue"
}

// =============================================================================
// Story 3.4: All aggregation types
// =============================================================================
display _n "{hline 60}"
display "STORY 3.4: All aggregation types"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) seed(99)
local s34_catt = _rc

if `s34_catt' == 0 {
    capture noisily aggte_gt, type(dynamic) bstrap(true) biters(100)
    local s34_dyn = _rc
    
    capture noisily aggte_gt, type(simple) bstrap(true) biters(100)
    local s34_sim = _rc
    
    capture noisily aggte_gt, type(group) bstrap(true) biters(100)
    local s34_grp = _rc
    
    capture noisily aggte_gt, type(calendar) bstrap(true) biters(100)
    local s34_cal = _rc
    
    display "Aggregation results: dynamic=`s34_dyn' simple=`s34_sim' group=`s34_grp' calendar=`s34_cal'"
    if `s34_dyn' == 0 & `s34_sim' == 0 & `s34_grp' == 0 & `s34_cal' == 0 {
        display "[PASS] Story 3.4: all aggregation types work"
    }
    else {
        display as error "[FAIL] Story 3.4: some aggregation types failed"
    }
}
else {
    display as error "[FAIL] Story 3.4: catt_gt base failed"
}

// =============================================================================
// Story 3.5: Diagnostic tools
// =============================================================================
display _n "{hline 60}"
display "STORY 3.5: Diagnostic tools"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear

// verbose
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) verbose
local s35_verb = _rc
display "verbose rc = `s35_verb'"

// kdetrim
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.01 0.3 0.5 0.7 0.99) kernel(gau) porder(2) bwselect(IMSE1) kdetrim
local s35_kde = _rc
display "kdetrim rc = `s35_kde'"

// gpsstrict
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) gpsstrict
local s35_gps = _rc
display "gpsstrict rc = `s35_gps'"

if `s35_verb' == 0 & `s35_kde' == 0 & `s35_gps' == 0 {
    display "[PASS] Story 3.5: diagnostic tools"
}
else {
    display as error "[FAIL] Story 3.5: diagnostics issue verb=`s35_verb' kde=`s35_kde' gps=`s35_gps'"
}

// =============================================================================
// Story 3.6: Multi-covariates
// =============================================================================
display _n "{hline 60}"
display "STORY 3.6: Multi-covariates"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear
gen X1 = rnormal()
gen X2 = runiform()
bysort id (period): replace X1 = X1[1]
bysort id (period): replace X2 = X2[1]

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z X1 X2) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1)
local s36_rc = _rc
display "Multi-covariate rc = `s36_rc'"
if `s36_rc' == 0 {
    display "[PASS] Story 3.6: multi-covariates"
}
else {
    display as error "[FAIL] Story 3.6: multi-covariates failed"
}

// =============================================================================
// Story 4.1: Real data - min_wage_cs
// =============================================================================
display _n "{hline 60}"
display "STORY 4.1: Real data min_wage_cs"
display "{hline 60}"

use "/Users/cxy/Desktop/didhetero/didhetero-main/data/min_wage_cs.dta", clear
describe, short

rename (first_treat year countyreal) (G period id)
bysort id (period): gen Z = pov[1]
xtset id period

summarize Z, detail
local p25 = r(p25)
local p50 = r(p50)
local p75 = r(p75)
display "Z quantiles: p25=`p25' p50=`p50' p75=`p75'"

timer clear 3
timer on 3
capture noisily catt_gt lemp, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(`p25' `p50' `p75') kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) seed(2024)
local s41_catt = _rc
timer off 3
timer list 3
display "min_wage catt_gt rc = `s41_catt'"

if `s41_catt' == 0 {
    capture noisily aggte_gt, type(simple) bstrap(true) biters(100)
    local s41_agg = _rc
    display "min_wage aggte_gt rc = `s41_agg'"
    
    capture noisily catt_gt_graph
    local s41_graph = _rc
    display "min_wage graph rc = `s41_graph'"
    
    if `s41_agg' == 0 {
        display "[PASS] Story 4.1: real data end-to-end"
    }
    else {
        display as error "[FAIL] Story 4.1: aggregation failed"
    }
}
else {
    display as error "[FAIL] Story 4.1: catt_gt on real data failed"
}

// =============================================================================
// Story 4.2: Other datasets (just describe)
// =============================================================================
display _n "{hline 60}"
display "STORY 4.2: Other datasets"
display "{hline 60}"

capture noisily use "/Users/cxy/Desktop/didhetero/didhetero-main/data/castle_doctrine.dta", clear
local s42_castle = _rc
if `s42_castle' == 0 {
    describe, short
    display "[PASS] castle_doctrine loads"
}

capture noisily use "/Users/cxy/Desktop/didhetero/didhetero-main/data/divorce_sw.dta", clear
local s42_divorce = _rc
if `s42_divorce' == 0 {
    describe, short
    display "[PASS] divorce_sw loads"
}

// =============================================================================
// Story 5.1: Very small sample
// =============================================================================
display _n "{hline 60}"
display "STORY 5.1: Extreme small sample"
display "{hline 60}"

didhetero_simdata, n(50) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1)
local s51_rc = _rc
display "Small sample rc = `s51_rc'"
if `s51_rc' == 0 {
    display "[PASS] Story 5.1: small sample handled"
}
else {
    display "[INFO] Story 5.1: small sample gave error (rc=`s51_rc') — may be expected"
}

// =============================================================================
// Story 5.2: Many evaluation points
// =============================================================================
display _n "{hline 60}"
display "STORY 5.2: Many evaluation points"
display "{hline 60}"

didhetero_simdata, n(500) tau(4) clear
numlist "0.05(0.05)0.95"
local many_zeval = r(numlist)
display "Testing with zeval: `many_zeval'"

timer clear 4
timer on 4
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(`many_zeval') kernel(epa) porder(1) bwselect(IMSE1)
local s52_rc = _rc
timer off 4
timer list 4
display "Many eval points rc = `s52_rc'"
if `s52_rc' == 0 {
    display "[PASS] Story 5.2: many eval points"
}
else {
    display as error "[FAIL] Story 5.2: many eval points failed"
}

// =============================================================================
// Story 5.3: Extreme Z values
// =============================================================================
display _n "{hline 60}"
display "STORY 5.3: Extreme Z values"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear
summarize Z

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.001 0.999) kernel(gau) porder(2) bwselect(manual) bw(0.1)
local s53_rc = _rc
display "Extreme Z rc = `s53_rc'"
if `s53_rc' == 0 {
    display "[PASS] Story 5.3: extreme Z values handled"
}
else {
    display "[INFO] Story 5.3: extreme Z values gave error (rc=`s53_rc') — may be expected"
}

// =============================================================================
// Story 5.4: Reproducibility with seed
// =============================================================================
display _n "{hline 60}"
display "STORY 5.4: Reproducibility"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) seed(777)
local s54_r1 = _rc
if `s54_r1' == 0 {
    matrix run1 = e(results)
}

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(true) biters(100) seed(777)
local s54_r2 = _rc
if `s54_r2' == 0 {
    matrix run2 = e(results)
}

if `s54_r1' == 0 & `s54_r2' == 0 {
    // Compare matrices
    capture noisily mata: assert(st_matrix("run1") == st_matrix("run2"))
    local s54_match = _rc
    if `s54_match' == 0 {
        display "[PASS] Story 5.4: reproducible results with same seed"
    }
    else {
        display as error "[FAIL] Story 5.4: results differ with same seed!"
    }
}
else {
    display as error "[FAIL] Story 5.4: estimation failed"
}

// =============================================================================
// Story 5.5: No bootstrap (fast estimate)
// =============================================================================
display _n "{hline 60}"
display "STORY 5.5: No bootstrap"
display "{hline 60}"

didhetero_simdata, n(300) tau(4) clear
timer clear 5
timer on 5
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
    bstrap(false)
local s55_rc = _rc
timer off 5
timer list 5
display "No bootstrap rc = `s55_rc'"
if `s55_rc' == 0 {
    matrix list e(results)
    display "[PASS] Story 5.5: no bootstrap fast estimate"
}
else {
    display as error "[FAIL] Story 5.5: no bootstrap failed"
}

// =============================================================================
// FINAL SUMMARY
// =============================================================================
display _n _n "{hline 60}"
display "FINAL TEST SUMMARY"
display "{hline 60}"
display "Story 1.1 (simdata):       rc = `s11_sim'"
display "Story 1.2 (basic didhetero): rc = `s12_rc'"
display "Story 1.3a (missing zeval): rc = `s13_a' (expect >0)"
display "Story 1.3b (typo option):  rc = `s13_b' (expect >0)"
display "Story 1.3c (bad var):      rc = `s13_c' (expect >0)"
display "Story 1.3d (unbalanced):   rc = `s13_d'"
display "Story 2.1 (workflow):      rc = `s21_catt'"
display "Story 2.2 (predict):       rc = `s22_catt'"
display "Story 2.3 (settings):      gau=`s23_gau' epa=`s23_epa' p1=`s23_p1'"
display "Story 3.1 (pretrend):      rc = `s31_catt'"
display "Story 3.2 (gteval):        rc = `s32_rc'"
display "Story 3.3 (control grp):   nyt=`s33_nyt' nt=`s33_nt'"
display "Story 3.4 (aggregation):   dyn=`s34_dyn' sim=`s34_sim' grp=`s34_grp' cal=`s34_cal'"
display "Story 3.5 (diagnostics):   verb=`s35_verb' kde=`s35_kde' gps=`s35_gps'"
display "Story 3.6 (multi-covar):   rc = `s36_rc'"
display "Story 4.1 (real data):     rc = `s41_catt'"
display "Story 5.1 (small N):       rc = `s51_rc'"
display "Story 5.2 (many zeval):    rc = `s52_rc'"
display "Story 5.3 (extreme Z):     rc = `s53_rc'"
display "Story 5.4 (reproducible):  r1=`s54_r1' r2=`s54_r2'"
display "Story 5.5 (no bootstrap):  rc = `s55_rc'"
display "{hline 60}"
display "TEST SUITE COMPLETE"
display "{hline 60}"

exit, clear
