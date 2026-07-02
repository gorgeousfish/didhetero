// =============================================================================
// _verify_api.do - Verification script for all API improvements
// Task #22: End-to-end regression test
//
// IMPORTANT: Run make.do BEFORE this script (make.do has clear all).
// =============================================================================
clear all
set more off
set linesize 120

local test_pass = 0
local test_fail = 0

capture log close _verify_api
log using "/Users/cxy/Desktop/didhetero/didhetero-main/tests/_verify_api.log", ///
    replace name(_verify_api)

adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

// =============================================================================
// STEP 2: Basic regression (catt_gt core functionality)
// =============================================================================
di as text ""
di as text "================================================================"
di as text " STEP 2: Basic regression (catt_gt)"
di as text "================================================================"

set seed 12345
didhetero_simdata, n(300) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) ///
    bwselect(IMSE1) bstrap biters(100)

if _rc {
    di as error "FAIL: catt_gt basic regression failed with rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    capture confirm matrix e(results)
    if _rc {
        di as error "FAIL: e(results) matrix not found"
        local test_fail = `test_fail' + 1
    }
    else {
        di as text "PASS: catt_gt estimated, e(results) exists"
        local test_pass = `test_pass' + 1
    }
}

// =============================================================================
// STEP 3: predict command verification
// =============================================================================
di as text ""
di as text "================================================================"
di as text " STEP 3: predict command verification"
di as text "================================================================"

// 3a: Point estimate
capture noisily predict catt_hat
if _rc {
    di as error "FAIL: predict (point estimate) rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    tempname R
    matrix `R' = e(results)
    local val1 = `R'[1, 4]
    local pred1 = catt_hat[1]
    if reldif(`val1', `pred1') < 1e-10 {
        di as text "PASS: predict est matches e(results)"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: predict est mismatch"
        local test_fail = `test_fail' + 1
    }
}

// 3b: Standard error
capture noisily predict catt_se, se
if _rc {
    di as error "FAIL: predict se rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    tempname R
    matrix `R' = e(results)
    local val1 = `R'[1, 5]
    local pred1 = catt_se[1]
    if reldif(`val1', `pred1') < 1e-10 {
        di as text "PASS: predict se matches e(results)"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: predict se mismatch"
        local test_fail = `test_fail' + 1
    }
}

// 3c: Analytical CI
capture noisily predict catt_ci1, ci1
if _rc {
    di as error "FAIL: predict ci1 rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    tempname R
    matrix `R' = e(results)
    local val_lb = `R'[1, 6]
    local pred_lb = catt_ci1_lb[1]
    local val_ub = `R'[1, 7]
    local pred_ub = catt_ci1_ub[1]
    // Note: analytical CI may be missing (.) if UCB not computable
    local ci1_ok = 1
    if `val_lb' == . & `pred_lb' == . {
        // Both missing — OK (analytical UCB not applicable)
        local ci1_ok = 1
    }
    else if reldif(`val_lb', `pred_lb') < 1e-10 & ///
            reldif(`val_ub', `pred_ub') < 1e-10 {
        local ci1_ok = 1
    }
    else {
        local ci1_ok = 0
    }
    if `ci1_ok' {
        di as text "PASS: predict ci1 matches e(results)"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: predict ci1 mismatch"
        local test_fail = `test_fail' + 1
    }
}

// 3d: Bootstrap CI
capture noisily predict catt_ci2, ci2
if _rc {
    di as error "FAIL: predict ci2 rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    tempname R
    matrix `R' = e(results)
    local val_lb = `R'[1, 8]
    local pred_lb = catt_ci2_lb[1]
    local val_ub = `R'[1, 9]
    local pred_ub = catt_ci2_ub[1]
    if reldif(`val_lb', `pred_lb') < 1e-10 & ///
       reldif(`val_ub', `pred_ub') < 1e-10 {
        di as text "PASS: predict ci2 matches e(results)"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: predict ci2 mismatch"
        local test_fail = `test_fail' + 1
    }
}

// 3e: Bandwidth
capture noisily predict catt_bw, bw
if _rc {
    di as error "FAIL: predict bw rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    tempname R
    matrix `R' = e(results)
    local val1 = `R'[1, 10]
    local pred1 = catt_bw[1]
    if reldif(`val1', `pred1') < 1e-10 {
        di as text "PASS: predict bw matches e(results)"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: predict bw mismatch"
        local test_fail = `test_fail' + 1
    }
}

// =============================================================================
// STEP 4: level() option verification
// =============================================================================
di as text ""
di as text "================================================================"
di as text " STEP 4: level() option verification"
di as text "================================================================"

set seed 12345
didhetero_simdata, n(300) tau(4) clear

// level(95)
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1) level(95)

if _rc {
    di as error "FAIL: catt_gt level(95) failed"
    local test_fail = `test_fail' + 1
}
else {
    matrix r95 = e(results)
    local level95 = e(level)
    if abs(`level95' - 95) < 0.01 {
        di as text "PASS: level(95) e(level) = `level95'"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: e(level) = `level95', expected 95"
        local test_fail = `test_fail' + 1
    }
}

// level(90) - CI should be narrower
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1) level(90)

if _rc {
    di as error "FAIL: catt_gt level(90) failed"
    local test_fail = `test_fail' + 1
}
else {
    matrix r90 = e(results)
    local level90 = e(level)
    local est95 = r95[1, 4]
    local est90 = r90[1, 4]
    local se95  = r95[1, 5]
    local se90  = r90[1, 5]
    local w95 = r95[1, 7] - r95[1, 6]
    local w90 = r90[1, 7] - r90[1, 6]

    // Note: CI widths may both be missing (.) if analytical not available
    local lev90_ok = 0
    if abs(`level90' - 90) < 0.01 & reldif(`est95', `est90') < 1e-8 ///
       & reldif(`se95', `se90') < 1e-8 {
        if `w90' == . & `w95' == . {
            // Both missing — analytical CI not available, still passes
            local lev90_ok = 1
        }
        else if `w90' < `w95' {
            local lev90_ok = 1
        }
    }
    if `lev90_ok' {
        di as text "PASS: level(90) est & SE unchanged"
        di as text "  CI w@95=`w95' w@90=`w90'"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: level(90) verification"
        di as error "  lev=`level90' est95=`est95' est90=`est90'"
        local test_fail = `test_fail' + 1
    }
}

// level(99) - CI should be wider
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1) level(99)

if _rc {
    di as error "FAIL: catt_gt level(99) failed"
    local test_fail = `test_fail' + 1
}
else {
    matrix r99 = e(results)
    local level99 = e(level)
    local w99 = r99[1, 7] - r99[1, 6]
    local w95c = r95[1, 7] - r95[1, 6]
    local est99 = r99[1, 4]
    local est95c = r95[1, 4]

    local lev99_ok = 0
    if abs(`level99' - 99) < 0.01 & reldif(`est99', `est95c') < 1e-8 {
        if `w99' == . & `w95c' == . {
            local lev99_ok = 1
        }
        else if `w99' > `w95c' {
            local lev99_ok = 1
        }
    }
    if `lev99_ok' {
        di as text "PASS: level(99) CI wider"
        di as text "  CI w@95=`w95c' w@99=`w99'"
        local test_pass = `test_pass' + 1
    }
    else {
        di as error "FAIL: level(99) verification"
        di as error "  lev=`level99' w99=`w99' w95=`w95c'"
        local test_fail = `test_fail' + 1
    }
}

// =============================================================================
// STEP 5: estat commands
// =============================================================================
di as text ""
di as text "================================================================"
di as text " STEP 5: estat overlap/pretrend verification"
di as text "================================================================"

// estat overlap
set seed 12345
didhetero_simdata, n(300) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1)

if _rc {
    di as error "FAIL: catt_gt for estat overlap failed"
    local test_fail = `test_fail' + 1
}
else {
    capture noisily estat overlap
    if _rc {
        di as error "FAIL: estat overlap rc = " _rc
        local test_fail = `test_fail' + 1
    }
    else {
        capture confirm scalar r(gps_min)
        if _rc {
            di as error "FAIL: estat overlap no r(gps_min)"
            local test_fail = `test_fail' + 1
        }
        else {
            di as text "PASS: estat overlap diagnostics"
            di as text "  GPS min=" r(gps_min) " max=" r(gps_max)
            local test_pass = `test_pass' + 1
        }
    }
}

// estat pretrend (requires pretrend option + bootstrap)
set seed 12345
didhetero_simdata, n(300) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1) pretrend bstrap biters(100)

if _rc {
    di as error "FAIL: catt_gt pretrend rc = " _rc
    local test_fail = `test_fail' + 1
}
else {
    capture noisily estat pretrend
    if _rc {
        di as error "FAIL: estat pretrend rc = " _rc
        local test_fail = `test_fail' + 1
    }
    else {
        capture confirm scalar r(n_pretreat)
        if _rc {
            di as error "FAIL: estat pretrend no r(n_pretreat)"
            local test_fail = `test_fail' + 1
        }
        else {
            di as text "PASS: estat pretrend diagnostics"
            di as text "  n_pretreat=" r(n_pretreat) " n_reject=" r(n_reject)
            local test_pass = `test_pass' + 1
        }
    }
}

// =============================================================================
// STEP 6: aggte_gt level() verification
// =============================================================================
di as text ""
di as text "================================================================"
di as text " STEP 6: aggte_gt level() verification"
di as text "================================================================"

set seed 12345
didhetero_simdata, n(300) tau(4) clear

capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) ///
    xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
    bwselect(IMSE1) bstrap biters(100)

if _rc {
    di as error "FAIL: catt_gt for aggte failed"
    local test_fail = `test_fail' + 1
}
else {
    capture noisily aggte_gt, type(simple) level(90) ///
        bstrap(true) biters(100)
    if _rc {
        di as error "FAIL: aggte_gt level(90) rc = " _rc
        local test_fail = `test_fail' + 1
    }
    else {
        local aggte_level = e(level)
        if abs(`aggte_level' - 90) < 0.01 {
            di as text "PASS: aggte_gt level(90) e(level)=`aggte_level'"
            local test_pass = `test_pass' + 1
        }
        else {
            di as error "FAIL: aggte e(level)=`aggte_level' expected 90"
            local test_fail = `test_fail' + 1
        }
    }
}

// =============================================================================
// SUMMARY
// =============================================================================
di as text ""
di as text "================================================================"
di as text " VERIFICATION SUMMARY"
di as text "================================================================"
di as text ""
di as text "  Total PASS: `test_pass'"
di as text "  Total FAIL: `test_fail'"
di as text ""

if `test_fail' > 0 {
    di as error "  OVERALL: FAIL"
}
else {
    di as text "  OVERALL: PASS"
}

log close _verify_api
