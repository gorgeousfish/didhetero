// =============================================================================
// End-to-end test for 4 bug fixes
// =============================================================================
clear all
set more off
set seed 12345

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

// =============================================================================
// TEST 1: Basic regression test — core functionality not broken
// =============================================================================
display _n "{hline 60}"
display "TEST 1: Basic regression test (core functionality)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        nobstrap
}
if _rc == 0 {
    display as result "  [PASS] catt_gt basic estimation"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] catt_gt basic estimation (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// Test aggregation: dynamic (no bootstrap for speed)
capture noisily {
    aggte_gt, type("dynamic") bstrap("false")
}
if _rc == 0 {
    display as result "  [PASS] aggte_gt dynamic aggregation"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] aggte_gt dynamic aggregation (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// Test aggregation: simple
capture noisily {
    aggte_gt, type("simple") bstrap("false")
}
if _rc == 0 {
    display as result "  [PASS] aggte_gt simple aggregation"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] aggte_gt simple aggregation (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 2: Bug #1 — GPS convergence (gpsstrict option)
// =============================================================================
display _n "{hline 60}"
display "TEST 2: Bug #1 — GPS convergence (gpsstrict option)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) gpsstrict kernel(gau) porder(2) bwselect(IMSE1) ///
        nobstrap
}
if _rc == 0 {
    display as result "  [PASS] gpsstrict option accepted, normal data converges"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] gpsstrict with normal data (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 3: Bug #2 — Analytical critical value (narrow zeval + large bw)
// =============================================================================
display _n "{hline 60}"
display "TEST 3: Bug #2 — Analytical critical value (narrow zeval)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.49 0.50 0.51) kernel(gau) porder(2) bwselect(manual) ///
        bw(2.0) nobstrap
}
if _rc == 0 {
    display as result "  [PASS] Narrow zeval + large bw handled gracefully"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Narrow zeval + large bw (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 4: Bug #3 — KDE truncation (boundary points)
// =============================================================================
display _n "{hline 60}"
display "TEST 4: Bug #3 — KDE truncation (boundary points)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.05 0.5 0.95) kernel(gau) porder(2) ///
        bwselect(IMSE1) nobstrap
}
if _rc == 0 {
    display as result "  [PASS] Boundary zeval points handled (no crash)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Boundary zeval points (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// Test kdetrim option
capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.05 0.5 0.95) kernel(gau) porder(2) ///
        bwselect(IMSE1) kdetrim nobstrap
}
if _rc == 0 {
    display as result "  [PASS] kdetrim option accepted"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] kdetrim option (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 5: Bug #4 — Event time (non-equidistant panels)
// =============================================================================
display _n "{hline 60}"
display "TEST 5: Bug #4 — Event time (non-equidistant panels)"
display "{hline 60}"

capture noisily {
    clear
    // Create non-equidistant panel: 100 units x 4 periods
    set obs 400
    gen id = mod(_n-1, 100) + 1
    gen period_idx = floor((_n-1) / 100) + 1
    // Non-equidistant periods: 2000, 2005, 2010, 2015
    gen period = cond(period_idx==1, 2000, cond(period_idx==2, 2005, ///
        cond(period_idx==3, 2010, 2015)))
    // Treatment group: first 50 units treated at 2010, rest never treated
    gen G = cond(id <= 50, 2010, 0)
    // Z must be time-invariant: generate per-unit value
    gen Z_tmp = rnormal() if period_idx == 1
    bysort id (period): replace Z_tmp = Z_tmp[1]
    rename Z_tmp Z
    gen Y = Z + (period >= G & G > 0) * Z * 0.5 + rnormal()
    xtset id period

    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(-0.5 0 0.5) kernel(gau) porder(1) bwselect(IMSE1) ///
        nobstrap
    aggte_gt, type("dynamic") bstrap("false")

    // Check results exist
    display "Event time labels from aggte_gt dynamic:"
    matrix list e(Estimate), format(%9.4f)
}
if _rc == 0 {
    display as result "  [PASS] Non-equidistant panel event time aggregation"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Non-equidistant panel event time (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "TEST SUMMARY"
display "{hline 60}"
display "  Total:  `_test_total'"
display "  Passed: `_test_pass'"
display "  Failed: `_test_fail'"
display "{hline 60}"
if `_test_fail' > 0 {
    display as error "SOME TESTS FAILED"
    exit 1
}
else {
    display as result "ALL TESTS PASSED"
}
