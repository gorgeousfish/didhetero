// =============================================================================
// test_rbc.do — Test suite for Simple RBC (Robust Bias Correction) option
//
// Tests:
//   1. rbc option runs without error
//   2. rbc reports correct e() scalars (porder=2, rbc=1, bwselect="IMSE1")
//   3. rbc conflict detection (bw(), bwselect(manual))
//   4. rbc compatibility with other options
//   5. didhetero command rbc support
//   6. aggte_gt inheritance from rbc estimation
// =============================================================================

clear all
set more off
set seed 12345

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

// Generate simulated data
didhetero_simdata, n(300) tau(4) clear

// =============================================================================
// TEST 1: rbc option runs without error
// =============================================================================
display _n "{hline 60}"
display "TEST 1: rbc option basic execution"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc nobstrap
}

if _rc == 0 {
    display as result "  [PASS] rbc option executed without error"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] rbc option returned rc = " _rc
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 2: rbc reports correct e() scalars
// =============================================================================
display _n "{hline 60}"
display "TEST 2: rbc e() scalar reporting"
display "{hline 60}"
local _test_total = `_test_total' + 1

local _t2_pass = 1

// e(porder) should be 2 (estimation polynomial order)
if e(porder) != 2 {
    display as error "  [FAIL] e(porder) = " e(porder) " (expected 2)"
    local _t2_pass = 0
}

// e(rbc) should be 1
if e(rbc) != 1 {
    display as error "  [FAIL] e(rbc) = " e(rbc) " (expected 1)"
    local _t2_pass = 0
}

// e(bwselect) should be "IMSE1"
if "`e(bwselect)'" != "IMSE1" {
    display as error "  [FAIL] e(bwselect) = `e(bwselect)' (expected IMSE1)"
    local _t2_pass = 0
}

if `_t2_pass' {
    display as result "  [PASS] e(porder)=2, e(rbc)=1, e(bwselect)=IMSE1"
    local _test_pass = `_test_pass' + 1
}
else {
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 3: rbc conflict detection — bw() + bwselect(manual)
// =============================================================================
display _n "{hline 60}"
display "TEST 3: rbc + bw() + bwselect(manual) conflict"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc bw(0.3) bwselect(manual) nobstrap
}

if _rc != 0 {
    display as result "  [PASS] rbc + bw() + bwselect(manual) correctly rejected (rc = " _rc ")"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] rbc + bw() + bwselect(manual) should have been rejected"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 4: rbc conflict detection — bwselect(manual) without bw
// =============================================================================
display _n "{hline 60}"
display "TEST 4: rbc + bwselect(manual) conflict"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc bwselect(manual) nobstrap
}

if _rc != 0 {
    display as result "  [PASS] rbc + bwselect(manual) correctly rejected (rc = " _rc ")"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] rbc + bwselect(manual) should have been rejected"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 5: rbc compatible with gpstrim
// =============================================================================
display _n "{hline 60}"
display "TEST 5: rbc + gpstrim compatibility"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc gpstrim(0.01 0.99) nobstrap
}

if _rc == 0 {
    display as result "  [PASS] rbc + gpstrim works"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] rbc + gpstrim failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 6: rbc + control_group(nevertreated) compatibility
// =============================================================================
display _n "{hline 60}"
display "TEST 6: rbc + control_group(nevertreated)"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc control_group(nevertreated) nobstrap
}

if _rc == 0 {
    display as result "  [PASS] rbc + control_group(nevertreated) works"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] rbc + control_group(nevertreated) failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 7: without rbc, e(rbc) = 0
// =============================================================================
display _n "{hline 60}"
display "TEST 7: without rbc, e(rbc) = 0"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) nobstrap
}

if _rc == 0 & e(rbc) == 0 {
    display as result "  [PASS] e(rbc) = 0 without rbc option"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] e(rbc) = " e(rbc) " (expected 0 without rbc)"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 8: didhetero command rbc support
// =============================================================================
display _n "{hline 60}"
display "TEST 8: didhetero command with rbc"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc bstrap(false)
}

if _rc == 0 {
    local _t8_pass = 1
    if e(rbc) != 1 {
        display as error "  [FAIL] didhetero e(rbc) = " e(rbc) " (expected 1)"
        local _t8_pass = 0
    }
    if e(porder) != 2 {
        display as error "  [FAIL] didhetero e(porder) = " e(porder) " (expected 2)"
        local _t8_pass = 0
    }
    if `_t8_pass' {
        display as result "  [PASS] didhetero with rbc works correctly"
        local _test_pass = `_test_pass' + 1
    }
    else {
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] didhetero with rbc failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 9: aggte_gt after rbc estimation
// =============================================================================
display _n "{hline 60}"
display "TEST 9: aggte_gt after rbc catt_gt"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc nobstrap
}

if _rc == 0 {
    capture noisily aggte_gt, type(simple)
    if _rc == 0 {
        display as result "  [PASS] aggte_gt(simple) works after rbc estimation"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] aggte_gt failed after rbc estimation (rc = " _rc ")"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] prerequisite catt_gt with rbc failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 10: didhetero with rbc + pretrend
// =============================================================================
display _n "{hline 60}"
display "TEST 10: didhetero with rbc + pretrend"
display "{hline 60}"
local _test_total = `_test_total' + 1

capture noisily {
    didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) rbc pretrend bstrap(false)
}

if _rc == 0 {
    display as result "  [PASS] didhetero with rbc + pretrend works"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] didhetero with rbc + pretrend failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// Summary
// =============================================================================
display _n "{hline 60}"
display "RBC Test Summary"
display "{hline 60}"
display "  Total tests:  `_test_total'"
display as result "  Passed:       `_test_pass'"
if `_test_fail' > 0 {
    display as error "  Failed:       `_test_fail'"
}
else {
    display "  Failed:       `_test_fail'"
}
display "{hline 60}"

if `_test_fail' > 0 {
    display as error _n "SOME TESTS FAILED"
    exit 9
}
else {
    display as result _n "ALL TESTS PASSED"
}
