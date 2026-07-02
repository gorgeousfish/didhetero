// =============================================================================
// test_undersmooth.do
// Tests for bandwidth undersmoothing adaptive adjustment (Assumption 4(iii))
//
// Tests:
//   1. Normal data: undersmooth option does not change result when BW in range
//   2. Small sample: undersmooth triggers upper bound adjustment
//   3. Manual large BW + undersmooth: clamped to upper bound
//   4. Manual small BW + undersmooth: clamped to lower bound
//   5. Without undersmooth: only warnings, BW unchanged (backward compat)
//   6. undersmooth + rbc combination works
//   7. Adjustment Note output is informational (non-error)
//   8. e(bw_adjusted) scalar correctness
// =============================================================================

clear all
set more off

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

// Ensure compiled library is available
capture noisily do "/Users/cxy/Desktop/didhetero/didhetero-main/make.do"
if _rc {
    display as error "make.do failed — cannot proceed with tests"
    exit _rc
}

local n_pass = 0
local n_fail = 0

// =============================================================================
// Generate test data using didhetero_simdata
// =============================================================================

// --- Standard data (n=500) ---
didhetero_simdata, n(500) tau(4) clear

display _n "{hline 70}"
display "TEST SUITE: Bandwidth Undersmoothing Adaptive Adjustment"
display "{hline 70}"

// =============================================================================
// TEST 1: Normal data — undersmooth should not change results if BW in range
// =============================================================================
display _n "TEST 1: Normal n=500 data — undersmooth should not alter in-range BW"

// Run without undersmooth
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-1 0 1) control_group(notyettreated) ///
    bstrap(false) uniformall(false)
matrix bw_no_us = e(bw)
scalar est_no_us = el(e(results), 1, 4)

// Run with undersmooth
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-1 0 1) control_group(notyettreated) ///
    bstrap(false) uniformall(false) undersmooth
matrix bw_with_us = e(bw)
scalar est_with_us = el(e(results), 1, 4)
scalar bw_adj_1 = e(bw_adjusted)

// For n=500, IMSE bandwidth should be within [n^{-0.49}, n^{-0.121}]
// n^{-0.49} ≈ 0.046, n^{-0.121} ≈ 0.517
// Typical IMSE bandwidth is ~0.1-0.4, so should be in range
scalar h_lower_500 = 500^(-0.49)
scalar h_upper_500 = 500^(-1.0/9.0 - 0.01)

// Check BW values are the same (no adjustment needed)
local bw1_no = el(bw_no_us, 1, 1)
local bw1_with = el(bw_with_us, 1, 1)

if abs(`bw1_no' - `bw1_with') < 1e-10 {
    display as result "  PASS: BW unchanged when already in range"
    local ++n_pass
}
else {
    // BW may differ if it was slightly outside range
    if `bw1_with' >= scalar(h_lower_500) & `bw1_with' <= scalar(h_upper_500) {
        display as result "  PASS: BW adjusted correctly to within bounds"
        local ++n_pass
    }
    else {
        display as error "  FAIL: unexpected BW change"
        display as error "    BW without undersmooth: `bw1_no'"
        display as error "    BW with undersmooth:    `bw1_with'"
        local ++n_fail
    }
}

// =============================================================================
// TEST 2: Small sample — check upper bound adjustment triggers
// =============================================================================
display _n "TEST 2: Small sample n=100 — upper bound adjustment"

clear
didhetero_simdata, n(100) tau(4) clear

// For n=100: upper bound = 100^{-0.121} ≈ 0.572
scalar h_upper_100 = 100^(-1.0/9.0 - 0.01)
scalar h_lower_100 = 100^(-0.49)

quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bstrap(false) uniformall(false) undersmooth
matrix bw_small = e(bw)
scalar bw_adj_2 = e(bw_adjusted)

// All bandwidths should be <= upper bound
local test2_pass = 1
forvalues k = 1/`=colsof(bw_small)' {
    local bw_k = el(bw_small, 1, `k')
    if `bw_k' > scalar(h_upper_100) + 1e-10 {
        local test2_pass = 0
        display as error "  BW[`k'] = `bw_k' > upper bound = " scalar(h_upper_100)
    }
}

if `test2_pass' {
    display as result "  PASS: All BW <= n^{-1/9-eps} upper bound"
    local ++n_pass
}
else {
    display as error "  FAIL: Some BW exceed upper bound"
    local ++n_fail
}

// =============================================================================
// TEST 3: Manual large BW + undersmooth → clamped to upper bound
// =============================================================================
display _n "TEST 3: Manual bw(0.8) + undersmooth → clamped"

// n=100, upper bound ≈ 0.572
// bw = 0.8 should be clamped to upper bound
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.8) ///
    bstrap(false) uniformall(false) undersmooth
matrix bw_manual_large = e(bw)
scalar bw_adj_3 = e(bw_adjusted)

local bw3 = el(bw_manual_large, 1, 1)
if abs(`bw3' - scalar(h_upper_100)) < 1e-10 {
    display as result "  PASS: BW clamped from 0.8 to upper bound " scalar(h_upper_100)
    local ++n_pass
}
else if `bw3' <= scalar(h_upper_100) + 1e-10 {
    display as result "  PASS: BW adjusted within bounds (= `bw3')"
    local ++n_pass
}
else {
    display as error "  FAIL: BW not clamped (= `bw3', expected <= " scalar(h_upper_100) ")"
    local ++n_fail
}

// e(bw_adjusted) should be 1
if scalar(bw_adj_3) == 1 {
    display as result "  PASS: e(bw_adjusted) = 1 (adjustment occurred)"
    local ++n_pass
}
else {
    display as error "  FAIL: e(bw_adjusted) = " scalar(bw_adj_3) " (expected 1)"
    local ++n_fail
}

// =============================================================================
// TEST 4: Manual small BW + undersmooth → clamped to lower bound
// =============================================================================
display _n "TEST 4: Manual bw(0.001) + undersmooth → clamped to lower bound"

// n=100, lower bound = 100^{-0.49} ≈ 0.105
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.001) ///
    bstrap(false) uniformall(false) undersmooth
matrix bw_manual_small = e(bw)
scalar bw_adj_4 = e(bw_adjusted)

local bw4 = el(bw_manual_small, 1, 1)
if `bw4' >= scalar(h_lower_100) - 1e-10 {
    display as result "  PASS: BW raised from 0.001 to lower bound (= `bw4')"
    local ++n_pass
}
else {
    display as error "  FAIL: BW not raised to lower bound (= `bw4', expected >= " scalar(h_lower_100) ")"
    local ++n_fail
}

// Also verify it's > 0.001 (was adjusted)
if `bw4' > 0.001 {
    display as result "  PASS: BW > original 0.001"
    local ++n_pass
}
else {
    display as error "  FAIL: BW not increased from 0.001"
    local ++n_fail
}

// =============================================================================
// TEST 5: Without undersmooth — backward compatibility (only warnings)
// =============================================================================
display _n "TEST 5: Without undersmooth — BW unchanged (backward compat)"

// Use same manual large BW without undersmooth
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.8) ///
    bstrap(false) uniformall(false)
matrix bw_no_adj = e(bw)

local bw5 = el(bw_no_adj, 1, 1)
if abs(`bw5' - 0.8) < 1e-10 {
    display as result "  PASS: BW remains 0.8 without undersmooth option"
    local ++n_pass
}
else {
    display as error "  FAIL: BW changed without undersmooth (= `bw5')"
    local ++n_fail
}

// e(bw_adjusted) should not exist (or be 0) when undersmooth is off
capture confirm scalar e(bw_adjusted)
if _rc != 0 {
    display as result "  PASS: e(bw_adjusted) not set (flag off)"
    local ++n_pass
}
else {
    if e(bw_adjusted) == 0 {
        display as result "  PASS: e(bw_adjusted) = 0 (no adjustment)"
        local ++n_pass
    }
    else {
        display as error "  FAIL: e(bw_adjusted) = " e(bw_adjusted) " (expected 0 or not set)"
        local ++n_fail
    }
}

// =============================================================================
// TEST 6: undersmooth + rbc combination
// =============================================================================
display _n "TEST 6: undersmooth + rbc combination"

capture noisily {
    quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(-0.5 0 0.5) control_group(notyettreated) ///
        bstrap(false) uniformall(false) rbc undersmooth
}
if _rc == 0 {
    display as result "  PASS: undersmooth + rbc combination runs without error"
    local ++n_pass
    // Verify rbc flag is set
    if e(rbc) == 1 {
        display as result "  PASS: e(rbc) = 1"
        local ++n_pass
    }
    else {
        display as error "  FAIL: e(rbc) != 1"
        local ++n_fail
    }
    // Verify undersmooth flag is set
    if e(undersmooth) == 1 {
        display as result "  PASS: e(undersmooth) = 1"
        local ++n_pass
    }
    else {
        display as error "  FAIL: e(undersmooth) != 1"
        local ++n_fail
    }
}
else {
    display as error "  FAIL: undersmooth + rbc combination failed (rc = " _rc ")"
    local ++n_fail
}

// =============================================================================
// TEST 7: Adjustment Note output is informational
// =============================================================================
display _n "TEST 7: Adjustment produces Note (not error)"

// Manual large BW with undersmooth — should produce Note output
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.8) ///
    bstrap(false) uniformall(false) undersmooth
if _rc == 0 {
    display as result "  PASS: Command completes successfully (Note is informational)"
    local ++n_pass
}
else {
    display as error "  FAIL: Command failed with rc = " _rc
    local ++n_fail
}

// =============================================================================
// TEST 8: e(bw_adjusted) scalar correctness
// =============================================================================
display _n "TEST 8: e(bw_adjusted) scalar"

// Case A: adjustment occurred
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.8) ///
    bstrap(false) uniformall(false) undersmooth
if e(bw_adjusted) == 1 {
    display as result "  PASS: e(bw_adjusted) = 1 when adjustment occurs"
    local ++n_pass
}
else {
    display as error "  FAIL: e(bw_adjusted) = " e(bw_adjusted) " (expected 1)"
    local ++n_fail
}

// Case B: no adjustment (BW in range)
clear
didhetero_simdata, n(500) tau(4) clear
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-0.5 0 0.5) control_group(notyettreated) ///
    bwselect(manual) bw(0.2) ///
    bstrap(false) uniformall(false) undersmooth

// For n=500: lower=500^{-0.49}≈0.046, upper=500^{-0.121}≈0.517
// bw=0.2 should be in range
if e(bw_adjusted) == 0 {
    display as result "  PASS: e(bw_adjusted) = 0 when no adjustment needed"
    local ++n_pass
}
else {
    display as error "  FAIL: e(bw_adjusted) = " e(bw_adjusted) " (expected 0)"
    local ++n_fail
}

// =============================================================================
// Summary
// =============================================================================
display _n "{hline 70}"
display "TEST RESULTS: `n_pass' passed, `n_fail' failed"
display "{hline 70}"

if `n_fail' > 0 {
    display as error "SOME TESTS FAILED"
    exit 9
}
else {
    display as result "ALL TESTS PASSED"
}
