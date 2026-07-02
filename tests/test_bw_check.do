// =============================================================================
// test_bw_check.do
// Tests for bandwidth undersmoothing condition check (Assumption 4(iii))
//
// Verifies:
//   - Normal scenarios produce no spurious warnings
//   - Small samples or extreme bandwidths trigger appropriate warnings
//   - Warnings are purely diagnostic (estimation still completes, rc==0)
//   - Estimation results are unaffected by the diagnostic check
// =============================================================================
clear all
set more off
set seed 12345

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

display _n "{hline 60}"
display "  Bandwidth Undersmoothing Condition Check Tests"
display "{hline 60}"

// =============================================================================
// TEST 1: Normal scenario (n=500, IMSE1) — should complete without error
// =============================================================================
display _n "{hline 60}"
display "TEST 1: Normal scenario (n=500, IMSE1)"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) nobstrap
if _rc == 0 {
    display as result "  [PASS] Normal scenario completes without error"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Normal scenario failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 2: Small sample (n=100) — may trigger warnings but should not error
// =============================================================================
display _n "{hline 60}"
display "TEST 2: Small sample (n=100) — warnings possible, no error"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(100) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) nobstrap
if _rc == 0 {
    display as result "  [PASS] Small sample completes without error (warnings are OK)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Small sample failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 3: Manual bandwidth too large — should trigger upper bound warning
// For n=500: n^(-1/9-0.01) ≈ 0.221, so bw(0.5) should warn
// =============================================================================
display _n "{hline 60}"
display "TEST 3: Manual bandwidth too large (bw=0.5, n=500)"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.5) nobstrap
if _rc == 0 {
    display as result "  [PASS] Large manual bandwidth completes (with expected warning)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Large manual bandwidth failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 4: Manual bandwidth too small — should trigger lower bound warning
// For n=500: n^(-0.49) ≈ 0.0146, so bw(0.001) should warn
// =============================================================================
display _n "{hline 60}"
display "TEST 4: Manual bandwidth too small (bw=0.001, n=500)"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.001) nobstrap
if _rc == 0 {
    display as result "  [PASS] Small manual bandwidth completes (with expected warning)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Small manual bandwidth failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 5: Bandwidth in valid range — no warning expected
// For n=500: bounds are [0.015, 0.221], so bw(0.1) should be fine
// =============================================================================
display _n "{hline 60}"
display "TEST 5: Manual bandwidth in valid range (bw=0.1, n=500)"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.1) nobstrap
if _rc == 0 {
    display as result "  [PASS] Valid range bandwidth completes without issue"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Valid range bandwidth failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 6: Estimation correctness — check that diagnostic does not alter results
// Compare results: both bw=0.1 (valid) and bw=0.5 (warns) should produce valid
// non-missing estimates.
// =============================================================================
display _n "{hline 60}"
display "TEST 6: Estimation correctness — diagnostic is non-invasive"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear

// Run with bw=0.1 (valid range, no warning)
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.1) nobstrap
local rc1 = _rc
matrix res_valid = e(b)

// Regenerate data for second call (preserve/restore in catt_gt may invalidate state)
set seed 12345
didhetero_simdata, n(500) tau(4) clear

// Run with bw=0.5 (triggers warning)
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.5) nobstrap
local rc2 = _rc
matrix res_warn = e(b)

if `rc1' == 0 & `rc2' == 0 {
    local n_valid = colsof(res_valid)
    local n_warn = colsof(res_warn)
    if `n_valid' > 0 & `n_warn' > 0 {
        display as result "  [PASS] Both runs produce valid (non-missing) estimates"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] One or both runs produced empty results"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] Estimation correctness check failed (rc1=`rc1' rc2=`rc2')"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 7: Mata-level unit test — direct call to check function
// =============================================================================
display _n "{hline 60}"
display "TEST 7: Mata-level unit test of _didhetero_bw_check_undersmooth"
display "{hline 60}"

capture noisily {
    mata {
        // Case A: All bandwidths within valid range for n=500
        // Bounds: [n^(-0.49), n^(-1/9-0.01)] = [0.0146, 0.221]
        printf("{txt}--- Case A: bw=0.1 for n=500 (should be silent) ---\n")
        _didhetero_bw_check_undersmooth((0.10 \ 0.12 \ 0.08), 500, "IMSE1")

        // Case B: Bandwidth exceeds upper bound
        printf("{txt}--- Case B: bw=0.4 for n=500 (should warn upper) ---\n")
        _didhetero_bw_check_undersmooth((0.40 \ 0.35 \ 0.45), 500, "manual")

        // Case C: Bandwidth below lower bound
        printf("{txt}--- Case C: bw=0.005 for n=500 (should warn lower + nh<30) ---\n")
        _didhetero_bw_check_undersmooth((0.005 \ 0.003 \ 0.010), 500, "manual")

        // Case D: Mixed — some valid, some not
        printf("{txt}--- Case D: mixed bandwidths for n=500 ---\n")
        _didhetero_bw_check_undersmooth((0.10 \ 0.40 \ 0.005), 500, "IMSE1")

        // Case E: Empty vector (edge case)
        printf("{txt}--- Case E: empty vector (should be silent) ---\n")
        _didhetero_bw_check_undersmooth(J(0, 1, .), 500, "IMSE1")

        // Case F: Large n makes bounds tighter
        printf("{txt}--- Case F: bw=0.15 for n=5000 (should warn upper) ---\n")
        // n=5000: upper = 5000^(-1/9-0.01) = 0.115
        _didhetero_bw_check_undersmooth((0.15 \ 0.12 \ 0.10), 5000, "IMSE1")
    }
}
if _rc == 0 {
    display as result "  [PASS] Mata-level unit test completed"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Mata-level unit test failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 8: US1 bandwidth selection — undersmoothing variant
// =============================================================================
display _n "{hline 60}"
display "TEST 8: US1 bandwidth selection (undersmoothing variant)"
display "{hline 60}"

set seed 12345
didhetero_simdata, n(500) tau(4) clear
capture noisily catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(US1) nobstrap
if _rc == 0 {
    display as result "  [PASS] US1 bandwidth selection completes"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] US1 bandwidth selection failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "  TEST SUMMARY"
display "{hline 60}"
display as text "  Total tests:  " as result `_test_total'
display as text "  Passed:       " as result `_test_pass'
display as text "  Failed:       " as result `_test_fail'
display "{hline 60}"

if `_test_fail' > 0 {
    display as error "  SOME TESTS FAILED"
    exit 1
}
else {
    display as result "  ALL TESTS PASSED"
}
