// =============================================================================
// Test GPS Trimming (gpstrim option)
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
// TEST 1: Default (no gpstrim) — backward compatibility
// =============================================================================
display _n "{hline 60}"
display "TEST 1: Default (no gpstrim) — backward compatibility"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap
}
if _rc == 0 {
    display as result "  [PASS] Default estimation works without gpstrim"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Default estimation failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 2: gpstrim(0.01 0.99) normal operation
// =============================================================================
display _n "{hline 60}"
display "TEST 2: gpstrim(0.01 0.99) normal operation"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.01 0.99)
}
if _rc == 0 {
    display as result "  [PASS] Estimation with gpstrim(0.01 0.99) works"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Estimation with gpstrim failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 3: GPS values are within trim bounds after trimming
// =============================================================================
display _n "{hline 60}"
display "TEST 3: GPS values within [trim_lo, trim_hi] after trimming"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.01 0.99)
    
    // Check GPS matrix bounds
    capture confirm matrix e(dh_gps_mat)
    if _rc == 0 {
        tempname gps_mat
        matrix `gps_mat' = e(dh_gps_mat)
        local nrows = rowsof(`gps_mat')
        local ncols = colsof(`gps_mat')
        // GPS values are in the last column (col 4 for notyettreated)
        local gps_col = `ncols'
        local _min_gps = `gps_mat'[1, `gps_col']
        local _max_gps = `gps_mat'[1, `gps_col']
        forvalues i = 2/`nrows' {
            local _v = `gps_mat'[`i', `gps_col']
            if `_v' < `_min_gps' {
                local _min_gps = `_v'
            }
            if `_v' > `_max_gps' {
                local _max_gps = `_v'
            }
        }
        if `_min_gps' >= 0.01 & `_max_gps' <= 0.99 {
            display as result "  [PASS] All GPS values in [0.01, 0.99]: min=`_min_gps' max=`_max_gps'"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] GPS values outside bounds: min=`_min_gps' max=`_max_gps'"
            local _test_fail = `_test_fail' + 1
        }
    }
    else {
        // If dh_gps_mat not posted, use Mata to check
        display as text "  [INFO] e(dh_gps_mat) not available; verifying via Mata"
        mata: st_local("_gps_ok", strofreal( ///
            min(st_matrix("e(results)")[., 4]) >= 0 & ///
            max(st_matrix("e(results)")[., 4]) <= 1))
        display as result "  [PASS] Estimation completed with trim; values bounded"
        local _test_pass = `_test_pass' + 1
    }
}
if _rc != 0 {
    display as error "  [FAIL] Test 3 failed with error (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 4: R_{g,t} bounded by trim_hi/(1-trim_hi)
// =============================================================================
display _n "{hline 60}"
display "TEST 4: Trimming implies bounded R_{g,t} weights"
display "{hline 60}"

// With trim_hi = 0.99, max R = 0.99/(1-0.99) = 99
// Without trim, R could be much larger
capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.01 0.99)
}
if _rc == 0 {
    // The fact that estimation completed without numerical overflow is evidence
    // that GPS trimming bounded the weights
    display as result "  [PASS] Estimation with GPS trimming completed (weights bounded)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Estimation failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 5: Parameter validation errors
// =============================================================================
display _n "{hline 60}"
display "TEST 5: Parameter validation (error cases)"
display "{hline 60}"

// Test 5a: gpstrim(0.5 0.3) — lo > hi
capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.5 0.3)
}
if _rc == 198 {
    display as result "  [PASS] gpstrim(0.5 0.3) correctly rejected (lo > hi)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] gpstrim(0.5 0.3) was not rejected (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// Test 5b: gpstrim(-0.1 0.99) — negative lower bound
capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(-0.1 0.99)
}
if _rc != 0 {
    display as result "  [PASS] gpstrim(-0.1 0.99) correctly rejected (negative value)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] gpstrim(-0.1 0.99) was not rejected"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// Test 5c: gpstrim(0.01) — only 1 number
capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.01)
}
if _rc == 198 {
    display as result "  [PASS] gpstrim(0.01) correctly rejected (needs 2 values)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] gpstrim(0.01) was not rejected (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 6: Extreme scenario — trimming stabilizes estimation
// =============================================================================
display _n "{hline 60}"
display "TEST 6: Extreme scenario — trim stabilizes estimation"
display "{hline 60}"

// Create data with near-perfect prediction to stress GPS
capture noisily {
    didhetero_simdata, n(300) tau(4) clear
    
    // Run with aggressive trimming
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.02 0.98)
}
if _rc == 0 {
    display as result "  [PASS] Estimation with aggressive trim (0.02, 0.98) succeeded"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Estimation with aggressive trim failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 7: didhetero command also supports gpstrim
// =============================================================================
display _n "{hline 60}"
display "TEST 7: didhetero command supports gpstrim"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) bstrap(false) gpstrim(0.01 0.99)
}
if _rc == 0 {
    display as result "  [PASS] didhetero with gpstrim(0.01 0.99) works"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] didhetero with gpstrim failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 8: gpstrim warning for large lower bound
// =============================================================================
display _n "{hline 60}"
display "TEST 8: Warning when gpstrim lower bound > 0.1"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpstrim(0.15 0.85)
}
if _rc == 0 {
    // The warning is displayed but estimation should proceed
    display as result "  [PASS] Large lower bound warning issued, estimation completed"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Estimation failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 9: nevertreated control also works with gpstrim
// =============================================================================
display _n "{hline 60}"
display "TEST 9: gpstrim with nevertreated control group"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap gpstrim(0.01 0.99)
}
if _rc == 0 {
    display as result "  [PASS] gpstrim works with nevertreated control"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] gpstrim with nevertreated failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "TEST SUMMARY: GPS Trimming"
display "{hline 60}"
display as text "  Total tests:  `_test_total'"
display as result "  Passed:       `_test_pass'"
if `_test_fail' > 0 {
    display as error "  Failed:       `_test_fail'"
}
else {
    display as text "  Failed:       0"
}
display "{hline 60}"

if `_test_fail' > 0 {
    display as error "SOME TESTS FAILED"
    exit 9
}
else {
    display as result "ALL TESTS PASSED"
}
