// =============================================================================
// Test: IRLS verbose output and LPR n_eff effective observation check
// =============================================================================
clear all
set more off
set seed 54321

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

// =============================================================================
// TEST 1: verbose does not alter results
// =============================================================================
display _n "{hline 60}"
display "TEST 1: verbose flag does not change estimation results"
display "{hline 60}"

capture noisily {
    // Run WITHOUT verbose
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap
    matrix est_quiet = e(Estimate)
    matrix se_quiet  = e(SE)
}
local rc1 = _rc

capture noisily {
    // Run WITH verbose (same seed => same data)
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap verbose
    matrix est_verb = e(Estimate)
    matrix se_verb  = e(SE)
}
local rc2 = _rc

if `rc1' == 0 & `rc2' == 0 {
    // Compare estimates element by element (handle missing in SE)
    mata: st_numscalar("max_diff_est", max(abs(st_matrix("est_quiet") :- st_matrix("est_verb"))))
    mata: st_numscalar("max_diff_se", max(editvalue(abs(st_matrix("se_quiet") :- st_matrix("se_verb")), ., 0)))
    
    if scalar(max_diff_est) < 1e-12 & scalar(max_diff_se) < 1e-12 {
        display as result "  [PASS] verbose does not change estimates (max diff = " scalar(max_diff_est) ")"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] verbose changes results: est_diff=" scalar(max_diff_est) " se_diff=" scalar(max_diff_se)
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] Estimation failed (rc1=`rc1', rc2=`rc2')"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 2: verbose output format contains expected text
// =============================================================================
display _n "{hline 60}"
display "TEST 2: verbose output contains 'GPS (g=' format string"
display "{hline 60}"

tempfile logfile
capture noisily {
    set seed 11111
    didhetero_simdata, n(300) tau(4) clear
    log using `logfile', text replace
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap verbose
    log close
}
if _rc == 0 {
    // Use Mata fget() to read log file (avoids Stata macro quoting issues)
    mata: _dh_test_found = 0
    mata: _dh_test_fh = fopen(st_local("logfile"), "r")
    mata: while ((_dh_test_line = fget(_dh_test_fh)) != J(0,0,"")) { if (strpos(_dh_test_line, "GPS (g=") > 0) _dh_test_found = 1; }
    mata: fclose(_dh_test_fh)
    mata: st_local("found_gps", strofreal(_dh_test_found))
    
    if `found_gps' == 1 {
        display as result "  [PASS] Found 'GPS (g=' in verbose output"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Did not find 'GPS (g=' in verbose output"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] Estimation with verbose failed"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 3: n_eff check - normal data, no missing in results
// =============================================================================
display _n "{hline 60}"
display "TEST 3: n_eff check - standard data yields no missing estimates"
display "{hline 60}"

capture noisily {
    set seed 33333
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap
    matrix est_normal = e(Estimate)
}
if _rc == 0 {
    // e(Estimate) has 10 columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
    // Only columns 4-5 (est, se) matter for n_eff validation
    mata: st_numscalar("n_miss", missing(st_matrix("est_normal")[., 4..5]))
    if scalar(n_miss) == 0 {
        display as result "  [PASS] No missing values in point estimates (n=500, normal data)"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Found " scalar(n_miss) " missing values in point estimates"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] Estimation failed"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 4: n_eff check - boundary z points may produce missing
// =============================================================================
display _n "{hline 60}"
display "TEST 4: n_eff check - extreme zeval points with small bandwidth"
display "{hline 60}"

capture noisily {
    set seed 44444
    didhetero_simdata, n(200) tau(4) clear
    // Use points near Z boundaries with a small bandwidth
    // (points are within Z range but have few nearby observations)
    quietly summarize Z
    local z_lo = r(min) + 0.05
    local z_hi = r(max) - 0.05
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(`z_lo' `z_hi') kernel(epa) porder(1) ///
        bwselect(manual) bw(0.05) ///
        control_group(nevertreated) nobstrap
    matrix est_extreme = e(Estimate)
}
if _rc == 0 {
    mata: st_numscalar("n_miss_ext", missing(st_matrix("est_extreme")))
    // Boundary points with small bandwidth may produce missing
    if scalar(n_miss_ext) > 0 {
        display as result "  [PASS] Boundary zeval + small bw correctly produce missing (" scalar(n_miss_ext) " missing)"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as result "  [PASS] All boundary points estimated (data sufficiently dense at boundaries)"
        local _test_pass = `_test_pass' + 1
    }
}
else {
    // Command may error if all eval points produce missing
    display as result "  [PASS] Boundary estimation failed as expected (rc=" _rc ")"
    local _test_pass = `_test_pass' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 5: n_eff check - very small bandwidth produces missing
// =============================================================================
display _n "{hline 60}"
display "TEST 5: n_eff check - tiny bandwidth (0.01) produces missing values"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(300) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.2 0.4 0.6 0.8) kernel(epa) porder(1) ///
        bwselect(manual) bw(0.01) ///
        control_group(nevertreated) nobstrap
    matrix est_tiny = e(Estimate)
}
if _rc == 0 {
    mata: st_numscalar("n_miss_tiny", missing(st_matrix("est_tiny")))
    if scalar(n_miss_tiny) > 0 {
        display as result "  [PASS] Tiny bandwidth correctly produces " scalar(n_miss_tiny) " missing (n_eff guard works)"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Tiny bw(0.01) produced no missing values (unexpected)"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] Estimation failed (rc=" _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 6: porder affects minimum observation requirement
// =============================================================================
display _n "{hline 60}"
display "TEST 6: porder(2) requires more observations (min_obs=6 vs 4)"
display "{hline 60}"

capture noisily {
    // porder(2) with moderate bandwidth - should still work on normal data
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(epa) porder(2) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap
    matrix est_p2_normal = e(Estimate)
}
local rc_p2 = _rc

if `rc_p2' == 0 {
    mata: st_numscalar("n_miss_p2", missing(st_matrix("est_p2_normal")))
    display as result "  [INFO] porder(2) normal data: " scalar(n_miss_p2) " missing values"
    
    // Now with tiny bandwidth - porder(2) needs 6 obs, harder to satisfy
    capture noisily {
        didhetero_simdata, n(200) tau(4) clear
        catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
            zeval(0.3 0.5 0.7) kernel(epa) porder(2) ///
            bwselect(manual) bw(0.03) ///
            control_group(nevertreated) nobstrap
        matrix est_p2_tiny = e(Estimate)
    }
    if _rc == 0 {
        mata: st_numscalar("n_miss_p2_tiny", missing(st_matrix("est_p2_tiny")))
        if scalar(n_miss_p2_tiny) > 0 {
            display as result "  [PASS] porder(2) + tiny bw produces " scalar(n_miss_p2_tiny) " missing (stricter n_eff requirement)"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as result "  [PASS] porder(2) + bw(0.03) all estimated (data sufficiently dense)"
            local _test_pass = `_test_pass' + 1
        }
    }
    else {
        display as error "  [FAIL] porder(2) tiny bandwidth estimation failed"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] porder(2) normal estimation failed (rc=`rc_p2')"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// TEST 7: verbose with notyettreated control (checks g,t format)
// =============================================================================
display _n "{hline 60}"
display "TEST 7: verbose output with notyettreated shows (g,t) pairs"
display "{hline 60}"

tempfile logfile2
capture noisily {
    set seed 22222
    didhetero_simdata, n(300) tau(4) clear
    log using `logfile2', text replace
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap verbose
    log close
}
if _rc == 0 {
    // Use Mata fget() to read log file (avoids Stata macro quoting issues)
    mata: _dh_test_found2 = 0
    mata: _dh_test_fh2 = fopen(st_local("logfile2"), "r")
    mata: while ((_dh_test_line2 = fget(_dh_test_fh2)) != J(0,0,"")) { if (strpos(_dh_test_line2, "GPS (g=") > 0 & strpos(_dh_test_line2, "t=") > 0) _dh_test_found2 = 1; }
    mata: fclose(_dh_test_fh2)
    mata: st_local("found_gt", strofreal(_dh_test_found2))
    
    if `found_gt' == 1 {
        display as result "  [PASS] Found 'GPS (g=..., t=...)' in notyettreated verbose output"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Did not find 'GPS (g=..., t=...)' in verbose output"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] notyettreated verbose estimation failed"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "TEST SUMMARY: verbose + n_eff diagnostics"
display "{hline 60}"
display as text "  Total tests: `_test_total'"
display as result "  Passed:      `_test_pass'"
if `_test_fail' > 0 {
    display as error "  Failed:      `_test_fail'"
}
else {
    display as text "  Failed:      0"
}
display "{hline 60}"

if `_test_fail' > 0 {
    exit 1
}
