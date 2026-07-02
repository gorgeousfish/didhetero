// =============================================================================
// test_profile.do
// Tests for the performance profiling framework (profile option)
//
// Tests:
//   1. profile does not change estimation results
//   2. profile output contains expected format string "Performance Profile"
//   3. No bootstrap → bootstrap time = 0.00s
//   4. didhetero command also supports profile
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

didhetero_simdata, n(300) tau(4) clear

display _n "{hline 70}"
display "TEST SUITE: Performance Profiling (profile option)"
display "{hline 70}"

// =============================================================================
// TEST 1: profile does not change estimation results
// =============================================================================
display _n "TEST 1: profile does not change estimation results"

// Run WITHOUT profile
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-1 0 1) control_group(notyettreated) ///
    nobstrap nouniformall bwselect(IMSE1) kernel(gau) porder(2)

matrix _noprof_results = e(results)
scalar _noprof_N = e(N)

// Run WITH profile
quietly catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-1 0 1) control_group(notyettreated) ///
    nobstrap nouniformall bwselect(IMSE1) kernel(gau) porder(2) profile

matrix _prof_results = e(results)
scalar _prof_N = e(N)

// Compare results
local test1_pass 1
if _noprof_N != _prof_N {
    local test1_pass 0
}

// Compare dimensions
local noprof_rows = rowsof(_noprof_results)
local prof_rows = rowsof(_prof_results)
if `noprof_rows' != `prof_rows' {
    local test1_pass 0
}

// Compare values (point estimates, column 4)
if `test1_pass' {
    forvalues i = 1/`noprof_rows' {
        if reldif(_noprof_results[`i', 4], _prof_results[`i', 4]) > 1e-10 {
            local test1_pass 0
        }
    }
}

if `test1_pass' {
    display as text "  PASS: profile option does not alter estimation results"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: profile option changed estimation results"
    local n_fail = `n_fail' + 1
}

// =============================================================================
// TEST 2: profile output contains "Performance Profile" and stage labels
// =============================================================================
display _n "TEST 2: profile output contains expected format"

// Write a small Stata script that runs profile and checks output via grep
tempfile logfile
tempfile grepresult

log using "`logfile'", text replace

catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
    zeval(-1 0 1) control_group(notyettreated) ///
    nobstrap nouniformall bwselect(IMSE1) kernel(gau) porder(2) profile

log close

// Use shell grep to check patterns (robust to Stata quoting issues)
shell grep -c "Performance Profile" "`logfile'" > "`grepresult'"
scalar _grep_profile = 0
capture {
    infix str v1 1-20 using "`grepresult'", clear
    scalar _grep_profile = real(v1[1])
}

didhetero_simdata, n(300) tau(4) clear

shell grep -c "GPS/OR estimation" "`logfile'" > "`grepresult'"
scalar _grep_gps = 0
capture {
    infix str v1 1-20 using "`grepresult'", clear
    scalar _grep_gps = real(v1[1])
}

didhetero_simdata, n(300) tau(4) clear

shell grep -c "CATT estimation" "`logfile'" > "`grepresult'"
scalar _grep_catt = 0
capture {
    infix str v1 1-20 using "`grepresult'", clear
    scalar _grep_catt = real(v1[1])
}

didhetero_simdata, n(300) tau(4) clear

shell grep -c "Total:" "`logfile'" > "`grepresult'"
scalar _grep_total = 0
capture {
    infix str v1 1-20 using "`grepresult'", clear
    scalar _grep_total = real(v1[1])
}

didhetero_simdata, n(300) tau(4) clear

if _grep_profile > 0 & _grep_gps > 0 & _grep_catt > 0 & _grep_total > 0 {
    display as text "  PASS: profile output contains all expected sections"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: missing profile output sections"
    display as error "    profile=" _grep_profile " gps=" _grep_gps " catt=" _grep_catt " total=" _grep_total
    local n_fail = `n_fail' + 1
}

// =============================================================================
// TEST 3: No bootstrap → bootstrap line shows 0.00s
// =============================================================================
display _n "TEST 3: No bootstrap → bootstrap time = 0.00s"

shell grep "Bootstrap:" "`logfile'" | grep -c "0.00s" > "`grepresult'"
scalar _grep_bszero = 0
capture {
    infix str v1 1-20 using "`grepresult'", clear
    scalar _grep_bszero = real(v1[1])
}

didhetero_simdata, n(300) tau(4) clear

if _grep_bszero > 0 {
    display as text "  PASS: bootstrap time is 0.00s when nobstrap specified"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: bootstrap time should be 0.00s with nobstrap"
    local n_fail = `n_fail' + 1
}

// =============================================================================
// TEST 4: didhetero command also supports profile
// =============================================================================
display _n "TEST 4: didhetero command supports profile option"

local test4_pass 1
capture noisily {
    quietly didhetero Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(-1 0 1) control_group(notyettreated) ///
        nobootstrap nouniformall bwselect(IMSE1) kernel(gau) porder(2) profile
}
if _rc {
    local test4_pass 0
}

if `test4_pass' {
    display as text "  PASS: didhetero command accepts profile option"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: didhetero command rejected profile option"
    local n_fail = `n_fail' + 1
}

// =============================================================================
// SUMMARY
// =============================================================================

display _n "{hline 70}"
display "TEST SUMMARY: `n_pass' passed, `n_fail' failed"
display "{hline 70}"

if `n_fail' > 0 {
    exit 1
}
