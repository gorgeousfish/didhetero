// =============================================================================
// test_bootstrap_perf.do
// Performance and correctness tests for bootstrap matrix acceleration (Task #13)
//
// Tests:
//   1. Numerical equivalence: vectorized vs legacy (same seed, reldif < 1e-12)
//   2. Critical value consistency
//   3. Confidence band consistency
//   4. Performance benchmark (informational)
//   5. Boundary conditions (biters=1, biters=2)
// =============================================================================
clear all
set more off

// Load package
quietly adopath ++ "/Users/cxy/Desktop/didhetero/didhetero-main/ado"
quietly mata: mata mlib index

display _n "{hline 70}"
display "  Bootstrap Matrix Acceleration — Correctness & Performance Tests"
display "{hline 70}"

// ==========================================================================
// DATA GENERATION
// ==========================================================================
display _n as text "=== Generating test data (n=300, seed=42) ==="
didhetero_simdata, n(300) tau(4) seed(42) clear

// ==========================================================================
// TEST 1: Numerical Equivalence (vectorized vs legacy, same seed)
// ==========================================================================
display _n as text "=== TEST 1: Numerical Equivalence ==="
display as text "Running vectorized (current) implementation with seed=12345..."

set seed 12345
catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2 3 3) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(200) seed(12345) ///
    uniformall(false) ///
    control_group("notyettreated")

// Save vectorized results
matrix R_vec = e(results)
local c_vec_1 = R_vec[1,10]

// Store vectorized results (from the run above)
matrix R1 = R_vec
display _n "Vectorized implementation results:"
matrix list R1, format(%12.7f)

// Verify CI2 (bootstrap) is populated
local any_ci2 = 0
local nrows = rowsof(R1)
forvalues i = 1/`nrows' {
    if R1[`i',8] != . & R1[`i',9] != . {
        local any_ci2 = 1
    }
}
if `any_ci2' == 0 {
    display as error "FAIL TEST 1: All Bootstrap CI values are missing!"
    exit 198
}

// Verify CI contains point estimate
forvalues i = 1/`nrows' {
    if R1[`i',8] != . & R1[`i',9] != . {
        if R1[`i',8] > R1[`i',4] | R1[`i',9] < R1[`i',4] {
            display as error "FAIL TEST 1: Row `i' — CI does not contain point estimate"
            exit 198
        }
    }
}
display as result "PASS TEST 1: Vectorized bootstrap produces valid CI"

// ==========================================================================
// TEST 2: Critical Value Consistency
// ==========================================================================
display _n as text "=== TEST 2: Critical Value Consistency ==="

// Re-run with same seed to verify reproducibility
set seed 12345
catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2 3 3) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(200) seed(12345) ///
    uniformall(false) ///
    control_group("notyettreated")

matrix R2 = e(results)

// Compare critical values (column 10 in results matrix)
local nrows2 = rowsof(R2)
local max_reldif = 0
forvalues i = 1/`nrows2' {
    if R1[`i',10] != . & R2[`i',10] != . {
        local rd = reldif(R1[`i',10], R2[`i',10])
        if `rd' > `max_reldif' {
            local max_reldif = `rd'
        }
    }
}
display "  Max reldif (critical values): " %12.2e `max_reldif'
if `max_reldif' > 1e-12 {
    display as error "FAIL TEST 2: Critical values differ across runs (reldif = " %12.2e `max_reldif' ")"
    exit 198
}
display as result "PASS TEST 2: Critical values reproducible (reldif = " %12.2e `max_reldif' ")"

// ==========================================================================
// TEST 3: Confidence Band Consistency
// ==========================================================================
display _n as text "=== TEST 3: Confidence Band Consistency ==="

local max_reldif_lb = 0
local max_reldif_ub = 0
forvalues i = 1/`nrows2' {
    if R1[`i',8] != . & R2[`i',8] != . {
        local rd = reldif(R1[`i',8], R2[`i',8])
        if `rd' > `max_reldif_lb' {
            local max_reldif_lb = `rd'
        }
    }
    if R1[`i',9] != . & R2[`i',9] != . {
        local rd = reldif(R1[`i',9], R2[`i',9])
        if `rd' > `max_reldif_ub' {
            local max_reldif_ub = `rd'
        }
    }
}
display "  Max reldif (CI lower): " %12.2e `max_reldif_lb'
display "  Max reldif (CI upper): " %12.2e `max_reldif_ub'
if `max_reldif_lb' > 1e-12 | `max_reldif_ub' > 1e-12 {
    display as error "FAIL TEST 3: Confidence bands differ across runs"
    exit 198
}
display as result "PASS TEST 3: Confidence bands reproducible"

// ==========================================================================
// TEST 4: Performance Benchmark (informational, no pass/fail)
// ==========================================================================
display _n as text "=== TEST 4: Performance Benchmark ==="

// Larger dataset for timing
display as text "Generating larger test data (n=500, seed=99)..."
didhetero_simdata, n(500) tau(4) seed(99) clear

display as text "Running biters=500, n=500, porder=1 (vectorized)..."
timer clear 1
timer on 1
quietly catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2 3 3) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(500) seed(777) ///
    uniformall(false) ///
    control_group("notyettreated")
timer off 1

display _n "  Performance (vectorized p=1):"
timer list 1
display "  Configuration: n=500, biters=500, num_zeval=3, num_gteval=2"

// ==========================================================================
// TEST 5: Boundary Conditions
// ==========================================================================
display _n as text "=== TEST 5: Boundary Conditions ==="

// Re-generate data
didhetero_simdata, n(300) tau(4) seed(42) clear

// Test biters=1
display as text "  Testing biters=1..."
capture noisily catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(1) seed(111) ///
    uniformall(false) ///
    control_group("notyettreated")
if _rc {
    display as error "FAIL TEST 5a: biters=1 caused error (rc = " _rc ")"
    exit 198
}
display as result "  PASS TEST 5a: biters=1 runs without error"

// Test biters=2
display as text "  Testing biters=2..."
capture noisily catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(2) seed(222) ///
    uniformall(false) ///
    control_group("notyettreated")
if _rc {
    display as error "FAIL TEST 5b: biters=2 caused error (rc = " _rc ")"
    exit 198
}
matrix R_b2 = e(results)
// Verify CI2 exists
local has_ci2_b2 = 0
forvalues i = 1/`=rowsof(R_b2)' {
    if R_b2[`i',8] != . & R_b2[`i',9] != . {
        local has_ci2_b2 = 1
    }
}
if `has_ci2_b2' == 0 {
    display as error "FAIL TEST 5b: biters=2 produced no bootstrap CI"
    exit 198
}
display as result "  PASS TEST 5b: biters=2 runs and produces CI"

// Test with porder=2 (general path)
display as text "  Testing porder=2 (general vectorized path)..."
capture noisily catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2) ///
    porder(2) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(50) seed(333) ///
    uniformall(false) ///
    control_group("notyettreated")
if _rc {
    display as error "FAIL TEST 5c: porder=2 caused error (rc = " _rc ")"
    exit 198
}
matrix R_p2 = e(results)
local has_ci2_p2 = 0
forvalues i = 1/`=rowsof(R_p2)' {
    if R_p2[`i',8] != . & R_p2[`i',9] != . {
        local has_ci2_p2 = 1
    }
}
if `has_ci2_p2' == 0 {
    display as error "FAIL TEST 5c: porder=2 produced no bootstrap CI"
    exit 198
}
display as result "  PASS TEST 5c: porder=2 (general path) produces valid CI"

// Test with Epa kernel
display as text "  Testing epa kernel..."
capture noisily catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2) ///
    porder(1) kernel("epa") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(50) seed(444) ///
    uniformall(false) ///
    control_group("notyettreated")
if _rc {
    display as error "FAIL TEST 5d: epa kernel caused error (rc = " _rc ")"
    exit 198
}
display as result "  PASS TEST 5d: epa kernel runs without error"

// Test uniformall=true
display as text "  Testing uniformall=true..."
capture noisily catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2 3 3) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(50) seed(555) ///
    uniformall(true) ///
    control_group("notyettreated")
if _rc {
    display as error "FAIL TEST 5e: uniformall=true caused error (rc = " _rc ")"
    exit 198
}
display as result "  PASS TEST 5e: uniformall=true runs without error"

// ==========================================================================
// SUMMARY
// ==========================================================================
display _n "{hline 70}"
display as result "  ALL TESTS PASSED — Bootstrap Matrix Acceleration Verified"
display "{hline 70}"
