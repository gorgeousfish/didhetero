// =============================================================================
// Test GPS Diagnostics Matrix: e(gps_diagnostics)
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
// TEST 1: Matrix existence and dimensions (not-yet-treated control)
// =============================================================================
display _n "{hline 60}"
display "TEST 1: GPS diagnostics matrix dimensions (notyettreated)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap
}
if _rc == 0 {
    // Check matrix exists
    capture confirm matrix e(gps_diagnostics)
    if _rc == 0 {
        local _nrow = rowsof(e(gps_diagnostics))
        local _ncol = colsof(e(gps_diagnostics))
        local _num_gt = e(num_gteval)

        // Verify columns = 6
        if `_ncol' == 6 {
            display as result "  [PASS] gps_diagnostics has 6 columns"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] gps_diagnostics has `_ncol' columns (expected 6)"
            local _test_fail = `_test_fail' + 1
        }
        local _test_total = `_test_total' + 1

        // Verify rows = num_gteval (for notyettreated)
        if `_nrow' == `_num_gt' {
            display as result "  [PASS] gps_diagnostics rows (`_nrow') = num_gteval (`_num_gt')"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as result "  [FAIL] gps_diagnostics rows (`_nrow') != num_gteval (`_num_gt')"
            local _test_fail = `_test_fail' + 1
        }
        local _test_total = `_test_total' + 1
    }
    else {
        display as error "  [FAIL] e(gps_diagnostics) matrix does not exist"
        local _test_fail = `_test_fail' + 1
        local _test_total = `_test_total' + 1
    }
}
else {
    display as error "  [FAIL] catt_gt estimation failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
    local _test_total = `_test_total' + 1
}

// =============================================================================
// TEST 2: All models converge with standard simulated data
// =============================================================================
display _n "{hline 60}"
display "TEST 2: All GPS models converge (standard data)"
display "{hline 60}"

capture confirm matrix e(gps_diagnostics)
if _rc == 0 {
    matrix __diag = e(gps_diagnostics)
    local _nrow = rowsof(__diag)
    local _all_converged = 1
    forvalues i = 1/`_nrow' {
        if el(__diag, `i', 1) != 1 {
            local _all_converged = 0
            display as error "  Row `i': converged = " el(__diag, `i', 1)
        }
    }
    if `_all_converged' == 1 {
        display as result "  [PASS] All `_nrow' GPS models converged"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Some GPS models did not converge"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [SKIP] e(gps_diagnostics) not available"
    local _test_total = `_test_total' + 1
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 3: Iterations in reasonable range (1-25)
// =============================================================================
display _n "{hline 60}"
display "TEST 3: Iteration counts in valid range [1, 25]"
display "{hline 60}"

capture confirm matrix e(gps_diagnostics)
if _rc == 0 {
    matrix __diag = e(gps_diagnostics)
    local _nrow = rowsof(__diag)
    local _iters_ok = 1
    forvalues i = 1/`_nrow' {
        local _iter_i = el(__diag, `i', 2)
        if `_iter_i' < 1 | `_iter_i' > 25 {
            local _iters_ok = 0
            display as error "  Row `i': iterations = `_iter_i' (out of [1, 25])"
        }
    }
    if `_iters_ok' == 1 {
        display as result "  [PASS] All iterations in [1, 25]"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Some iterations outside valid range"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [SKIP] e(gps_diagnostics) not available"
    local _test_total = `_test_total' + 1
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 4: Extreme GPS count (n_extreme) is non-negative
// =============================================================================
display _n "{hline 60}"
display "TEST 4: n_extreme values are non-negative"
display "{hline 60}"

capture confirm matrix e(gps_diagnostics)
if _rc == 0 {
    matrix __diag = e(gps_diagnostics)
    local _nrow = rowsof(__diag)
    local _extreme_ok = 1
    forvalues i = 1/`_nrow' {
        local _ne_i = el(__diag, `i', 5)
        if `_ne_i' < 0 {
            local _extreme_ok = 0
            display as error "  Row `i': n_extreme = `_ne_i' (negative!)"
        }
    }
    if `_extreme_ok' == 1 {
        display as result "  [PASS] All n_extreme >= 0"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Negative n_extreme detected"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [SKIP] e(gps_diagnostics) not available"
    local _test_total = `_test_total' + 1
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 5: Condition number is positive and finite
// =============================================================================
display _n "{hline 60}"
display "TEST 5: Condition number positive and finite"
display "{hline 60}"

capture confirm matrix e(gps_diagnostics)
if _rc == 0 {
    matrix __diag = e(gps_diagnostics)
    local _nrow = rowsof(__diag)
    local _cond_ok = 1
    forvalues i = 1/`_nrow' {
        local _cn_i = el(__diag, `i', 6)
        if `_cn_i' <= 0 | missing(`_cn_i') {
            local _cond_ok = 0
            display as error "  Row `i': cond_number = `_cn_i'"
        }
    }
    if `_cond_ok' == 1 {
        display as result "  [PASS] All condition numbers positive and finite"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] Invalid condition number detected"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [SKIP] e(gps_diagnostics) not available"
    local _test_total = `_test_total' + 1
    local _test_fail = `_test_fail' + 1
}

// =============================================================================
// TEST 6: Matrix dimensions with never-treated control
// =============================================================================
display _n "{hline 60}"
display "TEST 6: GPS diagnostics matrix dimensions (nevertreated)"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(nevertreated) nobstrap
}
if _rc == 0 {
    capture confirm matrix e(gps_diagnostics)
    if _rc == 0 {
        local _nrow = rowsof(e(gps_diagnostics))
        local _ncol = colsof(e(gps_diagnostics))

        // For nevertreated: rows = number of unique treatment groups
        if `_ncol' == 6 {
            display as result "  [PASS] gps_diagnostics has 6 columns (nevertreated)"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] gps_diagnostics has `_ncol' columns (expected 6)"
            local _test_fail = `_test_fail' + 1
        }
        local _test_total = `_test_total' + 1

        // Rows should be > 0 (number of treatment groups)
        if `_nrow' > 0 {
            display as result "  [PASS] gps_diagnostics has `_nrow' rows (nevertreated groups)"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] gps_diagnostics has 0 rows"
            local _test_fail = `_test_fail' + 1
        }
        local _test_total = `_test_total' + 1
    }
    else {
        display as error "  [FAIL] e(gps_diagnostics) matrix does not exist"
        local _test_fail = `_test_fail' + 1
        local _test_total = `_test_total' + 1
    }
}
else {
    display as error "  [FAIL] catt_gt estimation failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
    local _test_total = `_test_total' + 1
}

// =============================================================================
// TEST 7: GPS diagnostics preserved after aggte_gt
// =============================================================================
display _n "{hline 60}"
display "TEST 7: GPS diagnostics inherited by aggte_gt"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap
    // Save diagnostics from catt_gt
    matrix __diag_catt = e(gps_diagnostics)
    local _nrow_catt = rowsof(__diag_catt)

    // Run aggte
    aggte_gt, type("dynamic") bstrap("false")
}
if _rc == 0 {
    capture confirm matrix e(gps_diagnostics)
    if _rc == 0 {
        local _nrow_agg = rowsof(e(gps_diagnostics))
        if `_nrow_agg' == `_nrow_catt' {
            display as result "  [PASS] aggte_gt inherits gps_diagnostics (`_nrow_agg' rows)"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] aggte rows = `_nrow_agg', catt rows = `_nrow_catt'"
            local _test_fail = `_test_fail' + 1
        }
    }
    else {
        display as error "  [FAIL] e(gps_diagnostics) not found after aggte_gt"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [FAIL] aggte_gt failed (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
    local _test_total = `_test_total' + 1
}

// =============================================================================
// TEST 8: gpsstrict mode interaction
// When GPS models converge (standard data), gpsstrict should not error.
// =============================================================================
display _n "{hline 60}"
display "TEST 8: gpsstrict mode — no error when all converge"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        control_group(notyettreated) nobstrap gpsstrict
}
if _rc == 0 {
    capture confirm matrix e(gps_diagnostics)
    if _rc == 0 {
        matrix __diag = e(gps_diagnostics)
        local _nrow = rowsof(__diag)
        local _all_converged = 1
        forvalues i = 1/`_nrow' {
            if el(__diag, `i', 1) != 1 {
                local _all_converged = 0
            }
        }
        if `_all_converged' == 1 {
            display as result "  [PASS] gpsstrict + all converged = no error"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] gpsstrict: some models reported non-convergence"
            local _test_fail = `_test_fail' + 1
        }
    }
    else {
        display as error "  [FAIL] e(gps_diagnostics) not available with gpsstrict"
        local _test_fail = `_test_fail' + 1
    }
    local _test_total = `_test_total' + 1
}
else {
    display as error "  [FAIL] gpsstrict errored even with standard data (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
    local _test_total = `_test_total' + 1
}

// =============================================================================
// TEST 9: Log-likelihood values are finite and negative
// =============================================================================
display _n "{hline 60}"
display "TEST 9: Log-likelihood values finite and negative"
display "{hline 60}"

capture noisily {
    didhetero_simdata, n(500) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(1) bwselect(IMSE1) ///
        nobstrap
}
if _rc == 0 {
    capture confirm matrix e(gps_diagnostics)
    if _rc == 0 {
        matrix __diag = e(gps_diagnostics)
        local _nrow = rowsof(__diag)
        local _ll_ok = 1
        forvalues i = 1/`_nrow' {
            local _ll_i = el(__diag, `i', 4)
            if `_ll_i' >= 0 | missing(`_ll_i') {
                local _ll_ok = 0
                display as error "  Row `i': ll_final = `_ll_i' (expected < 0)"
            }
        }
        if `_ll_ok' == 1 {
            display as result "  [PASS] All log-likelihood values < 0 and finite"
            local _test_pass = `_test_pass' + 1
        }
        else {
            display as error "  [FAIL] Invalid log-likelihood values detected"
            local _test_fail = `_test_fail' + 1
        }
        local _test_total = `_test_total' + 1
    }
    else {
        display as error "  [SKIP] e(gps_diagnostics) not available"
        local _test_total = `_test_total' + 1
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] catt_gt estimation failed"
    local _test_fail = `_test_fail' + 1
    local _test_total = `_test_total' + 1
}

// =============================================================================
// TEST 10: Display full diagnostics matrix for inspection
// =============================================================================
display _n "{hline 60}"
display "TEST 10: Display diagnostics matrix (visual inspection)"
display "{hline 60}"

capture confirm matrix e(gps_diagnostics)
if _rc == 0 {
    matrix list e(gps_diagnostics), format(%12.4f)
    display as result "  [PASS] Diagnostics matrix displayed successfully"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Cannot display e(gps_diagnostics)"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "GPS DIAGNOSTICS TEST SUMMARY"
display "{hline 60}"
display as result "  Total:  `_test_total'"
display as result "  Passed: `_test_pass'"
if `_test_fail' > 0 {
    display as error "  Failed: `_test_fail'"
}
else {
    display as result "  Failed: 0"
}
display "{hline 60}"

if `_test_fail' > 0 {
    exit 198
}
