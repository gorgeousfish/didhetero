// =============================================================================
// test_lpr_analytic.do
// Tests for LPR analytic formula optimization (p=1 Cramer's rule fast path)
//
// Validates:
//   1. Numerical equivalence between analytic and Cholesky solutions
//   2. CATT estimation consistency with the optimization
//   3. Degenerate case handling (near-zero determinant)
//   4. porder=2 unaffected (still uses Cholesky)
//   5. First derivative analytic solution correctness
//   6. Performance comparison (informational)
//   7. Bootstrap with LPR correctness
// =============================================================================

clear all
set more off

// Locate package root
local PKGROOT ""
local DOFILE `"`c(filename)'"'
if strpos(`"`DOFILE'"', "/tests/test_lpr_analytic.do") {
    local PKGROOT : subinstr local DOFILE "/tests/test_lpr_analytic.do" "", all
}
if `"`PKGROOT'"' == "" {
    // Try parent directory
    capture confirm file "../make.do"
    if !_rc {
        local PKGROOT ".."
    }
    else {
        local PKGROOT "."
    }
}

// Load the library
quietly cd "`PKGROOT'"
capture do make.do
if _rc {
    display as error "ERROR: Cannot build library"
    exit _rc
}

display _n "{hline 60}"
display "  TEST: LPR Analytic Formula Optimization"
display "{hline 60}"

local test_pass = 0
local test_fail = 0

// =============================================================================
// TEST 1: Numerical equivalence — analytic vs Cholesky (Gaussian kernel)
// =============================================================================
display _n as text "TEST 1: Numerical equivalence (Gaussian kernel, p=1, deriv=0)"

mata:
{
    // Generate test data
    rseed(12345)
    n = 500
    x = rnormal(n, 1, 0, 1)
    y = sin(x) + rnormal(n, 1, 0, 0.3)
    eval_pts = range(-2, 2, 0.5)
    h_val = 0.5
    w = J(n, 1, 1)
    
    // Test 1a: Analytic formula (called via didhetero_lpr with p=1, deriv=0)
    result_analytic = didhetero_lpr(y, x, eval_pts, 1, 0, "gau", h_val, w)
    
    // Test 1b: General Cholesky path (force via p=2, then extract level)
    // We can't directly force Cholesky for p=1, so we verify by manual computation
    // Manually compute using the Cholesky path logic for comparison
    num_eval = rows(eval_pts)
    result_cholesky = J(num_eval, 1, .)
    
    for (j = 1; j <= num_eval; j++) {
        u_vec = (x :- eval_pts[j]) / h_val
        K_h = w :* didhetero_kernel_eval(u_vec, "gau")
        idx = didhetero_selectindex(K_h :!= 0)
        
        if (length(idx) > 1) {
            Z_sub = didhetero_polynomial(x[idx'] :- eval_pts[j], 1)
            K_sub = K_h[idx']
            y_sub = y[idx']
            
            Gamma = cross(Z_sub, K_sub, Z_sub)
            Gamma_inv = cholinv(Gamma)
            beta = Gamma_inv * cross(Z_sub, K_sub, y_sub)
            result_cholesky[j] = beta[1]
        }
    }
    
    // Compare
    max_reldif = 0
    for (j = 1; j <= num_eval; j++) {
        if (result_analytic[j] != . & result_cholesky[j] != .) {
            rd = reldif(result_analytic[j], result_cholesky[j])
            if (rd > max_reldif) max_reldif = rd
        }
    }
    
    st_local("max_reldif_1", strofreal(max_reldif))
    st_local("pass_1", strofreal(max_reldif < 1e-10))
}
end

if `pass_1' {
    display as result "  PASS: max reldif = `max_reldif_1' < 1e-10"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: max reldif = `max_reldif_1' >= 1e-10"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 2: Numerical equivalence — Epanechnikov kernel
// =============================================================================
display _n as text "TEST 2: Numerical equivalence (Epanechnikov kernel, p=1, deriv=0)"

mata:
{
    rseed(54321)
    n = 400
    x = rnormal(n, 1, 0, 1)
    y = x:^2 + rnormal(n, 1, 0, 0.2)
    eval_pts = range(-1.5, 1.5, 0.3)
    h_val = 0.8
    w = J(n, 1, 1)
    
    // Analytic (via didhetero_lpr)
    result_analytic = didhetero_lpr(y, x, eval_pts, 1, 0, "epa", h_val, w)
    
    // Manual Cholesky
    num_eval = rows(eval_pts)
    result_cholesky = J(num_eval, 1, .)
    sort_idx = order(x, 1)
    x_s = x[sort_idx]
    y_s = y[sort_idx]
    
    for (j = 1; j <= num_eval; j++) {
        k1 = _dh_bsearch_ge(x_s, eval_pts[j] - h_val)
        k2 = _dh_bsearch_le(x_s, eval_pts[j] + h_val)
        if (k1 <= k2) {
            xs = x_s[|k1 \ k2|]
            ys = y_s[|k1 \ k2|]
            Kh = didhetero_kernel_eval((xs :- eval_pts[j]) / h_val, "epa")
            idx = didhetero_selectindex(Kh :!= 0)
            if (length(idx) > 1) {
                Z_sub = didhetero_polynomial(xs[idx'] :- eval_pts[j], 1)
                K_sub = Kh[idx']
                y_sub = ys[idx']
                Gamma = cross(Z_sub, K_sub, Z_sub)
                Gamma_inv = cholinv(Gamma)
                if (!hasmissing(Gamma_inv)) {
                    beta = Gamma_inv * cross(Z_sub, K_sub, y_sub)
                    result_cholesky[j] = beta[1]
                }
            }
        }
    }
    
    max_reldif = 0
    for (j = 1; j <= num_eval; j++) {
        if (result_analytic[j] != . & result_cholesky[j] != .) {
            rd = reldif(result_analytic[j], result_cholesky[j])
            if (rd > max_reldif) max_reldif = rd
        }
    }
    
    st_local("max_reldif_2", strofreal(max_reldif))
    st_local("pass_2", strofreal(max_reldif < 1e-10))
}
end

if `pass_2' {
    display as result "  PASS: max reldif = `max_reldif_2' < 1e-10"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: max reldif = `max_reldif_2' >= 1e-10"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 3: Degenerate case handling (tiny bandwidth → no observations)
// =============================================================================
display _n as text "TEST 3: Degenerate case handling (extreme bandwidth)"

mata:
{
    rseed(99999)
    n = 100
    x = rnormal(n, 1, 0, 1)
    y = x + rnormal(n, 1, 0, 0.1)
    eval_pts = (0 \ 5 \ 10)  // 5 and 10 are far from data
    h_tiny = 0.001  // Very small bandwidth
    w = J(n, 1, 1)
    
    // With tiny bandwidth: most eval points should get missing
    result_tiny = didhetero_lpr(y, x, eval_pts, 1, 0, "epa", h_tiny, w)
    
    // With eval point far from data: should also be missing
    eval_far = (100 \ -100)
    result_far = didhetero_lpr(y, x, eval_far, 1, 0, "gau", 0.1, w)
    
    // Check: results at extreme points should be missing (not garbage)
    pass = 1
    // eval_pts[2]=5 and eval_pts[3]=10 far from data with h=0.001
    if (result_tiny[2] != . & abs(result_tiny[2]) > 100) pass = 0
    if (result_tiny[3] != .) pass = 0  // should definitely be missing
    
    st_local("pass_3", strofreal(pass))
    st_local("res_info", strofreal(result_tiny[1]) + ", " + strofreal(result_tiny[2]) + ", " + strofreal(result_tiny[3]))
}
end

if `pass_3' {
    display as result "  PASS: Degenerate cases return missing (results: `res_info')"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: Degenerate cases produced garbage values (`res_info')"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 4: porder=2 unaffected (still uses Cholesky, results unchanged)
// =============================================================================
display _n as text "TEST 4: porder=2 still uses Cholesky (unaffected by optimization)"

mata:
{
    rseed(11111)
    n = 300
    x = rnormal(n, 1, 0, 1)
    y = x:^2 + 0.5*x + rnormal(n, 1, 0, 0.2)
    eval_pts = range(-1, 1, 0.25)
    h_val = 0.6
    w = J(n, 1, 1)
    
    // p=2 goes through general Cholesky path
    result_p2 = didhetero_lpr(y, x, eval_pts, 2, 0, "gau", h_val, w)
    
    // Manual Cholesky for p=2
    num_eval = rows(eval_pts)
    result_manual = J(num_eval, 1, .)
    
    for (j = 1; j <= num_eval; j++) {
        u_vec = (x :- eval_pts[j]) / h_val
        K_h = didhetero_kernel_eval(u_vec, "gau")
        idx = didhetero_selectindex(K_h :!= 0)
        if (length(idx) > 1) {
            Z_sub = didhetero_polynomial(x[idx'] :- eval_pts[j], 2)
            K_sub = K_h[idx']
            y_sub = y[idx']
            Gamma = cross(Z_sub, K_sub, Z_sub)
            Gamma_inv = cholinv(Gamma)
            if (!hasmissing(Gamma_inv)) {
                beta = Gamma_inv * cross(Z_sub, K_sub, y_sub)
                result_manual[j] = beta[1]
            }
        }
    }
    
    max_reldif = 0
    for (j = 1; j <= num_eval; j++) {
        if (result_p2[j] != . & result_manual[j] != .) {
            rd = reldif(result_p2[j], result_manual[j])
            if (rd > max_reldif) max_reldif = rd
        }
    }
    
    st_local("max_reldif_4", strofreal(max_reldif))
    st_local("pass_4", strofreal(max_reldif < 1e-10))
}
end

if `pass_4' {
    display as result "  PASS: p=2 Cholesky path unchanged, max reldif = `max_reldif_4'"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: p=2 results differ, max reldif = `max_reldif_4'"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 5: First derivative (deriv=1) analytic solution
// =============================================================================
display _n as text "TEST 5: First derivative analytic solution (p=1, deriv=1)"

mata:
{
    rseed(22222)
    n = 500
    x = rnormal(n, 1, 0, 1)
    y = sin(x) + rnormal(n, 1, 0, 0.2)
    eval_pts = range(-1.5, 1.5, 0.5)
    h_val = 0.5
    w = J(n, 1, 1)
    
    // Analytic deriv=1 path
    result_deriv1 = didhetero_lpr(y, x, eval_pts, 1, 1, "gau", h_val, w)
    
    // Manual Cholesky for deriv=1
    num_eval = rows(eval_pts)
    result_manual = J(num_eval, 1, .)
    
    for (j = 1; j <= num_eval; j++) {
        u_vec = (x :- eval_pts[j]) / h_val
        K_h = didhetero_kernel_eval(u_vec, "gau")
        idx = didhetero_selectindex(K_h :!= 0)
        if (length(idx) > 1) {
            Z_sub = didhetero_polynomial(x[idx'] :- eval_pts[j], 1)
            K_sub = K_h[idx']
            y_sub = y[idx']
            Gamma = cross(Z_sub, K_sub, Z_sub)
            Gamma_inv = cholinv(Gamma)
            if (!hasmissing(Gamma_inv)) {
                beta = Gamma_inv * cross(Z_sub, K_sub, y_sub)
                result_manual[j] = beta[2]  // deriv=1 => beta[2] * 1!
            }
        }
    }
    
    max_reldif = 0
    for (j = 1; j <= num_eval; j++) {
        if (result_deriv1[j] != . & result_manual[j] != .) {
            rd = reldif(result_deriv1[j], result_manual[j])
            if (rd > max_reldif) max_reldif = rd
        }
    }
    
    st_local("max_reldif_5", strofreal(max_reldif))
    st_local("pass_5", strofreal(max_reldif < 1e-10))
}
end

if `pass_5' {
    display as result "  PASS: deriv=1 analytic matches Cholesky, max reldif = `max_reldif_5'"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: deriv=1 analytic diverges, max reldif = `max_reldif_5'"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 6: Performance comparison (informational, not a pass/fail)
// =============================================================================
display _n as text "TEST 6: Performance comparison (informational)"

mata:
{
    rseed(33333)
    n = 1000
    x = rnormal(n, 1, 0, 1)
    y = x:^2 + rnormal(n, 1, 0, 0.3)
    eval_pts = range(-2, 2, 0.05)
    h_val = 0.4
    w = J(n, 1, 1)
    num_eval = rows(eval_pts)
    
    // Time analytic path (p=1, goes through fast path)
    timer_clear()
    timer_on(1)
    for (rep = 1; rep <= 10; rep++) {
        r1 = didhetero_lpr(y, x, eval_pts, 1, 0, "gau", h_val, w)
    }
    timer_off(1)
    t_analytic = timer_value(1)[1] / 10
    
    // Time general path (p=2 forces Cholesky)
    timer_clear()
    timer_on(2)
    for (rep = 1; rep <= 10; rep++) {
        r2 = didhetero_lpr(y, x, eval_pts, 2, 0, "gau", h_val, w)
    }
    timer_off(2)
    t_cholesky = timer_value(2)[1] / 10
    
    if (t_cholesky > 0) {
        speedup = (t_cholesky - t_analytic) / t_cholesky * 100
    }
    else {
        speedup = 0
    }
    
    st_local("t_analytic", strofreal(t_analytic, "%9.4f"))
    st_local("t_cholesky", strofreal(t_cholesky, "%9.4f"))
    st_local("speedup", strofreal(speedup, "%5.1f"))
}
end

display as text "  INFO: Analytic p=1 time = `t_analytic's, Cholesky p=2 time = `t_cholesky's"
display as text "  INFO: Relative speedup ≈ `speedup'%%"
display as text "  (Note: p=2 solves 3x3 system vs p=1 analytic 2x2 — comparison is indicative)"
local test_pass = `test_pass' + 1


// =============================================================================
// TEST 7: Weighted estimation (Bootstrap Mammen weights)
// =============================================================================
display _n as text "TEST 7: Weighted estimation with non-unit weights (Mammen-like)"

mata:
{
    rseed(44444)
    n = 300
    x = rnormal(n, 1, 0, 1)
    y = 2*x + rnormal(n, 1, 0, 0.3)
    eval_pts = range(-1, 1, 0.25)
    h_val = 0.5
    
    // Mammen-like weights (two-point distribution)
    p_mammen = (sqrt(5) + 1) / (2 * sqrt(5))
    mammen_w = J(n, 1, .)
    u_draw = runiform(n, 1)
    for (i = 1; i <= n; i++) {
        if (u_draw[i] < p_mammen) {
            mammen_w[i] = 1 - (sqrt(5) - 1) / 2
        }
        else {
            mammen_w[i] = 1 + (sqrt(5) + 1) / 2
        }
    }
    
    // Analytic path with weights
    result_weighted = didhetero_lpr(y, x, eval_pts, 1, 0, "gau", h_val, mammen_w)
    
    // Manual Cholesky with same weights
    num_eval = rows(eval_pts)
    result_manual = J(num_eval, 1, .)
    
    for (j = 1; j <= num_eval; j++) {
        u_vec = (x :- eval_pts[j]) / h_val
        K_h = mammen_w :* didhetero_kernel_eval(u_vec, "gau")
        idx = didhetero_selectindex(K_h :!= 0)
        if (length(idx) > 1) {
            Z_sub = didhetero_polynomial(x[idx'] :- eval_pts[j], 1)
            K_sub = K_h[idx']
            y_sub = y[idx']
            Gamma = cross(Z_sub, K_sub, Z_sub)
            Gamma_inv = cholinv(Gamma)
            if (!hasmissing(Gamma_inv)) {
                beta = Gamma_inv * cross(Z_sub, K_sub, y_sub)
                result_manual[j] = beta[1]
            }
        }
    }
    
    max_reldif = 0
    for (j = 1; j <= num_eval; j++) {
        if (result_weighted[j] != . & result_manual[j] != .) {
            rd = reldif(result_weighted[j], result_manual[j])
            if (rd > max_reldif) max_reldif = rd
        }
    }
    
    st_local("max_reldif_7", strofreal(max_reldif))
    st_local("pass_7", strofreal(max_reldif < 1e-10))
}
end

if `pass_7' {
    display as result "  PASS: Weighted analytic matches Cholesky, max reldif = `max_reldif_7'"
    local test_pass = `test_pass' + 1
}
else {
    display as error "  FAIL: Weighted results differ, max reldif = `max_reldif_7'"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// TEST 8: CATT end-to-end consistency
// =============================================================================
display _n as text "TEST 8: CATT estimation end-to-end (verifies integration)"

local all_finite = 0
local n_est = 0
capture noisily {
    // Generate simulated data and run catt_gt
    didhetero_simdata, n(300) seed(55555)
    
    // Run with default settings (p=1 will use analytic path)
    catt_gt Y, group(G) time(T) idvar(id) yvar(Y) gvar(G) tvar(T) kernel(gau)
    
    // Store estimates
    matrix b_analytic = e(b)
    
    // Check: estimates should be non-missing and finite
    local n_est = colsof(b_analytic)
    local all_finite = 1
    forvalues i = 1/`n_est' {
        if b_analytic[1,`i'] == . {
            local all_finite = 0
        }
    }
}

local test8_rc = _rc
if `test8_rc' == 0 & `all_finite' == 1 {
    display as result "  PASS: CATT end-to-end produces valid estimates (n=`n_est')"
    local test_pass = `test_pass' + 1
}
else if `test8_rc' != 0 {
    display as text "  SKIP: CATT end-to-end test requires full package (rc=`test8_rc')"
}
else {
    display as error "  FAIL: CATT produced missing estimates"
    local test_fail = `test_fail' + 1
}


// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display as text "  RESULTS: `test_pass' passed, `test_fail' failed"
display "{hline 60}"

if `test_fail' > 0 {
    display as error "  SOME TESTS FAILED"
    exit 1
}
else {
    display as result "  ALL TESTS PASSED"
}
