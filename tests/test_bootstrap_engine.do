// =============================================================================
// test_bootstrap_engine.do
// Unit tests for bootstrap engine infrastructure (Task #15)
//
// Tests:
//   1. Mammen weights: E[V*]=1 and Var[V*]=1
//   2. Quantile function accuracy (known distribution)
//   3. Sup operation: uniformall=TRUE vs FALSE behavior
//   4. Mammen constants correctness
//   5. Weight generation: dimension and method dispatch
// =============================================================================
clear all
set more off

// Load package
quietly adopath ++ "/Users/cxy/Desktop/didhetero/didhetero-main/ado"
quietly mata: mata mlib index

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

display _n "{hline 70}"
display "  Bootstrap Engine — Infrastructure Unit Tests"
display "{hline 70}"

// ==========================================================================
// TEST 1: Mammen Weights Mean=1, Variance=1
// ==========================================================================
display _n as text "=== TEST 1: Mammen Weights E[V*]=1, Var[V*]=1 ==="
local _test_total = `_test_total' + 1

mata:
    rseed(12345)
    _test_W = dh_boot_generate_weights(10000, 500, "mammen")
    _test_ok = (rows(_test_W) == 10000) & (cols(_test_W) == 500)
    _test_mean = mean(vec(_test_W))
    _test_var = variance(vec(_test_W))
    printf("  Dimensions: %g x %g (expected 10000 x 500)\n", rows(_test_W), cols(_test_W))
    printf("  Mean of weights: %12.8f (expected 1.0)\n", _test_mean)
    printf("  Var  of weights: %12.8f (expected 1.0)\n", _test_var)
    _test_ok = _test_ok & (abs(_test_mean - 1) < 0.02) & (abs(_test_var - 1) < 0.02)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_W _test_ok _test_mean _test_var
end

if "`_result'" == "PASS" {
    display as result "  PASS: Mammen weights satisfy E[V*]=1, Var[V*]=1"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Mammen weights do not satisfy expected moments"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 2: Quantile Function Accuracy (Known Distribution)
// ==========================================================================
display _n as text "=== TEST 2: Quantile Function Accuracy ==="
local _test_total = `_test_total' + 1

mata:
    _test_x = (1::1000)
    _test_q50 = dh_boot_quantile(_test_x, 0.5)
    _test_q25 = dh_boot_quantile(_test_x, 0.25)
    _test_q75 = dh_boot_quantile(_test_x, 0.75)
    _test_q99 = dh_boot_quantile(_test_x, 0.99)
    printf("  q(0.25) = %12.4f (expected ~250.75)\n", _test_q25)
    printf("  q(0.50) = %12.4f (expected ~500.50)\n", _test_q50)
    printf("  q(0.75) = %12.4f (expected ~750.25)\n", _test_q75)
    printf("  q(0.99) = %12.4f (expected ~990.01)\n", _test_q99)
    _test_ok = (abs(_test_q50 - 500.5) < 0.5) & (abs(_test_q25 - 250.75) < 0.5) & ///
        (abs(_test_q75 - 750.25) < 0.5) & (abs(_test_q99 - 990.01) < 0.5)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_x _test_q50 _test_q25 _test_q75 _test_q99 _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Quantile function matches expected values"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Quantile function produced incorrect values"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 2b: Quantile Function with Missing Values
// ==========================================================================
display _n as text "=== TEST 2b: Quantile Function with Missing Values ==="
local _test_total = `_test_total' + 1

mata:
    _test_x = (1 \ . \ 3 \ . \ 5 \ . \ 7 \ . \ 9)
    _test_q = dh_boot_quantile(_test_x, 0.5)
    _test_qm = dh_boot_quantile(J(10, 1, .), 0.5)
    printf("  q(0.5) with missing: %12.4f (expected 5.0)\n", _test_q)
    printf("  q(0.5) all missing:  %12.4f (expected .)\n", _test_qm)
    _test_ok = (_test_q == 5) & (_test_qm >= .)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_x _test_q _test_qm _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Quantile function handles missing values correctly"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Quantile function does not handle missing correctly"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 3: Sup Operation — uniformall=TRUE vs FALSE
// ==========================================================================
display _n as text "=== TEST 3: Sup Operation (uniformall behavior) ==="
local _test_total = `_test_total' + 1

mata:
    rseed(999)
    _test_T = abs(rnormal(100, 3, 0, 1))
    _test_cu = dh_boot_compute_sup(_test_T, 1, 0.05)
    _test_cp = dh_boot_compute_sup(_test_T, 0, 0.05)
    printf("  uniformall=1: c = (%9.4f, %9.4f, %9.4f)\n", ///
        _test_cu[1], _test_cu[2], _test_cu[3])
    printf("  uniformall=0: c = (%9.4f, %9.4f, %9.4f)\n", ///
        _test_cp[1], _test_cp[2], _test_cp[3])
    _test_ok = (_test_cu[1] == _test_cu[2]) & (_test_cu[2] == _test_cu[3])
    _test_ok = _test_ok & (_test_cu[1] >= max(_test_cp) - 1e-10)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_T _test_cu _test_cp _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Sup operation correct for uniformall TRUE/FALSE"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Sup operation inconsistent"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 4: Mammen Constants Correctness
// ==========================================================================
display _n as text "=== TEST 4: Mammen Distribution Constants ==="
local _test_total = `_test_total' + 1

mata:
    _test_kappa = .
    _test_prob = .
    _test_vlo = .
    _test_vhi = .
    dh_boot_mammen_params(_test_kappa, _test_prob, _test_vlo, _test_vhi)
    _test_EV = _test_prob * _test_vlo + (1 - _test_prob) * _test_vhi
    _test_VarV = _test_prob * _test_vlo^2 + (1 - _test_prob) * _test_vhi^2 - _test_EV^2
    printf("  kappa   = %18.15f (golden ratio)\n", _test_kappa)
    printf("  v_lo    = %18.15f\n", _test_vlo)
    printf("  v_hi    = %18.15f\n", _test_vhi)
    printf("  prob_lo = %18.15f\n", _test_prob)
    printf("  E[V*]   = %18.15f (expected 1.0)\n", _test_EV)
    printf("  Var[V*] = %18.15f (expected 1.0)\n", _test_VarV)
    _test_ok = (abs(_test_kappa - (sqrt(5) + 1) / 2) < 1e-14) & ///
        (abs(_test_EV - 1) < 1e-14) & (abs(_test_VarV - 1) < 1e-14)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_kappa _test_prob _test_vlo _test_vhi _test_EV _test_VarV _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Mammen constants satisfy E=1, Var=1 exactly"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Mammen constants incorrect"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 5: Weight Generation — Rademacher Method
// ==========================================================================
display _n as text "=== TEST 5: Weight Generation — Rademacher Method ==="
local _test_total = `_test_total' + 1

mata:
    rseed(54321)
    _test_W = dh_boot_generate_weights(5000, 200, "rademacher")
    _test_mean = mean(vec(_test_W))
    _test_var = variance(vec(_test_W))
    _test_binary = (sum(_test_W :== 0) + sum(_test_W :== 2) == rows(_test_W) * cols(_test_W))
    printf("  Mean of Rademacher weights: %12.8f (expected 1.0)\n", _test_mean)
    printf("  Var  of Rademacher weights: %12.8f (expected 1.0)\n", _test_var)
    printf("  All values in {0, 2}:       %g\n", _test_binary)
    _test_ok = (abs(_test_mean - 1) < 0.03) & (abs(_test_var - 1) < 0.03) & _test_binary
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_W _test_mean _test_var _test_binary _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Rademacher weights correct (E=1, Var=1, binary)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Rademacher weights incorrect"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 6: Backward Compatibility — Legacy API wrappers
// ==========================================================================
display _n as text "=== TEST 6: Legacy API Backward Compatibility ==="
local _test_total = `_test_total' + 1

mata:
    rseed(111)
    _test_w1 = _didhetero_mammen_weights(1000)
    rseed(222)
    _test_wb = _didhetero_mammen_weights_batch(50, 200)
    _test_m1 = mean(_test_w1)
    _test_mb = mean(vec(_test_wb))
    printf("  _didhetero_mammen_weights(1000): %g x %g, mean=%9.6f\n", ///
        rows(_test_w1), cols(_test_w1), _test_m1)
    printf("  _didhetero_mammen_weights_batch(50,200): %g x %g, mean=%9.6f\n", ///
        rows(_test_wb), cols(_test_wb), _test_mb)
    _test_ok = (rows(_test_w1) == 1000) & (cols(_test_w1) == 1) & ///
        (rows(_test_wb) == 50) & (cols(_test_wb) == 200) & ///
        (abs(_test_m1 - 1) < 0.1) & (abs(_test_mb - 1) < 0.1)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_w1 _test_wb _test_m1 _test_mb _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: Legacy API wrappers maintain correct signatures"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: Legacy API wrappers broken"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// TEST 7: Max Non-Missing Function
// ==========================================================================
display _n as text "=== TEST 7: dh_boot_max_nonmissing() ==="
local _test_total = `_test_total' + 1

mata:
    _test_m1 = dh_boot_max_nonmissing((1 \ 5 \ 3 \ . \ 2))
    _test_m2 = dh_boot_max_nonmissing(J(5, 1, .))
    _test_m3 = dh_boot_max_nonmissing((42))
    printf("  max(1,5,3,.,2) = %g (expected 5)\n", _test_m1)
    printf("  max(all .)     = %g (expected .)\n", _test_m2)
    printf("  max(42)        = %g (expected 42)\n", _test_m3)
    _test_ok = (_test_m1 == 5) & (_test_m2 >= .) & (_test_m3 == 42)
    st_local("_result", _test_ok ? "PASS" : "FAIL")
    mata drop _test_m1 _test_m2 _test_m3 _test_ok
end

if "`_result'" == "PASS" {
    display as result "  PASS: dh_boot_max_nonmissing() correct"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  FAIL: dh_boot_max_nonmissing() incorrect"
    local _test_fail = `_test_fail' + 1
}

// ==========================================================================
// SUMMARY
// ==========================================================================
display _n "{hline 70}"
display as text "  Bootstrap Engine Test Summary"
display "{hline 70}"
display as text "  Total:  `_test_total'"
display as result "  Passed: `_test_pass'"
if `_test_fail' > 0 {
    display as error "  Failed: `_test_fail'"
}
else {
    display as text "  Failed: 0"
}
display "{hline 70}"

if `_test_fail' > 0 {
    display as error _n "SOME TESTS FAILED"
    exit 9
}
else {
    display as result _n "ALL TESTS PASSED"
}
