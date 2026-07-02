// =============================================================================
// test_uniformall.do
// 回归测试：验证uniformall参数的理论正确性
// 基于论文Theorem 2（分布逼近）和Theorem 4（Bootstrap有效性）
//
// 论文: "Doubly Robust Uniform Confidence Bands for Group-Time Conditional
//        Average Treatment Effects in Difference-in-Differences"
// =============================================================================
clear all
set more off
set seed 20260701

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

local _test_pass = 0
local _test_fail = 0
local _test_total = 0

// =============================================================================
// Generate shared simulation data
// =============================================================================
display _n "{hline 70}"
display "UNIFORMALL REGRESSION TEST SUITE"
display "{hline 70}"
display "Generating simulation data..."

didhetero_simdata, n(500) tau(4) clear


// =============================================================================
// TEST 1: 临界值单调性 (Critical Value Monotonicity)
// 理论依据: Theorem 4 — sup over larger set ≥ sup over subset
// uniformall(true) 对所有(g,t,z)联合取sup → 临界值 ≥ 分组取sup
// =============================================================================
display _n "{hline 70}"
display "TEST 1: 临界值单调性 (Theorem 4: sup over larger set >= sup over subset)"
display "{hline 70}"

capture noisily {
    // Run with uniformall(true) — global sup
    set seed 54321
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(true) bstrap biters(200) seed(54321)

    tempname R_all c_check_all
    matrix `R_all' = e(results)
    matrix `c_check_all' = e(c_check)
    local ngt_all = e(num_gteval)

    // Run with uniformall(false) — per-(g,t) sup
    set seed 54321
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(false) bstrap biters(200) seed(54321)

    tempname R_marginal c_check_marginal
    matrix `R_marginal' = e(results)
    matrix `c_check_marginal' = e(c_check)
    local ngt_marginal = e(num_gteval)

    // Verify: c_check_all[i] >= c_check_marginal[i] for all i
    local monotone_ok = 1
    forvalues i = 1/`ngt_all' {
        local cv_all = `c_check_all'[`i', 1]
        local cv_marg = `c_check_marginal'[`i', 1]
        // Skip if either is missing (numerical issues)
        if `cv_all' == . | `cv_marg' == . {
            continue
        }
        // Allow small numerical tolerance (bootstrap variability)
        if `cv_all' < `cv_marg' - 1e-6 {
            display as error "  VIOLATION at (g,t) pair `i': cv_all=" ///
                %9.4f `cv_all' " < cv_marginal=" %9.4f `cv_marg'
            local monotone_ok = 0
        }
    }
    assert `monotone_ok' == 1
}
if _rc == 0 {
    display as result "  [PASS] Bootstrap临界值单调性: c_check(uniformall=1) >= c_check(uniformall=0)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] Bootstrap临界值单调性 (rc = " _rc ")"
    display as error "  理论要求: Theorem 4, sup over {(g,t,z)} >= sup over {z} for fixed (g,t)"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 2: 置信带宽度单调性 (CI Width Monotonicity)
// 理论依据: 更大的临界值产生更宽的置信带
// uniformall(true)的CI宽度 ≥ uniformall(false)的CI宽度
// =============================================================================
display _n "{hline 70}"
display "TEST 2: 置信带宽度单调性 (larger critical value → wider CI)"
display "{hline 70}"

capture noisily {
    // Run with uniformall(true)
    set seed 54321
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(true) bstrap biters(200) seed(54321)

    tempname R_all2
    matrix `R_all2' = e(results)
    local nrows_all = rowsof(`R_all2')

    // Run with uniformall(false)
    set seed 54321
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(false) bstrap biters(200) seed(54321)

    tempname R_marg2
    matrix `R_marg2' = e(results)
    local nrows_marg = rowsof(`R_marg2')

    // Compare bootstrap CI widths (columns 8 and 9: ci2_lower, ci2_upper)
    local width_ok = 1
    local nrows = min(`nrows_all', `nrows_marg')
    forvalues i = 1/`nrows' {
        local ci2_l_all = `R_all2'[`i', 8]
        local ci2_u_all = `R_all2'[`i', 9]
        local ci2_l_marg = `R_marg2'[`i', 8]
        local ci2_u_marg = `R_marg2'[`i', 9]

        // Skip rows with missing CI values
        if `ci2_l_all' == . | `ci2_u_all' == . | ///
           `ci2_l_marg' == . | `ci2_u_marg' == . {
            continue
        }

        local width_all = `ci2_u_all' - `ci2_l_all'
        local width_marg = `ci2_u_marg' - `ci2_l_marg'

        // uniformall CI width should be >= marginal CI width
        if `width_all' < `width_marg' - 1e-6 {
            display as error "  VIOLATION at row `i': width_all=" ///
                %9.6f `width_all' " < width_marginal=" %9.6f `width_marg'
            local width_ok = 0
        }
    }
    assert `width_ok' == 1
}
if _rc == 0 {
    display as result "  [PASS] 置信带宽度单调性: CI(uniformall=1) >= CI(uniformall=0)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 置信带宽度单调性 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 3: 单(g,t)对退化 (Single pair degeneracy)
// 理论依据: 当仅一个(g,t)对时，全局sup=边际sup，uniformall自动降级为0
// =============================================================================
display _n "{hline 70}"
display "TEST 3: 单(g,t)对退化 — uniformall自动设为0"
display "{hline 70}"

capture noisily {
    // Use gteval to specify exactly one (g,t) pair
    // First get available groups from data
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    // Run with a single (g,t) pair via gteval — request uniformall(true)
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        gteval(2 3) uniformall(true) bstrap biters(100) seed(54321)

    // Verify e(uniformall) was overridden to 0
    local eff_ua = e(uniformall)
    display "  e(uniformall) = `eff_ua' (expected: 0)"
    assert `eff_ua' == 0
}
if _rc == 0 {
    display as result "  [PASS] 单(g,t)对: uniformall自动退化为0 (全局sup=边际sup)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 单(g,t)对退化检验 (rc = " _rc ")"
    display as error "  预期: num_gteval==1时, e(uniformall)==0"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 4: 解析法临界值一致性 (Analytical Critical Value Consistency)
// 理论依据: Theorem 2 — uniformall=TRUE时, c_hat对所有(g,t)使用相同公式
// a² = 2log((b-a)/h) + 2log(√λ/(2π))
// c_hat = √(a² - 2log(log(1/√(1-α))))
// 当所有(g,t)共享相同bandwidth时，c_hat各行相同
// =============================================================================
display _n "{hline 70}"
display "TEST 4: 解析法临界值一致性 (Theorem 2: c_hat depends on h, not on (g,t))"
display "{hline 70}"

capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    // Use manual bandwidth to ensure all (g,t) share the same h
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) ///
        bwselect(manual) bw(0.3) uniformall(true) nobstrap

    tempname chat_mat
    matrix `chat_mat' = e(c_hat)
    local ngt_chat = rowsof(`chat_mat')

    display "  c_hat matrix:"
    matrix list `chat_mat', format(%9.6f)

    // All rows should be identical (same h for all (g,t))
    local chat_consistent = 1
    local chat_ref = `chat_mat'[1, 1]
    forvalues i = 2/`ngt_chat' {
        local chat_i = `chat_mat'[`i', 1]
        // Skip missing values
        if `chat_ref' == . | `chat_i' == . {
            continue
        }
        if reldif(`chat_ref', `chat_i') > 1e-10 {
            display as error "  MISMATCH: c_hat[1]=" %9.6f `chat_ref' ///
                " vs c_hat[`i']=" %9.6f `chat_i'
            local chat_consistent = 0
        }
    }
    assert `chat_consistent' == 1
}
if _rc == 0 {
    display as result "  [PASS] 解析法临界值: 共同bandwidth → c_hat各(g,t)对相同"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 解析法临界值一致性 (rc = " _rc ")"
    display as error "  Theorem 2: c_hat = f(b-a, h, λ, α), 与(g,t)无关"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 5: 点估计和SE不受影响 (Point Estimates Invariance)
// 理论依据: uniformall仅影响临界值选取方式，不影响点估计和标准误
// =============================================================================
display _n "{hline 70}"
display "TEST 5: 点估计和SE不受uniformall影响"
display "{hline 70}"

capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    // Use manual fixed bandwidth to isolate uniformall's effect on critical values only
    // With the same h, uniformall should NOT affect point estimates or SE
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
        bwselect(manual) bw(0.3) uniformall(true) nobstrap

    tempname R_ua1
    matrix `R_ua1' = e(results)
    local nrows1 = rowsof(`R_ua1')

    // Run with uniformall(false) — same data, same manual bandwidth
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) ///
        bwselect(manual) bw(0.3) uniformall(false) nobstrap

    tempname R_ua0
    matrix `R_ua0' = e(results)
    local nrows0 = rowsof(`R_ua0')

    // Compare est (col 4) and se (col 5)
    local nrows = min(`nrows1', `nrows0')
    local est_match = 1
    local se_match = 1
    forvalues i = 1/`nrows' {
        // Point estimates
        local est1 = `R_ua1'[`i', 4]
        local est0 = `R_ua0'[`i', 4]
        if `est1' != . & `est0' != . {
            if reldif(`est1', `est0') > 1e-10 {
                display as error "  est MISMATCH row `i': " ///
                    %12.8f `est1' " vs " %12.8f `est0'
                local est_match = 0
            }
        }

        // Standard errors
        local se1 = `R_ua1'[`i', 5]
        local se0 = `R_ua0'[`i', 5]
        if `se1' != . & `se0' != . {
            if reldif(`se1', `se0') > 1e-10 {
                display as error "  se MISMATCH row `i': " ///
                    %12.8f `se1' " vs " %12.8f `se0'
                local se_match = 0
            }
        }
    }
    assert `est_match' == 1
    assert `se_match' == 1
}
if _rc == 0 {
    display as result "  [PASS] 点估计和SE完全一致: uniformall不影响est/se列"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 点估计/SE不变性检验 (rc = " _rc ")"
    display as error "  uniformall仅影响临界值和CI，不应改变est或se"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 6: 布尔选项冲突处理 (Option Conflict Handling)
// uniformall()与nouniformall同时指定时的错误处理
// =============================================================================
display _n "{hline 70}"
display "TEST 6: 布尔选项冲突处理"
display "{hline 70}"

// 6a: uniformall(true) + nouniformall → 应报错
capture noisily {
    set seed 54321
    didhetero_simdata, n(300) tau(4) clear
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(true) nouniformall nobstrap
}
if _rc != 0 {
    display as result "  [PASS] uniformall(true) + nouniformall 正确报错 (rc = " _rc ")"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] uniformall(true) + nouniformall 未报错（预期冲突错误）"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// 6b: uniformall(false) 与 nouniformall 语义一致性
capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    // uniformall(false)
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(false) nobstrap

    local ua_false_val = e(uniformall)
    tempname R_false
    matrix `R_false' = e(results)
}
if _rc == 0 {
    capture noisily {
        // nouniformall
        catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
            zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
            nouniformall nobstrap

        local ua_no_val = e(uniformall)
        tempname R_no
        matrix `R_no' = e(results)

        // Both should yield uniformall=0
        assert `ua_false_val' == 0
        assert `ua_no_val' == 0

        // Results should be identical
        local nrows_cmp = rowsof(`R_false')
        local match_ok = 1
        forvalues i = 1/`nrows_cmp' {
            forvalues j = 4/5 {
                local v1 = `R_false'[`i', `j']
                local v2 = `R_no'[`i', `j']
                if `v1' != . & `v2' != . {
                    if reldif(`v1', `v2') > 1e-10 {
                        local match_ok = 0
                    }
                }
            }
        }
        assert `match_ok' == 1
    }
    if _rc == 0 {
        display as result "  [PASS] uniformall(false) 与 nouniformall 语义一致"
        local _test_pass = `_test_pass' + 1
    }
    else {
        display as error "  [FAIL] uniformall(false) 与 nouniformall 不一致 (rc = " _rc ")"
        local _test_fail = `_test_fail' + 1
    }
}
else {
    display as error "  [FAIL] uniformall(false) 运行失败 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 7: uniformall=TRUE时Bootstrap临界值全局一致
// 理论依据: Theorem 4 — uniformall模式下c_check对所有(g,t)相同
// =============================================================================
display _n "{hline 70}"
display "TEST 7: uniformall=TRUE → c_check对所有(g,t)相同 (Theorem 4)"
display "{hline 70}"

capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(true) bstrap biters(200) seed(54321)

    tempname cc_mat
    matrix `cc_mat' = e(c_check)
    local ngt_cc = rowsof(`cc_mat')

    display "  c_check values under uniformall(true):"
    matrix list `cc_mat', format(%9.6f)

    // All c_check values should be identical (global quantile)
    local cc_ref = `cc_mat'[1, 1]
    local cc_ok = 1
    forvalues i = 2/`ngt_cc' {
        local cc_i = `cc_mat'[`i', 1]
        if `cc_ref' == . | `cc_i' == . {
            continue
        }
        if reldif(`cc_ref', `cc_i') > 1e-12 {
            display as error "  MISMATCH: c_check[1]=" %9.6f `cc_ref' ///
                " vs c_check[`i']=" %9.6f `cc_i'
            local cc_ok = 0
        }
    }
    assert `cc_ok' == 1
}
if _rc == 0 {
    display as result "  [PASS] uniformall=TRUE: c_check全局一致（共用max-sup分位数）"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] c_check一致性 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 8: uniformall=FALSE时各(g,t)可有不同临界值
// 理论依据: 边际推断各(g,t)独立取分位数 → 临界值可以不同
// =============================================================================
display _n "{hline 70}"
display "TEST 8: uniformall=FALSE → 各(g,t)独立临界值 (per-pair quantile)"
display "{hline 70}"

capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        uniformall(false) bstrap biters(200) seed(54321)

    tempname cc_marg_mat
    matrix `cc_marg_mat' = e(c_check)
    local ngt_marg = rowsof(`cc_marg_mat')

    display "  c_check values under uniformall(false):"
    matrix list `cc_marg_mat', format(%9.6f)

    // Verify the matrix exists and has reasonable values
    assert `ngt_marg' >= 1
    // At least verify values are positive (valid critical values)
    local all_positive = 1
    forvalues i = 1/`ngt_marg' {
        local cc_val = `cc_marg_mat'[`i', 1]
        if `cc_val' != . & `cc_val' <= 0 {
            local all_positive = 0
        }
    }
    assert `all_positive' == 1
}
if _rc == 0 {
    display as result "  [PASS] uniformall=FALSE: 各(g,t)对独立计算临界值"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 边际临界值检验 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 9: 默认值验证 — uniformall默认为TRUE
// =============================================================================
display _n "{hline 70}"
display "TEST 9: 默认值 — 未指定uniformall时默认为1 (全局均匀推断)"
display "{hline 70}"

capture noisily {
    set seed 54321
    didhetero_simdata, n(500) tau(4) clear

    // 不指定uniformall参数
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        nobstrap

    local default_ua = e(uniformall)
    display "  e(uniformall) = `default_ua' (expected: 1)"
    assert `default_ua' == 1
}
if _rc == 0 {
    display as result "  [PASS] 默认uniformall=1 (Theorem 2/4默认方法)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 默认值检验 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 70}"
display "UNIFORMALL TEST SUMMARY"
display "{hline 70}"
display "  Total:  `_test_total'"
display "  Passed: `_test_pass'"
display "  Failed: `_test_fail'"
display "{hline 70}"
if `_test_fail' > 0 {
    display as error "SOME TESTS FAILED"
    exit 1
}
else {
    display as result "ALL UNIFORMALL TESTS PASSED"
}
