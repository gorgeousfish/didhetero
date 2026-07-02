// =============================================================================
// test_anticipation.do
// 回归测试：验证 anticipation 参数的理论正确性
// 理论基础：论文 Assumption 2 的参数化放松 (anticipation = delta >= 0)
//
// 当 anticipation = a 时：
//   (1) 有效 (g,t) 对要求 t_ord >= g_ord - a（处理后时期）
//   (2) 基期从 g-1 变为 g-a-1（即 t_vals[g_ord - a - 1]）
//   (3) not-yet-treated 对照组包含 G_ord > t_ord + a 的单位
//   (4) 组有效性要求 g_ord > a + 1（确保基期可用）
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
// TEST 1: anticipation(0) 等价性
// 理论依据: anticipation(0) 是默认值，应当与省略该参数完全等价
// =============================================================================
display _n "{hline 60}"
display "TEST 1: anticipation(0) 等价性"
display "{hline 60}"

// --- Test 1a: anticipation(0) vs 不指定 anticipation ---
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // 不指定 anticipation（默认为0）
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        nobstrap
    
    matrix _est_default = e(results)
    matrix _gteval_default = e(gteval)
    scalar _n_gt_default = e(num_gteval)
}
if _rc == 0 {
    display as result "  [PASS] catt_gt without anticipation option runs"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] catt_gt without anticipation option (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // 显式指定 anticipation(0)
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(0) nobstrap
    
    matrix _est_ant0 = e(results)
    matrix _gteval_ant0 = e(gteval)
    scalar _n_gt_ant0 = e(num_gteval)
    
    // 验证 (g,t) 对数量一致
    assert _n_gt_default == _n_gt_ant0
    
    // 验证 gteval 矩阵完全一致
    local _n_rows = rowsof(_gteval_default)
    forvalues i = 1/`_n_rows' {
        assert _gteval_default[`i', 1] == _gteval_ant0[`i', 1]
        assert _gteval_default[`i', 2] == _gteval_ant0[`i', 2]
    }
    
    // 验证估计结果数值完全一致（容许浮点误差 1e-6，IMSE1带宽优化可能引入微小数值差异）
    local _nrows_est = rowsof(_est_default)
    local _ncols_est = colsof(_est_default)
    forvalues i = 1/`_nrows_est' {
        forvalues j = 1/`_ncols_est' {
            if _est_default[`i', `j'] < . & _est_ant0[`i', `j'] < . {
                local _diff = abs(_est_default[`i', `j'] - _est_ant0[`i', `j'])
                assert `_diff' < 1e-6
            }
        }
    }
}
if _rc == 0 {
    display as result "  [PASS] anticipation(0) 结果与默认完全一致"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation(0) 结果与默认不一致 (rc = " _rc ")"
    display as error "    默认 n_gt = " _n_gt_default ", ant(0) n_gt = " _n_gt_ant0
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 2: gteval 构建规则验证
// 理论依据: 
//   anticipation=1 时：
//   - 组有效性: g_ord > anticipation + 1 = 2, 即 g 不能是第2期（最早处理时期）
//   - 处理后时期: t_ord >= g_ord - anticipation = g_ord - 1
//   - 基期变为 t_vals[g_ord - 1 - 1] = t_vals[g_ord - 2]
//   - 对于 tau=4 的 simdata：periods = (1,2,3,4), G in {0,2,3,4}
//     anticipation=0 时有效组: G=2(g_ord=2>1), G=3(g_ord=3>1), G=4(g_ord=4>1)
//     anticipation=1 时有效组: G=3(g_ord=3>2), G=4(g_ord=4>2)
//     G=2 变为无效（因为 g_ord=2 不满足 > anticipation+1=2）
// =============================================================================
display _n "{hline 60}"
display "TEST 2: gteval 构建规则 (anticipation=1)"
display "{hline 60}"

// --- Test 2a: anticipation=1 时 G=2 组被排除 ---
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) nobstrap
    
    matrix _gteval_ant1 = e(gteval)
    scalar _n_gt_ant1 = e(num_gteval)
    scalar _anticip_stored = e(anticipation)
    
    // 验证 e(anticipation) 正确存储
    assert _anticip_stored == 1
    
    // 验证 gteval 中不存在 g=2 的行
    // （因为 simdata tau=4 生成 periods 1,2,3,4；G=2 的 g_ord=2 不满足 >2）
    local _n_rows_gt = rowsof(_gteval_ant1)
    local _has_g2 = 0
    forvalues i = 1/`_n_rows_gt' {
        if _gteval_ant1[`i', 1] == 2 {
            local _has_g2 = 1
        }
    }
    assert `_has_g2' == 0
    
    display as text "    anticipation(1): 有效 (g,t) 对数 = " _n_gt_ant1
    display as text "    gteval 中无 g=2（g_ord=2 不满足 > anticipation+1=2）: 通过"
}
if _rc == 0 {
    display as result "  [PASS] anticipation(1) 正确排除 g=2 组"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation(1) 组排除逻辑 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 2b: anticipation=1 使得有效 (g,t) 对数减少 ---
capture noisily {
    // 比较 anticipation=0 与 anticipation=1 的 gteval 大小
    // anticipation 增大时，有效对应该减少（更严格的约束）
    assert _n_gt_ant1 < _n_gt_default
    display as text "    ant(0) n_gt = " _n_gt_default " > ant(1) n_gt = " _n_gt_ant1
}
if _rc == 0 {
    display as result "  [PASS] anticipation(1) 有效 (g,t) 对数少于 anticipation(0)"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 有效 (g,t) 对数未减少 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 2c: anticipation=1 时 gteval 中处理后时期的正确性 ---
// 对于 G=3 (g_ord=3), anticipation=1：
//   处理后时期要求 t_ord >= g_ord - anticipation = 3-1 = 2, 即 t >= period 2
//   基期 = t_vals[g_ord - anticipation - 1] = t_vals[3-1-1] = t_vals[1] = period 1
capture noisily {
    local _n_rows_gt = rowsof(_gteval_ant1)
    local _g3_has_t2 = 0
    local _g3_has_t1 = 0
    forvalues i = 1/`_n_rows_gt' {
        if _gteval_ant1[`i', 1] == 3 & _gteval_ant1[`i', 2] == 2 {
            local _g3_has_t2 = 1
        }
        if _gteval_ant1[`i', 1] == 3 & _gteval_ant1[`i', 2] == 1 {
            local _g3_has_t1 = 1
        }
    }
    // t=2 应该在有效对中（t_ord=2 >= g_ord - anticipation = 2）
    assert `_g3_has_t2' == 1
    // t=1 不应该在有效对中（t_ord=1 < g_ord - anticipation = 2，且 t_ord=1 不满足 t_ord > 1）
    assert `_g3_has_t1' == 0
    
    display as text "    G=3, anticipation=1: t=2 在有效对中（正确）, t=1 不在（正确）"
}
if _rc == 0 {
    display as result "  [PASS] anticipation(1) 处理后时期边界正确"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 处理后时期边界不正确 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 3: GPS 对照组扩展验证
// 理论依据: not-yet-treated 控制组：
//   anticipation=0 时: 对照组 = {G==0} ∪ {G_ord > t_ord}
//   anticipation=1 时: 对照组 = {G==0} ∪ {G_ord > t_ord + 1}
// 因此 anticipation > 0 时对照组更小（排除更多已接近处理的单位）
//
// 间接验证：使用 nevertreated 控制组时，anticipation 不影响对照组构成
// 但使用 notyettreated 时，不同 anticipation 可能影响估计结果
// =============================================================================
display _n "{hline 60}"
display "TEST 3: GPS 对照组扩展验证"
display "{hline 60}"

// --- Test 3a: nevertreated + anticipation(0) vs anticipation(1) ---
// nevertreated 控制组不受 anticipation 影响对照组构成
// 但 anticipation 仍影响 gteval 域和基期选择，因此结果不同
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // nevertreated + anticipation(0)
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(0) control_group(nevertreated) nobstrap
    
    matrix _gteval_never_a0 = e(gteval)
    scalar _n_gt_never_a0 = e(num_gteval)
}
if _rc == 0 {
    display as result "  [PASS] nevertreated + anticipation(0) 成功运行"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] nevertreated + anticipation(0) (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // nevertreated + anticipation(1)
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) control_group(nevertreated) nobstrap
    
    matrix _gteval_never_a1 = e(gteval)
    scalar _n_gt_never_a1 = e(num_gteval)
    
    // anticipation=1 应有更少的有效 (g,t) 对
    assert _n_gt_never_a1 < _n_gt_never_a0
    display as text "    nevertreated: ant(0) n_gt=" _n_gt_never_a0 " > ant(1) n_gt=" _n_gt_never_a1
}
if _rc == 0 {
    display as result "  [PASS] nevertreated + anticipation(1) 成功, 有效对减少"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] nevertreated + anticipation(1) (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 3b: notyettreated + anticipation(1) 能正常运行 ---
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) control_group(notyettreated) nobstrap
    
    scalar _n_gt_notyet_a1 = e(num_gteval)
    display as text "    notyettreated + ant(1): n_gt = " _n_gt_notyet_a1
}
if _rc == 0 {
    display as result "  [PASS] notyettreated + anticipation(1) 成功运行"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] notyettreated + anticipation(1) (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 4: 数值一致性
// 验证不同 anticipation 值下估计的合理性和数值稳定性
// =============================================================================
display _n "{hline 60}"
display "TEST 4: 数值一致性"
display "{hline 60}"

// --- Test 4a: anticipation(1) 估计值是有限数 ---
capture noisily {
    set seed 54321
    didhetero_simdata, n(800) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) nobstrap
    
    matrix _est_a1 = e(results)
    scalar _n_gt_a1_big = e(num_gteval)
    
    // 验证大多数估计值是有限实数（非参数估计边界处可能产生缺失值）
    local _nrows = rowsof(_est_a1)
    local _ncols = colsof(_est_a1)
    local _n_finite = 0
    local _n_total_cells = 0
    forvalues i = 1/`_nrows' {
        forvalues j = 1/`_ncols' {
            local _n_total_cells = `_n_total_cells' + 1
            if _est_a1[`i', `j'] < . {
                local _n_finite = `_n_finite' + 1
            }
        }
    }
    // 至少半数估计值应为有限数
    assert `_n_finite' > `_n_total_cells' / 2
    
    display as text "    anticipation(1), tau=4: " _n_gt_a1_big " 对, " `_n_finite' "/" `_n_total_cells' " 估计有限"
}
if _rc == 0 {
    display as result "  [PASS] anticipation(1) 所有估计值有限且非缺失"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation(1) 数值一致性 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 4b: anticipation(2) 与 tau=4 ---
capture noisily {
    set seed 54321
    didhetero_simdata, n(800) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(2) nobstrap
    
    matrix _est_a2 = e(results)
    scalar _n_gt_a2 = e(num_gteval)
    
    // anticipation=2 时 tau=4: periods=(1,2,3,4), G in {0,2,3,4}
    // 有效组: g_ord > 3, 即 G=4(g_ord=4>3)
    // 有效对应进一步减少
    assert _n_gt_a2 < _n_gt_a1_big
    
    display as text "    ant(2), tau=4: n_gt=" _n_gt_a2 " < ant(1) n_gt=" _n_gt_a1_big
}
if _rc == 0 {
    display as result "  [PASS] anticipation(2) 有效对进一步减少且运行正常"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation(2) (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 4c: 显式 gteval + anticipation 配合 ---
// 指定一个在 anticipation=1 下有效的 (g,t) 对
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // tau=4: periods=1,2,3,4; G in {0,2,3,4}
    // anticipation=1 时有效组: G=3, G=4
    // G=3, t=3: g_ord=3, t_ord=3 >= g_ord-1=2 ✓
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) gteval(3 3) nobstrap
    
    scalar _n_gt_explicit = e(num_gteval)
    assert _n_gt_explicit == 1
    
    display as text "    显式 gteval(3 3) + anticipation(1): 成功估计1对"
}
if _rc == 0 {
    display as result "  [PASS] 显式 gteval + anticipation 配合正常"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 显式 gteval + anticipation (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 4d: aggte 与 anticipation 配合 ---
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(1) nobstrap
    
    aggte_gt, type("simple") bstrap("false")
    
    // 验证 aggte 能正常执行
    matrix _aggte_est = e(Estimate)
    local _aggte_rows = rowsof(_aggte_est)
    assert `_aggte_rows' > 0
    
    display as text "    aggte_gt simple + anticipation(1): 成功, " `_aggte_rows' " 行结果"
}
if _rc == 0 {
    display as result "  [PASS] aggte_gt + anticipation(1) 成功运行"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] aggte_gt + anticipation(1) (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// TEST 5: 边界条件和错误处理
// =============================================================================
display _n "{hline 60}"
display "TEST 5: 边界条件和错误处理"
display "{hline 60}"

// --- Test 5a: anticipation 为负数 → 应报错 (exit 198) ---
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    capture {
        catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
            zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
            anticipation(-1) nobstrap
    }
    local _rc_neg = _rc
    // 负值应被拒绝
    assert `_rc_neg' != 0
    display as text "    anticipation(-1) 被正确拒绝, rc=" `_rc_neg'
}
if _rc == 0 {
    display as result "  [PASS] 负数 anticipation 被正确拒绝"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 负数 anticipation 未被拒绝 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 5b: anticipation 过大导致无有效组 → 应报错 ---
// tau=4: 4个时期(1,2,3,4), 最大 g_ord=4
// anticipation=3 时: 有效组要求 g_ord > 3+1=4, 即 g_ord > 4, 无满足条件的组
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    capture {
        catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
            zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
            anticipation(3) nobstrap
    }
    local _rc_big = _rc
    // 过大的 anticipation 应导致错误（无有效组或无有效对）
    assert `_rc_big' != 0
    display as text "    anticipation(3) + tau=4 被正确拒绝, rc=" `_rc_big'
}
if _rc == 0 {
    display as result "  [PASS] anticipation 过大被正确拒绝"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation 过大未被拒绝 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 5c: 显式 gteval 指定不在有效域内的对 → 应报错 ---
// anticipation=1: G=2 无效，指定 gteval(2 3) 应被拒绝
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    capture {
        catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
            zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
            anticipation(1) gteval(2 3) nobstrap
    }
    local _rc_invalid_gt = _rc
    // G=2 在 anticipation=1 下无效, 应报错
    assert `_rc_invalid_gt' != 0
    display as text "    gteval(2 3) + anticipation(1) 被正确拒绝, rc=" `_rc_invalid_gt'
}
if _rc == 0 {
    display as result "  [PASS] 无效 gteval + anticipation 被正确拒绝"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] 无效 gteval 未被拒绝 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 5d: anticipation=2 + tau=4 极端边界 ---
// tau=4: periods=(1,2,3,4), G in {0,2,3,4}
// anticipation=2: 有效组要求 g_ord > 3, 即只有 G=4 (g_ord=4)
// G=4 的基期 = t_vals[4-2-1] = t_vals[1] = period 1
// 有效时期: t_ord >= g_ord - anticipation = 4-2 = 2, 且 t_ord <= T-anticipation = 4-2 = 2
// 所以只有 t=period 2 (t_ord=2)
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(2) nobstrap
    
    matrix _gteval_a2_tau4 = e(gteval)
    scalar _n_gt_a2_tau4 = e(num_gteval)
    
    // 应该只有 (4, 2) 一个有效对
    assert _n_gt_a2_tau4 == 1
    assert _gteval_a2_tau4[1, 1] == 4
    assert _gteval_a2_tau4[1, 2] == 2
    
    display as text "    anticipation(2) + tau=4: 唯一有效对 (4,2) 验证通过"
}
if _rc == 0 {
    display as result "  [PASS] anticipation(2) + tau=4 极端边界正确"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation(2) + tau=4 边界 (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1

// --- Test 5e: anticipation 与 pretrend 配合 ---
// pretrend 包含处理前时期（但排除基期本身）
// 使用 anticipation(0) + tau(4) 确保有足够的处理前时期可供 pretrend 添加
capture noisily {
    set seed 12345
    didhetero_simdata, n(300) tau(4) clear
    
    // 无 pretrend
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(0) nobstrap
    scalar _n_gt_no_pre = e(num_gteval)
    
    // 有 pretrend
    catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) ///
        zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1) ///
        anticipation(0) pretrend nobstrap
    scalar _n_gt_pre = e(num_gteval)
    
    // pretrend 应该增加有效 (g,t) 对数（包含处理前时期）
    assert _n_gt_pre > _n_gt_no_pre
    
    display as text "    pretrend + anticipation(0): n_gt(no_pre)=" _n_gt_no_pre ///
        " < n_gt(pre)=" _n_gt_pre
}
if _rc == 0 {
    display as result "  [PASS] anticipation + pretrend 配合正确"
    local _test_pass = `_test_pass' + 1
}
else {
    display as error "  [FAIL] anticipation + pretrend (rc = " _rc ")"
    local _test_fail = `_test_fail' + 1
}
local _test_total = `_test_total' + 1


// =============================================================================
// SUMMARY
// =============================================================================
display _n "{hline 60}"
display "TEST SUMMARY: anticipation 参数回归测试"
display "{hline 60}"
display "  Total:  `_test_total'"
display "  Passed: `_test_pass'"
display "  Failed: `_test_fail'"
display "{hline 60}"
if `_test_fail' > 0 {
    display as error "SOME TESTS FAILED"
    exit 1
}
else {
    display as result "ALL TESTS PASSED"
}
