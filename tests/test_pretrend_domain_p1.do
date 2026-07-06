// =============================================================================
// test_pretrend_domain_p1.do
// 回归测试：Appendix F 预处理 DR 估计量的域判定 + P1 缺陷修复
//
// 理论基础（Imai, Qin, and Yanagi, 2026）：
//   - Assumption 2.5'（第2888-2894行）对预处理与后处理对平等要求对照群
//     支撑条件 s < gbar。
//   - 预处理对（t < g - delta）对照池锚定在 t = g（R_{g,g}, Eq. 2.31,
//     第2904行），基期 t = g - delta - 1 被排除（Eq. 2.32-2.33）。
//
// P1 缺陷（已修复）：pretrend 分支中后处理对曾无条件 return(1)，绕过
//   notyettreated 无 never-treated 时的 gbar 支撑检查。修复后后处理对与
//   主分支一致执行 t_ord < gbar_ord - anticipation。
//
// 测试策略：直接对 Mata 纯函数 didhetero_pair_in_domain() 进行单元测试。
//   该函数是无副作用的域判定核心（8 参数），完全绕开 catt_gt 数据管线与
//   GPS 估计，因此可对每个 (g,t) 组合在 has_never=0/1、pretrend=0/1 下
//   稳定、确定地枚举返回值，与论文 Appendix F 手工推导逐一对照。
//
//   时间支撑 t_vals = (1\2\3\4)（对应 simdata tau=4），anticipation=0。
//   NO_never: has_never=0, gbar=4（最大处理组）。
//   WITH_never: has_never=1, gbar=.（缺失，never-treated 提供对照）。
// =============================================================================
clear all
set more off

// 用 ++ 前置开发目录 ado/，确保加载本仓库最新编译的 ldidhetero.mlib，
// 覆盖 PLUS 目录 (ado/plus/l/) 中可能存在的旧安装副本（同名 mlib 首个命中优先）。
adopath ++ "/Users/cxy/Desktop/didhetero/didhetero-main/ado"
quietly mata: mata mlib index

mata:
// -------------------------------------------------------------------------
// _p1_chk(): 比较 didhetero_pair_in_domain 返回值与期望值，返回 1=通过。
// -------------------------------------------------------------------------
real scalar _p1_chk(real scalar g, real scalar t, real scalar ant,
    string scalar cg, real scalar pre, real scalar hasnv,
    real scalar gbar, real colvector tv,
    real scalar expect, string scalar label)
{
    real scalar got
    got = didhetero_pair_in_domain(g, t, ant, cg, pre, hasnv, gbar, tv)
    if (got == expect) {
        printf("  [PASS] %s: in_domain=%g (expect %g)\n", label, got, expect)
        return(1)
    }
    printf("  {err}[FAIL] %s: in_domain=%g (expect %g){txt}\n", label, got, expect)
    return(0)
}

// -------------------------------------------------------------------------
// _p1_run(): 枚举所有 (g,t) 域判定用例，返回失败计数。
// -------------------------------------------------------------------------
real scalar _p1_run()
{
    real scalar np, nf
    real colvector tv
    tv = (1 \ 2 \ 3 \ 4)
    np = 0
    nf = 0

    printf("\n{hline 70}\n")
    printf("GROUP 1: P1 核心 — NO_never + pretrend + 后处理对无有效对照 → 拒纳\n")
    printf("{hline 70}\n")
    // 后处理对 (t_ord >= g_ord): 修复后要求 t_ord < gbar_ord(=4)。
    // (4,4)/(2,4)/(3,4): t_ord=4 !< 4 → 拒纳(0)。修复前 P1 缺陷会误纳(1)。
    np = np + _p1_chk(4, 4, 0, "notyettreated", 1, 0, 4, tv, 0, "(4,4) NO_never 后处理对 t_ord=gbar_ord")
    np = np + _p1_chk(2, 4, 0, "notyettreated", 1, 0, 4, tv, 0, "(2,4) NO_never 后处理对 t_ord=gbar_ord")
    np = np + _p1_chk(3, 4, 0, "notyettreated", 1, 0, 4, tv, 0, "(3,4) NO_never 后处理对 t_ord=gbar_ord")

    printf("\n{hline 70}\n")
    printf("GROUP 2: NO_never + pretrend 有效后处理对通过域判定（未被 P1 误拒）\n")
    printf("{hline 70}\n")
    // t_ord < gbar_ord=4 → 纳入(1)。
    np = np + _p1_chk(2, 2, 0, "notyettreated", 1, 0, 4, tv, 1, "(2,2) NO_never 有效后处理对")
    np = np + _p1_chk(2, 3, 0, "notyettreated", 1, 0, 4, tv, 1, "(2,3) NO_never 有效后处理对")
    np = np + _p1_chk(3, 3, 0, "notyettreated", 1, 0, 4, tv, 1, "(3,3) NO_never 有效后处理对")

    printf("\n{hline 70}\n")
    printf("GROUP 3: 基期排除 t = g - delta - 1 在 pretrend 下保持 (Eq. 2.32-2.33)\n")
    printf("{hline 70}\n")
    // (3,2): t_ord=2 == g_ord-ant-1 = 2 → 基期，拒纳(0)，即便含 never。
    np = np + _p1_chk(3, 2, 0, "notyettreated", 1, 1, ., tv, 0, "(3,2) WITH_never 基期")
    np = np + _p1_chk(3, 2, 0, "notyettreated", 1, 0, 4, tv, 0, "(3,2) NO_never  基期")
    // (4,3): t_ord=3 == 4-0-1 = 3 → 基期，拒纳(0)。
    np = np + _p1_chk(4, 3, 0, "notyettreated", 1, 1, ., tv, 0, "(4,3) WITH_never 基期")

    printf("\n{hline 70}\n")
    printf("GROUP 4: 预处理对 (t < g - delta) 锚定 R_{g,g}，控制群支撑条件\n")
    printf("{hline 70}\n")
    // (4,2): t_ord=2 < g_ord=4 预处理对, has_never=1 → 纳入(1)（never 提供对照）。
    np = np + _p1_chk(4, 2, 0, "notyettreated", 1, 1, ., tv, 1, "(4,2) WITH_never 预处理对")
    // (4,2) NO_never: 预处理对特判 gbar_ord(=4) > g_ord+ant(=4)? 假 → 拒纳(0)。
    np = np + _p1_chk(4, 2, 0, "notyettreated", 1, 0, 4, tv, 0, "(4,2) NO_never  预处理对 gbar<=g")
    // nevertreated 控制组，含 never → 纳入(1)；无 never → 拒纳(0)。
    np = np + _p1_chk(4, 2, 0, "nevertreated", 1, 1, ., tv, 1, "(4,2) nevertreated 有 never")
    np = np + _p1_chk(4, 2, 0, "nevertreated", 1, 0, ., tv, 0, "(4,2) nevertreated 无 never")

    printf("\n{hline 70}\n")
    printf("GROUP 5: 非 pretrend 主分支不受影响（数值不变性的逻辑保证）\n")
    printf("{hline 70}\n")
    np = np + _p1_chk(2, 2, 0, "notyettreated", 0, 0, 4, tv, 1, "(2,2) 非pretrend NO_never")
    np = np + _p1_chk(4, 4, 0, "notyettreated", 0, 0, 4, tv, 0, "(4,4) 非pretrend NO_never")
    // (4,2) 非pretrend: t_ord=2 < g_ord-ant=4 → 主分支 return(0) 预处理对不纳入。
    np = np + _p1_chk(4, 2, 0, "notyettreated", 0, 1, ., tv, 0, "(4,2) 非pretrend WITH_never 预处理对")
    np = np + _p1_chk(2, 3, 0, "notyettreated", 0, 1, ., tv, 1, "(2,3) 非pretrend WITH_never")

    printf("\n{hline 70}\n")
    printf("GROUP 6: anticipation>0 下 P1 修复与基期偏移一致\n")
    printf("{hline 70}\n")
    // ant=1: (3,3) 后处理对(t_ord=3>=g_ord-ant=2), t_ord<gbar_ord-ant=3? 假 → 拒纳(0)。
    np = np + _p1_chk(3, 3, 1, "notyettreated", 1, 0, 4, tv, 0, "(3,3) ant=1 NO_never 后处理对无对照")
    // (3,2) ant=1: 后处理对(t_ord=2>=g_ord-ant=2), t_ord<gbar_ord-ant=3 → 纳入(1)。
    np = np + _p1_chk(3, 2, 1, "notyettreated", 1, 0, 4, tv, 1, "(3,2) ant=1 NO_never 有效后处理对")

    nf = 19 - np
    printf("\n{hline 70}\n")
    printf("SUMMARY: pretrend Appendix F 域判定 + P1 修复（Mata 单元测试）\n")
    printf("{hline 70}\n")
    printf("  Passed: %g / 19\n", np)
    printf("  Failed: %g\n", nf)
    printf("{hline 70}\n")
    return(nf)
}

st_numscalar("_p1_nfail", _p1_run())
end

if scalar(_p1_nfail) > 0 {
    display as error "SOME TESTS FAILED"
    exit 9
}
display as result "ALL TESTS PASSED"
