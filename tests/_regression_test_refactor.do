// =============================================================================
// End-to-End Regression Test for Refactored Code Structure
// Tests: compilation, bwselect split, utils split, error codes
// =============================================================================
clear all
set more off
set seed 12345

// Step 1: Compile
display _n "{hline 60}"
display "  STEP 1: COMPILATION (make.do)"
display "{hline 60}"

cd "/Users/cxy/Desktop/didhetero/didhetero-main"
do make.do

display _n as result ">>> COMPILATION: PASSED <<<"

// Step 2: Setup adopath
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

// Step 3: Run functional tests
display _n "{hline 60}"
display "  STEP 2: FUNCTIONAL REGRESSION TESTS"
display "{hline 60}"

// ---------- Test 1: IMSE1 Gaussian p=2 ----------
display _n ">>> TEST 1: IMSE1 Gaussian kernel p=2 <<<"
set seed 12345
didhetero_simdata, n(300) tau(4) clear
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) bwselect(IMSE1) bstrap biters(100)
display as result ">>> TEST 1 (IMSE1): PASSED <<<"

// ---------- Test 2: IMSE2 Epanechnikov p=1 ----------
display _n ">>> TEST 2: IMSE2 Epanechnikov kernel p=1 <<<"
set seed 12345
didhetero_simdata, n(300) tau(4) clear
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(epa) porder(1) bwselect(IMSE2) bstrap biters(50)
display as result ">>> TEST 2 (IMSE2): PASSED <<<"

// ---------- Test 3: US1 undersmooth ----------
display _n ">>> TEST 3: US1 undersmooth bandwidth <<<"
set seed 12345
didhetero_simdata, n(300) tau(4) clear
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(US1) bstrap biters(50)
display as result ">>> TEST 3 (US1): PASSED <<<"

// ---------- Test 4: Manual bandwidth ----------
display _n ">>> TEST 4: Manual bandwidth <<<"
set seed 12345
didhetero_simdata, n(300) tau(4) clear
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(manual) bw(0.3)
display as result ">>> TEST 4 (manual): PASSED <<<"

// ---------- Test 5: Aggregation - dynamic ----------
display _n ">>> TEST 5: Aggregation - dynamic <<<"
aggte_gt, type(dynamic) bstrap(true) biters(50)
display as result ">>> TEST 5 (dynamic): PASSED <<<"

// ---------- Test 6: Aggregation - simple ----------
display _n ">>> TEST 6: Aggregation - simple <<<"
aggte_gt, type(simple) bstrap(true) biters(50)
display as result ">>> TEST 6 (simple): PASSED <<<"

// ---------- Test 7: Aggregation - group ----------
display _n ">>> TEST 7: Aggregation - group <<<"
aggte_gt, type(group) bstrap(true) biters(50)
display as result ">>> TEST 7 (group): PASSED <<<"

// ---------- Test 8: Aggregation - calendar ----------
display _n ">>> TEST 8: Aggregation - calendar <<<"
aggte_gt, type(calendar) bstrap(true) biters(50)
display as result ">>> TEST 8 (calendar): PASSED <<<"

// ---------- Test 9: Complex formula (utils_formula) ----------
display _n ">>> TEST 9: Complex formula parsing <<<"
set seed 12345
didhetero_simdata, n(300) tau(4) clear
// Generate time-invariant covariates at the unit level
bysort id: gen X1 = rnormal() if _n == 1
bysort id: replace X1 = X1[1]
bysort id: gen X2 = rnormal() if _n == 1
bysort id: replace X2 = X2[1]
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z X1 X2) zeval(0.3 0.5 0.7) kernel(gau) porder(2) bwselect(IMSE1)
display as result ">>> TEST 9 (formula): PASSED <<<"

// ---------- Test 10: KDE bandwidth (larger sample) ----------
display _n ">>> TEST 10: KDE bandwidth (n=500) <<<"
set seed 12345
didhetero_simdata, n(500) tau(4) clear
catt_gt Y, id(id) time(period) group(G) z(Z) xformula(Z) zeval(0.2 0.4 0.6 0.8) kernel(gau) porder(2) bwselect(IMSE1) bstrap biters(50)
display as result ">>> TEST 10 (KDE): PASSED <<<"

// =============================================================================
display _n "{hline 60}"
display as result "  ALL REGRESSION TESTS PASSED SUCCESSFULLY"
display "{hline 60}"
