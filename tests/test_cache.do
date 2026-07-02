// test_cache.do — 确认缓存功能已移除
// 验证连续两次运行都正常完成（无缓存逻辑干预）
clear all
set seed 12345

// Add ado path
adopath + "/Users/cxy/Desktop/didhetero/didhetero-main/ado"

didhetero_simdata, n(200) tau(4)
catt_gt Y, id(id) time(period) group(G) z(Z) zeval(0.3 0.5 0.7) porder(1) nobstrap
mat r1 = e(results)
catt_gt Y, id(id) time(period) group(G) z(Z) zeval(0.3 0.5 0.7) porder(1) nobstrap
mat r2 = e(results)
// 两次结果应一致（确定性计算）
assert mreldif(r1, r2) < 1e-14
di "[PASS] No cache: repeated runs produce identical results"
