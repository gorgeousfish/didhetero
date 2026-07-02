*! didhetero_errors.mata
*! Unified error code and message format definitions
*! Compile order: after types.mata, before other modules

mata:

// =============================================================================
// didhetero_errors.mata
// Unified error code constants and formatted error-throwing utilities
//
// Provides:
//   - DH_ERR_ESTIMATION_FAILED() // Error code 3498
//   - DH_ERR_INVALID_INPUT()     // Error code 198
//   - DH_ERR_DIMENSION_MISMATCH()// Error code 3001
//   - DH_ERR_LAPACK_FAILURE()    // Error code 430
//   - DH_ERR_DATA_STRUCTURE()    // Error code 3200
//   - _dh_error()                // Generic formatted error
//   - _dh_error_estimation()     // Estimation failure wrapper
//   - _dh_error_input()          // Input error wrapper
//   - _dh_error_dimension()      // Dimension mismatch wrapper
//   - _dh_error_data()           // Data structure error wrapper
//
// Requires:
//   - (none -- base error utilities)
// =============================================================================


// ─── 错误码常量 ─────────────────────────────────────────────────────────────
// Stata标准错误码，通过语义别名提升可读性

real scalar DH_ERR_ESTIMATION_FAILED()
{
    return(3498)
}

real scalar DH_ERR_INVALID_INPUT()
{
    return(198)
}

real scalar DH_ERR_DIMENSION_MISMATCH()
{
    return(3001)
}

real scalar DH_ERR_LAPACK_FAILURE()
{
    return(430)
}

real scalar DH_ERR_DATA_STRUCTURE()
{
    return(3200)
}


// ─── 统一错误抛出函数 ───────────────────────────────────────────────────────
// 格式: "[模块名]: 错误描述"

void _dh_error(real scalar code, string scalar module, string scalar msg)
{
    string scalar full_msg
    full_msg = "[" + module + "] " + msg
    _error(code, full_msg)
}


// ─── 便捷包装：估计失败 (3498) ─────────────────────────────────────────────

void _dh_error_estimation(string scalar module, string scalar msg)
{
    _dh_error(DH_ERR_ESTIMATION_FAILED(), module, msg)
}


// ─── 便捷包装：用户输入错误 (198) ──────────────────────────────────────────

void _dh_error_input(string scalar module, string scalar msg)
{
    _dh_error(DH_ERR_INVALID_INPUT(), module, msg)
}


// ─── 便捷包装：维度不匹配 (3001) ──────────────────────────────────────────

void _dh_error_dimension(string scalar module, string scalar msg)
{
    _dh_error(DH_ERR_DIMENSION_MISMATCH(), module, msg)
}


// ─── 便捷包装：数据结构错误 (3200) ────────────────────────────────────────

void _dh_error_data(string scalar module, string scalar msg)
{
    _dh_error(DH_ERR_DATA_STRUCTURE(), module, msg)
}

end
