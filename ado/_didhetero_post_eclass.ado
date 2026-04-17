*! _didhetero_post_eclass.ado
*! Internal helper: turn didhetero/catt_gt output into a storable eclass object

mata:
void __dh_capture_post_state(real scalar has_c_check)
{
    external real matrix __dh_post_results
    external real matrix __dh_post_estimate_b
    external real matrix __dh_post_gteval
    external real matrix __dh_post_zeval
    external real matrix __dh_post_bw
    external real matrix __dh_post_c_hat
    external real matrix __dh_post_B_g_t
    external real matrix __dh_post_G_g
    external real matrix __dh_post_Z
    external real matrix __dh_post_dh_Y_wide
    external real matrix __dh_post_dh_G_unit
    external real matrix __dh_post_dh_t_vals
    external real matrix __dh_post_dh_gps_mat
    external real matrix __dh_post_dh_or_mat
    external real matrix __dh_post_mu_G_g
    external real matrix __dh_post_catt_est
    external real matrix __dh_post_catt_se
    external real matrix __dh_post_kd0_Z
    external real matrix __dh_post_kd1_Z
    external real matrix __dh_post_Z_supp
    external real matrix __dh_post_c_check

    __dh_post_results    = st_matrix("e(results)")
    __dh_post_estimate_b = st_matrix("e(Estimate_b)")
    __dh_post_gteval     = st_matrix("e(gteval)")
    __dh_post_zeval      = st_matrix("e(zeval)")
    __dh_post_bw         = st_matrix("e(bw)")
    __dh_post_c_hat      = st_matrix("e(c_hat)")
    __dh_post_B_g_t      = st_matrix("e(B_g_t)")
    __dh_post_G_g        = st_matrix("e(G_g)")
    __dh_post_Z          = st_matrix("e(Z)")
    __dh_post_dh_Y_wide  = st_matrix("e(dh_Y_wide)")
    __dh_post_dh_G_unit  = st_matrix("e(dh_G_unit)")
    __dh_post_dh_t_vals  = st_matrix("e(dh_t_vals)")
    __dh_post_dh_gps_mat = st_matrix("e(dh_gps_mat)")
    __dh_post_dh_or_mat  = st_matrix("e(dh_or_mat)")
    __dh_post_mu_G_g     = st_matrix("e(mu_G_g)")
    __dh_post_catt_est   = st_matrix("e(catt_est)")
    __dh_post_catt_se    = st_matrix("e(catt_se)")
    __dh_post_kd0_Z      = st_matrix("e(kd0_Z)")
    __dh_post_kd1_Z      = st_matrix("e(kd1_Z)")
    __dh_post_Z_supp     = st_matrix("e(Z_supp)")

    if (has_c_check) {
        __dh_post_c_check = st_matrix("e(c_check)")
    }
    else {
        __dh_post_c_check = J(0, 0, .)
    }
}

void __dh_restore_post_state(
    real scalar has_c_check,
    real scalar gbar,
    real scalar gbar_isinf,
    real scalar n_units,
    real scalar num_gteval,
    real scalar num_zeval,
    real scalar T_num)
{
    external real matrix __dh_post_results
    external real matrix __dh_post_estimate_b
    external real matrix __dh_post_gteval
    external real matrix __dh_post_zeval
    external real matrix __dh_post_bw
    external real matrix __dh_post_c_hat
    external real matrix __dh_post_B_g_t
    external real matrix __dh_post_G_g
    external real matrix __dh_post_Z
    external real matrix __dh_post_dh_Y_wide
    external real matrix __dh_post_dh_G_unit
    external real matrix __dh_post_dh_t_vals
    external real matrix __dh_post_dh_gps_mat
    external real matrix __dh_post_dh_or_mat
    external real matrix __dh_post_mu_G_g
    external real matrix __dh_post_catt_est
    external real matrix __dh_post_catt_se
    external real matrix __dh_post_kd0_Z
    external real matrix __dh_post_kd1_Z
    external real matrix __dh_post_Z_supp
    external real matrix __dh_post_c_check

    st_matrix("e(results)", __dh_post_results)
    st_matrix("e(Estimate)", __dh_post_results)
    st_matrix("e(Estimate_b)", __dh_post_estimate_b)
    st_matrix("e(gteval)", __dh_post_gteval)
    st_matrix("e(zeval)", __dh_post_zeval)
    st_matrix("e(bw)", __dh_post_bw)
    st_matrix("e(c_hat)", __dh_post_c_hat)
    st_matrix("e(B_g_t)", __dh_post_B_g_t)
    st_matrix("e(G_g)", __dh_post_G_g)
    st_matrix("e(Z)", __dh_post_Z)
    st_matrix("e(dh_Y_wide)", __dh_post_dh_Y_wide)
    st_matrix("e(dh_G_unit)", __dh_post_dh_G_unit)
    st_matrix("e(dh_t_vals)", __dh_post_dh_t_vals)
    st_matrix("e(dh_gps_mat)", __dh_post_dh_gps_mat)
    st_matrix("e(dh_or_mat)", __dh_post_dh_or_mat)
    st_matrix("e(mu_G_g)", __dh_post_mu_G_g)
    st_matrix("e(catt_est)", __dh_post_catt_est)
    st_matrix("e(catt_se)", __dh_post_catt_se)
    st_matrix("e(kd0_Z)", __dh_post_kd0_Z)
    st_matrix("e(kd1_Z)", __dh_post_kd1_Z)
    st_matrix("e(Z_supp)", __dh_post_Z_supp)
    if (has_c_check) {
        st_matrix("e(c_check)", __dh_post_c_check)
    }
    st_numscalar("e(gbar)", gbar)
    st_numscalar("e(gbar_isinf)", gbar_isinf)
    st_numscalar("e(N)", n_units)
    st_numscalar("e(num_gteval)", num_gteval)
    st_numscalar("e(num_zeval)", num_zeval)
    st_numscalar("e(T)", T_num)
}
end

program define _didhetero_post_eclass, eclass
    version 16.0

    capture confirm matrix e(c_check)
    local has_c_check = (_rc == 0)
    mata: __dh_capture_post_state(`has_c_check')

    local gbar = e(gbar)
    capture confirm scalar e(gbar_isinf)
    if _rc == 0 local gbar_isinf = e(gbar_isinf)
    else local gbar_isinf = missing(e(gbar))
    local n_units = e(N)
    local num_gteval = e(num_gteval)
    local num_zeval = e(num_zeval)
    local T_num = e(T)

    tempname b
    matrix `b' = e(Estimate_b)

    // The package's authoritative inference lives in e(results)/e(catt_se).
    // Post only e(b) so estimates store/restore remains available without
    // advertising an unsupported covariance matrix for generic Wald commands.
    local stripe_names
    forvalues i = 1/`=colsof(`b')' {
        local b_i = el(`b', 1, `i')
        if missing(`b_i') {
            matrix `b'[1, `i'] = 0
        }
        local stripe_names `stripe_names' catt_`i'
    }

    matrix colnames `b' = `stripe_names'

    // Build a conformable zero V so that e(V) exists and estimates store works.
    // The package's authoritative inference lives in e(results)/e(catt_se);
    // the zero V signals that generic Wald tests on e(b) are not supported.
    tempname V
    local _k = colsof(`b')
    matrix `V' = J(`_k', `_k', 0)
    matrix colnames `V' = `stripe_names'
    matrix rownames `V' = `stripe_names'

    ereturn post `b' `V', obs(`n_units')
    mata: __dh_restore_post_state(`has_c_check', `gbar', `gbar_isinf', `n_units', `num_gteval', `num_zeval', `T_num')
    mata: st_matrixcolstripe("e(results)", ///
        (J(10, 1, ""), ("g" \ "t" \ "z" \ "est" \ "se" \ ///
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))
    mata: st_matrixcolstripe("e(Estimate)", ///
        (J(10, 1, ""), ("g" \ "t" \ "z" \ "est" \ "se" \ ///
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))
    capture confirm matrix e(gteval)
    if !_rc {
        mata: st_matrixcolstripe("e(gteval)", (J(2, 1, ""), ("g" \ "t")))
    }
end
