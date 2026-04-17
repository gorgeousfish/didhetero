*! _didhetero_post_aggte_eclass.ado
*! Internal helper: turn aggte_gt output into a coherent storable eclass object

mata:
void __dh_capture_aggte_post_state(real scalar has_c_check)
{
    external real matrix __dh_agg_results
    external real matrix __dh_agg_estimate
    external real matrix __dh_agg_aggte_est
    external real matrix __dh_agg_aggte_se
    external real matrix __dh_agg_aggte_ci1_lower
    external real matrix __dh_agg_aggte_ci1_upper
    external real matrix __dh_agg_aggte_ci2_lower
    external real matrix __dh_agg_aggte_ci2_upper
    external real matrix __dh_agg_aggte_bw
    external real matrix __dh_agg_aggte_eval
    external real matrix __dh_agg_aggte_zeval
    external real matrix __dh_agg_gteval
    external real matrix __dh_agg_zeval
    external real matrix __dh_agg_bw
    external real matrix __dh_agg_c_hat
    external real matrix __dh_agg_B_g_t
    external real matrix __dh_agg_G_g
    external real matrix __dh_agg_Z
    external real matrix __dh_agg_dh_Y_wide
    external real matrix __dh_agg_dh_G_unit
    external real matrix __dh_agg_dh_t_vals
    external real matrix __dh_agg_dh_gps_mat
    external real matrix __dh_agg_dh_or_mat
    external real matrix __dh_agg_mu_G_g
    external real matrix __dh_agg_catt_est
    external real matrix __dh_agg_catt_se
    external real matrix __dh_agg_kd0_Z
    external real matrix __dh_agg_kd1_Z
    external real matrix __dh_agg_Z_supp
    external real matrix __dh_agg_c_check

    __dh_agg_results         = st_matrix("e(results)")
    __dh_agg_estimate        = st_matrix("e(Estimate)")
    __dh_agg_aggte_est       = st_matrix("e(aggte_est)")
    __dh_agg_aggte_se        = st_matrix("e(aggte_se)")
    __dh_agg_aggte_ci1_lower = st_matrix("e(aggte_ci1_lower)")
    __dh_agg_aggte_ci1_upper = st_matrix("e(aggte_ci1_upper)")
    __dh_agg_aggte_ci2_lower = st_matrix("e(aggte_ci2_lower)")
    __dh_agg_aggte_ci2_upper = st_matrix("e(aggte_ci2_upper)")
    __dh_agg_aggte_bw        = st_matrix("e(aggte_bw)")
    __dh_agg_aggte_eval      = st_matrix("e(aggte_eval)")
    __dh_agg_aggte_zeval     = st_matrix("e(aggte_zeval)")
    __dh_agg_gteval          = st_matrix("e(gteval)")
    __dh_agg_zeval           = st_matrix("e(zeval)")
    __dh_agg_bw              = st_matrix("e(bw)")
    __dh_agg_c_hat           = st_matrix("e(c_hat)")
    __dh_agg_B_g_t           = st_matrix("e(B_g_t)")
    __dh_agg_G_g             = st_matrix("e(G_g)")
    __dh_agg_Z               = st_matrix("e(Z)")
    __dh_agg_dh_Y_wide       = st_matrix("e(dh_Y_wide)")
    __dh_agg_dh_G_unit       = st_matrix("e(dh_G_unit)")
    __dh_agg_dh_t_vals       = st_matrix("e(dh_t_vals)")
    __dh_agg_dh_gps_mat      = st_matrix("e(dh_gps_mat)")
    __dh_agg_dh_or_mat       = st_matrix("e(dh_or_mat)")
    __dh_agg_mu_G_g          = st_matrix("e(mu_G_g)")
    __dh_agg_catt_est        = st_matrix("e(catt_est)")
    __dh_agg_catt_se         = st_matrix("e(catt_se)")
    __dh_agg_kd0_Z           = st_matrix("e(kd0_Z)")
    __dh_agg_kd1_Z           = st_matrix("e(kd1_Z)")
    __dh_agg_Z_supp          = st_matrix("e(Z_supp)")

    if (has_c_check) {
        __dh_agg_c_check = st_matrix("e(c_check)")
    }
    else {
        __dh_agg_c_check = J(0, 0, .)
    }
}

void __dh_restore_aggte_post_state(
    real scalar has_c_check,
    real scalar gbar,
    real scalar gbar_isinf,
    real scalar n_units,
    real scalar num_gteval,
    real scalar num_zeval,
    real scalar T_num,
    string scalar b_matname)
{
    external real matrix __dh_agg_results
    external real matrix __dh_agg_estimate
    external real matrix __dh_agg_aggte_est
    external real matrix __dh_agg_aggte_se
    external real matrix __dh_agg_aggte_ci1_lower
    external real matrix __dh_agg_aggte_ci1_upper
    external real matrix __dh_agg_aggte_ci2_lower
    external real matrix __dh_agg_aggte_ci2_upper
    external real matrix __dh_agg_aggte_bw
    external real matrix __dh_agg_aggte_eval
    external real matrix __dh_agg_aggte_zeval
    external real matrix __dh_agg_gteval
    external real matrix __dh_agg_zeval
    external real matrix __dh_agg_bw
    external real matrix __dh_agg_c_hat
    external real matrix __dh_agg_B_g_t
    external real matrix __dh_agg_G_g
    external real matrix __dh_agg_Z
    external real matrix __dh_agg_dh_Y_wide
    external real matrix __dh_agg_dh_G_unit
    external real matrix __dh_agg_dh_t_vals
    external real matrix __dh_agg_dh_gps_mat
    external real matrix __dh_agg_dh_or_mat
    external real matrix __dh_agg_mu_G_g
    external real matrix __dh_agg_catt_est
    external real matrix __dh_agg_catt_se
    external real matrix __dh_agg_kd0_Z
    external real matrix __dh_agg_kd1_Z
    external real matrix __dh_agg_Z_supp
    external real matrix __dh_agg_c_check

    st_matrix("e(results)", __dh_agg_results)
    st_matrix("e(Estimate)", __dh_agg_estimate)
    st_matrix("e(Estimate_b)", st_matrix(b_matname))
    st_matrix("e(aggte_est)", __dh_agg_aggte_est)
    st_matrix("e(aggte_se)", __dh_agg_aggte_se)
    st_matrix("e(aggte_ci1_lower)", __dh_agg_aggte_ci1_lower)
    st_matrix("e(aggte_ci1_upper)", __dh_agg_aggte_ci1_upper)
    st_matrix("e(aggte_ci2_lower)", __dh_agg_aggte_ci2_lower)
    st_matrix("e(aggte_ci2_upper)", __dh_agg_aggte_ci2_upper)
    st_matrix("e(aggte_bw)", __dh_agg_aggte_bw)
    st_matrix("e(aggte_eval)", __dh_agg_aggte_eval)
    st_matrix("e(aggte_zeval)", __dh_agg_aggte_zeval)
    st_matrix("e(gteval)", __dh_agg_gteval)
    st_matrix("e(zeval)", __dh_agg_zeval)
    st_matrix("e(bw)", __dh_agg_bw)
    st_matrix("e(c_hat)", __dh_agg_c_hat)
    st_matrix("e(B_g_t)", __dh_agg_B_g_t)
    st_matrix("e(G_g)", __dh_agg_G_g)
    st_matrix("e(Z)", __dh_agg_Z)
    st_matrix("e(dh_Y_wide)", __dh_agg_dh_Y_wide)
    st_matrix("e(dh_G_unit)", __dh_agg_dh_G_unit)
    st_matrix("e(dh_t_vals)", __dh_agg_dh_t_vals)
    st_matrix("e(dh_gps_mat)", __dh_agg_dh_gps_mat)
    st_matrix("e(dh_or_mat)", __dh_agg_dh_or_mat)
    st_matrix("e(mu_G_g)", __dh_agg_mu_G_g)
    st_matrix("e(catt_est)", __dh_agg_catt_est)
    st_matrix("e(catt_se)", __dh_agg_catt_se)
    st_matrix("e(kd0_Z)", __dh_agg_kd0_Z)
    st_matrix("e(kd1_Z)", __dh_agg_kd1_Z)
    st_matrix("e(Z_supp)", __dh_agg_Z_supp)
    if (has_c_check) {
        st_matrix("e(c_check)", __dh_agg_c_check)
    }
    st_numscalar("e(gbar)", gbar)
    st_numscalar("e(gbar_isinf)", gbar_isinf)
    st_numscalar("e(N)", n_units)
    st_numscalar("e(num_gteval)", num_gteval)
    st_numscalar("e(num_zeval)", num_zeval)
    st_numscalar("e(T)", T_num)
}
end

program define _didhetero_post_aggte_eclass, eclass
    version 16.0

    capture confirm matrix e(c_check)
    local has_c_check = (_rc == 0)
    mata: __dh_capture_aggte_post_state(`has_c_check')

    local gbar = e(gbar)
    capture confirm scalar e(gbar_isinf)
    if _rc == 0 local gbar_isinf = e(gbar_isinf)
    else local gbar_isinf = missing(e(gbar))
    local n_units = e(N)
    local num_gteval = e(num_gteval)
    local num_zeval = e(num_zeval)
    local T_num = e(T)

    tempname b V agg_estimate
    matrix `agg_estimate' = e(Estimate)
    local n_agg = rowsof(`agg_estimate')
    matrix `b' = J(1, `n_agg', 0)
    matrix `V' = J(`n_agg', `n_agg', 0)

    local stripe_names
    forvalues i = 1/`n_agg' {
        local b_i = el(`agg_estimate', `i', 3)
        local se_i = el(`agg_estimate', `i', 4)
        if !missing(`b_i') {
            matrix `b'[1, `i'] = `b_i'
        }
        if !missing(`se_i') {
            matrix `V'[`i', `i'] = (`se_i')^2
        }
        local stripe_names `stripe_names' aggte_`i'
    }

    matrix colnames `b' = `stripe_names'
    matrix colnames `V' = `stripe_names'
    matrix rownames `V' = `stripe_names'

    ereturn post `b' `V', obs(`n_units')
    mata: __dh_restore_aggte_post_state(`has_c_check', `gbar', `gbar_isinf', `n_units', `num_gteval', `num_zeval', `T_num', "`b'")
    mata: st_matrixcolstripe("e(results)", ///
        (J(10, 1, ""), ("g" \ "t" \ "z" \ "est" \ "se" \ ///
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))
    mata: st_matrixcolstripe("e(Estimate)", ///
        (J(9, 1, ""), ("eval" \ "z" \ "est" \ "se" \ ///
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))
    mata: st_matrixcolstripe("e(aggte_eval)", (J(1, 1, ""), "eval"))
    capture confirm matrix e(gteval)
    if !_rc {
        mata: st_matrixcolstripe("e(gteval)", (J(2, 1, ""), ("g" \ "t")))
    }
end
