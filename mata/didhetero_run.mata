mata:

// =============================================================================
// didhetero_run.mata
// Full estimation pipeline orchestration for didhetero command
//
// Functions:
//   1. didhetero_bw_preloop()      - BW selection pre-loop
//   2. didhetero_run_from_ado()    - Full pipeline called from didhetero.ado
//   3. didhetero_post_results()    - Post results to Stata e() and r()
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Section 4 (CATT estimation pipeline)
// =============================================================================


// -----------------------------------------------------------------------------
// didhetero_bw_preloop()
// Bandwidth selection pre-loop: for each (g,t) pair, construct intermediate
// variables using p=1 LLR, build influence function B, then call the
// bandwidth selector.
//
// For each (g,t):
//     1. Construct G_ig, R_g, Y_diff via didhetero_intermediate_vars()
//     2. LLR (p=1) estimation of mu_G_g, mu_R_g, mu_E_g_t, mu_F_g_t
//     3. Construct A_{i,g,t}(z_r) and B_{i,g,t}(z_r)
//     4. Store B into data.B_g_t for bandwidth selection
//   Then call _didhetero_bwselect_all() on the B matrices.
//
// Args:
//   data          - DidHeteroData struct (modified in place)
//   gps_mat       - GPS result matrix from Stage 1
//   or_mat        - OR result matrix from Stage 1
//   kd0_Z         - num_zeval x 1 kernel density estimates
//   kd1_Z         - num_zeval x 1 kernel density derivative estimates
//   bwselect      - bandwidth selection method string
//   bw_manual     - manual bandwidth (scalar or vector), used when bwselect="manual"
//
// Returns:
//   K x 1 vector of bandwidths
//
// Paper ref: Section 4.2 (bandwidth selection)
// -----------------------------------------------------------------------------
real colvector didhetero_bw_preloop(
    struct DidHeteroData scalar data,
    real matrix gps_mat,
    real matrix or_mat,
    real colvector kd0_Z,
    real colvector kd1_Z,
    string scalar bwselect,
    real colvector bw_manual)
{
    real scalar K, n, num_zeval, id_gt, g1, t1, r
    real colvector bw_vec
    real colvector mu_G_col, mu_R_col, mu_G_bw, mu_R_bw
    real colvector E_g_t_vec, F_g_t_vec, mu_E_bw, mu_F_bw
    real matrix A_mat, B_mat
    struct DidHeteroIntermediate scalar intermed

    K = data.num_gteval
    n = data.n
    num_zeval = data.num_zeval

    // === If manual, skip the pre-loop entirely ===
    if (bwselect == "manual") {
        if (rows(bw_manual) == 1) {
            bw_vec = J(K, 1, bw_manual[1])
        }
        else if (rows(bw_manual) == K) {
            bw_vec = bw_manual
        }
        else {
            _error("bw must be scalar or vector of length " + strofreal(K))
        }
        // uniformall or scalar manual: use min
        if (data.uniformall | rows(bw_manual) == 1) {
            bw_vec = J(K, 1, min(bw_vec))
        }
        return(bw_vec)
    }

    // === Initialize arrays for BW selection ===
    // These will be overwritten by stage23 later
    data.A_g_t = J(1, K, NULL)
    data.B_g_t = J(1, K, NULL)
    if (rows(data.G_g) != n | cols(data.G_g) != K) {
        data.G_g = J(n, K, 0)
    }
    if (rows(data.mu_E_g_t) != num_zeval | cols(data.mu_E_g_t) != K) {
        data.mu_E_g_t = J(num_zeval, K, .)
    }
    if (rows(data.mu_F_g_t) != num_zeval | cols(data.mu_F_g_t) != K) {
        data.mu_F_g_t = J(num_zeval, K, .)
    }
    data.mu_G_g = J(num_zeval, K, .)

    // === Pre-loop: construct intermediate vars and B for each (g,t) ===
    for (id_gt = 1; id_gt <= K; id_gt++) {
        g1 = data.gteval[id_gt, 1]
        t1 = data.gteval[id_gt, 2]

        // Step 1: Intermediate variables
        intermed = didhetero_intermediate_vars(data, gps_mat, or_mat,
                       g1, t1, id_gt, data.control_group,
                       data.anticipation, data.G_g)

        // Step 2: LLR (p=1) estimation of mu_G, mu_R, mu_E, mu_F

        // mu_G_g bandwidth and estimate (p=1)
        mu_G_bw = _didhetero_lpbwselect_mse(data.G_g[., id_gt], data.Z,
                      data.zeval, 1, 0, data.kernel)
        data.mu_G_g[., id_gt] = didhetero_lpr(data.G_g[., id_gt], data.Z,
                      data.zeval, 1, 0, data.kernel, mu_G_bw)

        // mu_R bandwidth and estimate (p=1)
        mu_R_bw = _didhetero_lpbwselect_mse(intermed.R_g, data.Z,
                      data.zeval, 1, 0, data.kernel)
        mu_R_col = didhetero_lpr(intermed.R_g, data.Z,
                      data.zeval, 1, 0, data.kernel, mu_R_bw)

        // E_{g,t} = R_g * Y_diff
        E_g_t_vec = intermed.R_g :* intermed.Y_diff
        mu_E_bw = _didhetero_lpbwselect_mse(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, data.kernel)
        data.mu_E_g_t[., id_gt] = didhetero_lpr(E_g_t_vec, data.Z,
                      data.zeval, 1, 0, data.kernel, mu_E_bw)

        // F_{g,t} = G_ig * Y_diff
        F_g_t_vec = intermed.G_ig :* intermed.Y_diff
        mu_F_bw = _didhetero_lpbwselect_mse(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, data.kernel)
        data.mu_F_g_t[., id_gt] = didhetero_lpr(F_g_t_vec, data.Z,
                      data.zeval, 1, 0, data.kernel, mu_F_bw)

        // Step 3: Construct A and B
        A_mat = _didhetero_construct_A(intermed.G_ig, intermed.R_g,
                    intermed.Y_diff, data.mu_G_g[., id_gt], mu_R_col)

        B_mat = _didhetero_construct_B(A_mat, intermed.G_ig, intermed.R_g,
                    data.mu_G_g[., id_gt], mu_R_col,
                    data.mu_E_g_t[., id_gt], data.mu_F_g_t[., id_gt])

        // Store B (force independent copy to avoid pointer aliasing)
        data.A_g_t[id_gt] = &(1 * A_mat)
        data.B_g_t[id_gt] = &(1 * B_mat)
    }

    // === Call bandwidth selector on all (g,t) pairs ===
    bw_vec = _didhetero_bwselect_all(bwselect, data.B_g_t, data.Z,
                 data.zeval, kd0_Z, kd1_Z, data.kernel, data.uniformall)

    return(bw_vec)
}


// -----------------------------------------------------------------------------
// _didhetero_parse_user_gteval()
// Parse the user-supplied gteval numlist from Stata locals.
// -----------------------------------------------------------------------------
real matrix _didhetero_parse_user_gteval()
{
    string scalar gteval_str
    string rowvector gteval_tokens
    real scalar ngt_pairs, i
    real matrix gteval_user

    gteval_str = st_local("_gteval_user")
    ngt_pairs = strtoreal(st_local("_ngt_pairs"))

    gteval_tokens = tokens(gteval_str)
    gteval_user = J(ngt_pairs, 2, .)

    for (i = 1; i <= ngt_pairs; i++) {
        gteval_user[i, 1] = strtoreal(gteval_tokens[(i - 1) * 2 + 1])
        gteval_user[i, 2] = strtoreal(gteval_tokens[(i - 1) * 2 + 2])
    }

    return(gteval_user)
}


// -----------------------------------------------------------------------------
// _didhetero_validate_user_gteval()
// Reject user-supplied (g,t) pairs that are outside the automatically
// constructed identification domain for the current sample and options.
// -----------------------------------------------------------------------------
void _didhetero_validate_user_gteval()
{
    external struct DidHeteroData scalar _dh_data

    real scalar pretrend, i
    real matrix gteval_user
    string scalar duplicate_pairs, invalid_pairs, pair_label

    gteval_user = _didhetero_parse_user_gteval()
    pretrend = strtoreal(st_local("pretrend"))

    duplicate_pairs = didhetero_duplicate_gteval_pairs(gteval_user)
    invalid_pairs = ""
    for (i = 1; i <= rows(gteval_user); i++) {
        if (!didhetero_user_gteval_in_domain(
                _dh_data,
                gteval_user[i, 1],
                gteval_user[i, 2],
                _dh_data.anticipation,
                _dh_data.control_group,
                pretrend)) {
            pair_label = "(" + strofreal(gteval_user[i, 1], "%9.0g") + "," +
                strofreal(gteval_user[i, 2], "%9.0g") + ")"
            if (strlen(invalid_pairs) > 0) {
                invalid_pairs = invalid_pairs + " "
            }
            invalid_pairs = invalid_pairs + pair_label
        }
    }

    st_local("_gteval_duplicate", strofreal(strlen(duplicate_pairs) > 0, "%9.0g"))
    st_local("_gteval_duplicate_pairs", duplicate_pairs)
    st_local("_gteval_invalid", strofreal(strlen(invalid_pairs) > 0, "%9.0g"))
    st_local("_gteval_invalid_pairs", invalid_pairs)
}


// -----------------------------------------------------------------------------
// _didhetero_set_user_gteval()
// Parse user-specified gteval from Stata locals and override the auto-built
// gteval in the external _dh_data struct.
//
// Reads Stata locals:
//   _gteval_user  - space-separated numlist "g1 t1 g2 t2 ..."
//   _ngt_pairs    - number of (g,t) pairs
//
// Side effects:
//   Overrides _dh_data.gteval, _dh_data.num_gteval
//   May set _dh_data.uniformall = 0 when num_gteval == 1
//
// Paper ref: Section 4.2 (evaluation pair handling)
// -----------------------------------------------------------------------------
void _didhetero_set_user_gteval()
{
    external struct DidHeteroData scalar _dh_data

    real matrix gteval_user
    real scalar ngt_pairs

    ngt_pairs = strtoreal(st_local("_ngt_pairs"))
    gteval_user = _didhetero_parse_user_gteval()

    didhetero_validate_user_gteval(
        _dh_data,
        gteval_user,
        _dh_data.anticipation,
        _dh_data.control_group,
        strtoreal(st_local("pretrend")))

    _dh_data.gteval = gteval_user
    _dh_data.num_gteval = ngt_pairs

    // If single (g,t) pair, uniform over (g,t,z) degenerates to z-only
    if (ngt_pairs == 1) {
        _dh_data.uniformall = 0
    }
}


// -----------------------------------------------------------------------------
// didhetero_run_from_ado()
// Full estimation pipeline orchestration.
// Called from didhetero.ado after data preparation.
//
// Reads external globals _dh_data and _dh_kc set by didhetero_init_from_ado().
// Runs the complete pipeline:
//   1. Build gteval
//   2. Init core arrays
//   3. Stage 1 (GPS + OR + KDE)
//   4. BW selection pre-loop
//   5. Stage 2/3 (DR estimation + SE + bootstrap UCB)
//   6. Store results in external _dh_data
//
// Side effects:
//   Modifies external _dh_data with all estimation results.
//   Creates Stata matrices and scalars for result posting.
//
// Paper ref: Section 4.2 (full function)
// -----------------------------------------------------------------------------
void didhetero_run_from_ado()
{
    external struct DidHeteroData scalar _dh_data
    external struct DidHeteroKernelConsts scalar _dh_kc

    struct DidHeteroStage1Results scalar s1
    real colvector bw_vec, bw_manual
    real matrix est
    real scalar pretrend, uniformall, bstrap, seed
    string scalar bwselect, bw_str
    string rowvector bw_tokens
    real scalar i

    // === Read remaining Stata locals ===
    pretrend   = strtoreal(st_local("pretrend"))
    uniformall = strtoreal(st_local("uniform"))
    bstrap     = strtoreal(st_local("bstrap"))
    seed       = strtoreal(st_local("seed"))
    bwselect   = st_local("bwselect")

    // Parse manual bandwidth if specified
    bw_str = st_local("bw")
    if (bw_str != "") {
        bw_tokens = tokens(bw_str)
        bw_manual = J(cols(bw_tokens), 1, .)
        for (i = 1; i <= cols(bw_tokens); i++) {
            bw_manual[i] = strtoreal(bw_tokens[i])
        }
    }
    else {
        bw_manual = J(0, 1, .)
    }

    // =========================================================================
    // Step 1: Build (g,t) evaluation pairs
    // Paper ref: Section 4.2 (evaluation pair construction)
    // =========================================================================
    // Check if user already specified gteval (via _didhetero_set_user_gteval)
    if (_dh_data.num_gteval == . | _dh_data.num_gteval == 0) {
        printf("{txt}Building (g,t) evaluation pairs...\n")
        didhetero_build_gteval(_dh_data, _dh_data.anticipation,
            _dh_data.control_group, pretrend, uniformall)

        // Update uniformall in struct (ensure single (g,t) forces uniformall=0)
        // Note: Mata passes scalars by value, so any override within
        // didhetero_build_gteval() would not propagate back to this caller.
        // When only one (g,t) pair exists, uniform inference over (g,t,z)
        // degenerates to z-only, so explicitly force uniformall=0.
        _dh_data.uniformall = uniformall
        if (_dh_data.num_gteval == 1) {
            _dh_data.uniformall = 0
        }
    }
    else {
        didhetero_validate_user_gteval(
            _dh_data,
            _dh_data.gteval,
            _dh_data.anticipation,
            _dh_data.control_group,
            pretrend)

        printf("{txt}Using user-specified (g,t) evaluation pairs\n")
        // Still need to compute supp_g, supp_t, period1, geval, teval, gbar
        // for Stage 1 dispatch
        _dh_data.supp_g = sort(uniqrows(_dh_data.G), 1)
        _dh_data.supp_t = sort(uniqrows(_dh_data.t_vals), 1)
        _dh_data.period1 = _dh_data.t_vals[1]

        // geval: unique g values from gteval
        _dh_data.geval = sort(uniqrows(_dh_data.gteval[., 1]), 1)

        // gbar computation (Paper Section 2)
        if (sum(_dh_data.G :== 0) == 0) {
            _dh_data.gbar = max(_dh_data.supp_g)
        }
        else {
            _dh_data.gbar = .
        }
    }

    printf("{txt}  Found %g valid (g,t) pairs\n", _dh_data.num_gteval)

    // =========================================================================
    // Step 2: Initialize core estimation arrays
    // Paper ref: Section 4.2 (initialization)
    // =========================================================================
    didhetero_init_core_arrays(_dh_data, _dh_data.num_gteval)

    // =========================================================================
    // Step 3: Stage 1 — GPS + OR + KDE
    // Paper ref: Section 4.2.1 (Stage 1 parametric estimation)
    // =========================================================================
    printf("{txt}Stage 1: Parametric estimation (GPS + OR) + KDE...\n")
    s1 = didhetero_stage1_dispatch(_dh_data, _dh_data.gteval, _dh_data.geval,
             _dh_data.control_group, _dh_data.anticipation, _dh_data.zeval)

    // Store density estimates in data struct for SE computation
    _dh_data.kd0_Z = s1.kd0_Z

    printf("{txt}  GPS, OR, and density estimation complete\n")

    // Store Stage 1 results as externals for aggte_gt Pass 2
    external real matrix _dh_gps_mat
    external real matrix _dh_or_mat
    external real colvector _dh_kd1_Z

    _dh_gps_mat = s1.gps_mat
    _dh_or_mat  = s1.or_mat
    _dh_kd1_Z   = s1.kd1_Z

    // =========================================================================
    // Step 4: Bandwidth selection pre-loop
    // Paper ref: Section 4.2 (bandwidth selection)
    // =========================================================================
    printf("{txt}Bandwidth selection (%s)...\n", bwselect)
    bw_vec = didhetero_bw_preloop(_dh_data, s1.gps_mat, s1.or_mat,
                 s1.kd0_Z, s1.kd1_Z, bwselect, bw_manual)

    printf("{txt}  Bandwidths: ")
    for (i = 1; i <= rows(bw_vec); i++) {
        printf("%9.6f ", bw_vec[i])
        if (mod(i, 8) == 0 & i < rows(bw_vec)) printf("\n              ")
    }
    printf("\n")

    // =========================================================================
    // Step 5: Stage 2/3 — DR estimation + SE + analytical UCB + bootstrap UCB
    // Paper ref: Section 4.2 (DR estimation, SE, UCB)
    // =========================================================================
    printf("{txt}Stage 2/3: DR estimation")
    if (bstrap) {
        printf(" + bootstrap (%g iterations)", _dh_data.biters)
    }
    printf("...\n")

    est = didhetero_stage23(_dh_data, s1.gps_mat, s1.or_mat, bw_vec, bwselect, seed)

    printf("{txt}Estimation complete.\n")

    // =========================================================================
    // Step 6: Store results in Stata matrices
    // =========================================================================
    didhetero_post_results(_dh_data, est, bw_vec, bstrap)
}


// -----------------------------------------------------------------------------
// didhetero_post_results()
// Post estimation results to Stata e() and create result matrices.
//
// Creates:
//   Stata matrix e(Estimate_b) - point estimates (vectorized)
//   Stata matrix e(results)  - full results table (g, t, z, est, se, ci1_l, ci1_u, ci2_l, ci2_u, bw)
//   Stata matrix e(gteval)   - (g,t) evaluation pairs
//   Stata matrix e(zeval)    - z evaluation points
//   Stata matrix e(bw)       - bandwidths
//   Stata scalar e(N)        - sample size
//   Stata scalar e(num_gteval) - number of (g,t) pairs
//   Stata scalar e(num_zeval)  - number of z evaluation points
//
// Args:
//   data   - DidHeteroData struct with all results
//   est    - num_zeval x num_gteval point estimates
//   bw_vec - num_gteval x 1 bandwidths
//   bstrap - 1 if bootstrap was performed
//
// Paper ref: Section 4.2 (results assembly)
// -----------------------------------------------------------------------------
void didhetero_post_results(
    struct DidHeteroData scalar data,
    real matrix est,
    real colvector bw_vec,
    real scalar bstrap)
{
    real scalar K, R, id_gt, r, row
    real matrix results
    real scalar g1, t1

    K = data.num_gteval
    R = data.num_zeval

    // === Build results matrix ===
    // Columns: g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw
    // Rows: K * R (one row per (g,t,z) combination)
    results = J(K * R, 10, .)

    row = 0
    for (id_gt = 1; id_gt <= K; id_gt++) {
        g1 = data.gteval[id_gt, 1]
        t1 = data.gteval[id_gt, 2]

        for (r = 1; r <= R; r++) {
            row = row + 1
            results[row, 1]  = g1                           // g
            results[row, 2]  = t1                           // t
            results[row, 3]  = data.zeval[r]                // z
            results[row, 4]  = est[r, id_gt]                // est
            results[row, 5]  = data.se[r, id_gt]            // se
            results[row, 6]  = data.ci1_lower[r, id_gt]     // ci1_lower
            results[row, 7]  = data.ci1_upper[r, id_gt]     // ci1_upper
            if (bstrap) {
                results[row, 8]  = data.ci2_lower[r, id_gt] // ci2_lower
                results[row, 9]  = data.ci2_upper[r, id_gt] // ci2_upper
            }
            results[row, 10] = bw_vec[id_gt]                // bw
        }
    }

    // === Post to Stata ===
    // Results matrix
    st_matrix("e(results)", results)
    st_matrixcolstripe("e(results)",
        (J(10, 1, ""), ("g" \ "t" \ "z" \ "est" \ "se" \
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))

    // Keep a Stata-side alias `e(Estimate)` for backward compatibility.
    st_matrix("e(Estimate)", results)
    st_matrixcolstripe("e(Estimate)",
        (J(10, 1, ""), ("g" \ "t" \ "z" \ "est" \ "se" \
         "ci1_lower" \ "ci1_upper" \ "ci2_lower" \ "ci2_upper" \ "bw")))

    // Point estimates (vectorized: all z for gt1, then all z for gt2, ...)
    // Note: e(b) is reserved by Stata (requires ereturn post); use e(Estimate_b)
    st_matrix("e(Estimate_b)", vec(est)')

    // gteval matrix
    st_matrix("e(gteval)", data.gteval)
    st_matrixcolstripe("e(gteval)", (J(2, 1, ""), ("g" \ "t")))

    // zeval vector
    st_matrix("e(zeval)", data.zeval')

    // Bandwidths
    st_matrix("e(bw)", bw_vec')

    // Analytical critical values
    st_matrix("e(c_hat)", data.c_hat')

    // Bootstrap critical values (if available)
    if (bstrap & rows(data.c_check_bs) > 0) {
        st_matrix("e(c_check)", data.c_check_bs')
    }

    // === Store matrices needed by aggte_gt ===

    // B_g_t: flatten pointer array to n x (R * K) matrix
    // Each B_g_t[k] is n x R, concatenate horizontally
    {
        real matrix B_g_t_flat
        real scalar k2
        B_g_t_flat = J(data.n, data.num_zeval * K, .)
        for (k2 = 1; k2 <= K; k2++) {
            if (data.B_g_t[k2] != NULL) {
                B_g_t_flat[., ((k2-1)*R+1)..(k2*R)] = *data.B_g_t[k2]
            }
        }
        st_matrix("e(B_g_t)", B_g_t_flat)
    }

    // G_g: n x K group indicator matrix
    st_matrix("e(G_g)", data.G_g)

    // Z: n x 1 covariate vector (store as row for Stata compatibility)
    st_matrix("e(Z)", data.Z')

    // Persist the minimum Stage-0 / Stage-1 inputs needed by aggte_gt Pass 2.
    // These matrices make the post-estimation object self-contained even after
    // `mata clear` removes external globals from the current session.
    st_matrix("e(dh_Y_wide)", data.Y_wide)
    st_matrix("e(dh_G_unit)", data.G')
    st_matrix("e(dh_t_vals)", data.t_vals')

    {
        external real matrix _dh_gps_mat
        external real matrix _dh_or_mat
        st_matrix("e(dh_gps_mat)", _dh_gps_mat)
        st_matrix("e(dh_or_mat)", _dh_or_mat)
    }

    // mu_G_g: R x K conditional group density
    st_matrix("e(mu_G_g)", data.mu_G_g)

    // catt_est: R x K point estimates (same as est parameter)
    st_matrix("e(catt_est)", est)

    // catt_se: R x K standard errors
    st_matrix("e(catt_se)", data.se)

    // kd0_Z: R x 1 density estimates
    st_matrix("e(kd0_Z)", data.kd0_Z')

    // kd1_Z: R x 1 density derivative estimates (from external)
    {
        external real colvector _dh_kd1_Z
        st_matrix("e(kd1_Z)", _dh_kd1_Z')
    }

    // Z_supp: support points
    st_matrix("e(Z_supp)", data.Z_supp')

    // gbar: comparison-set upper bound on treatment timing.
    // Public contract: when never-treated units are present, Stata stores
    // missing (.) as the numeric sentinel for +Inf and exposes e(gbar_isinf)=1.
    st_numscalar("e(gbar)", data.gbar)
    st_numscalar("e(gbar_isinf)", missing(data.gbar))

    // Scalars
    st_numscalar("e(N)", data.n)
    st_numscalar("e(num_gteval)", data.num_gteval)
    st_numscalar("e(num_zeval)", data.num_zeval)
    st_numscalar("e(T)", data.T_num)
}

end
