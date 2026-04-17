// =============================================================================
// didhetero_simdata.mata — Core data generation for didhetero_simdata command
//
// Implements the DGP described in Section 6 of Imai, Qin, and Yanagi (2025)
//
// Reference: Imai, Qin, and Yanagi (2025), Section 6
// Group assignment builds on Callaway and Sant'Anna (2021).
// =============================================================================

mata:
mata set matastrict on

// -----------------------------------------------------------------------------
// didhetero_simdata_core()
//
// Args:
//   n          : number of cross-sectional units
//   tau        : number of time periods (tau > 2)
//   hc         : 1 = heteroscedastic errors, 0 = homoscedastic
//   dimx       : dimension of covariates X (>= 1)
//   dgpy       : DGP for treated potential outcome (1 or 2)
//   continuous : 1 = continuous Z, 0 = discrete Z in {-1, 0, 1}
//
// Returns:
//   real matrix (n*tau) x (4 + dimx) with columns:
//     id, period, Y, G, Z [, X1, X2, ...]
// -----------------------------------------------------------------------------
real matrix didhetero_simdata_core(
    real scalar n,
    real scalar tau,
    real scalar hc,
    real scalar dimx,
    real scalar dgpy,
    real scalar continuous)
{
    real colvector Z, G, eta, u_sd, v_sd, u, v
    real colvector Y_0, Y_g, Y, mgt, delta_e, treated
    real matrix X, group_P, data
    real rowvector group_supp, gamma_vec, beta_t0, cumprobs
    real scalar i, t, delta_t, k, j
    real colvector log_probs_i, probs_i
    real scalar U_draw

    // -----------------------------------------------------------------
    // Step 1: Pre-treatment covariate Z
    // -----------------------------------------------------------------
    if (continuous) {
        Z = rnormal(n, 1, 0, 1)
    }
    else {
        // Discrete Z in {-1, 0, 1} with equal probability
        real colvector U_disc
        U_disc = runiform(n, 1)
        Z = J(n, 1, 1)
        for (i = 1; i <= n; i++) {
            if (U_disc[i] < 1/3) {
                Z[i] = -1
            }
            else if (U_disc[i] < 2/3) {
                Z[i] = 0
            }
            // else Z[i] = 1 (already set)
        }
    }

    // -----------------------------------------------------------------
    // Step 2: All pre-treatment covariates X
    // -----------------------------------------------------------------
    X = Z
    if (dimx > 1) {
        for (k = 1; k <= (dimx - 1); k++) {
            X = X, rnormal(n, 1, 0, 1)
        }
    }

    // -----------------------------------------------------------------
    // Step 3: Group support and multinomial logit assignment
    // -----------------------------------------------------------------
    // group_supp = (0, 2, 3, ..., tau)
    group_supp = 0, (2..tau)

    // Compute group choice probabilities via multinomial logit
    // gamma = 0.5 * g / tau for each g in group_supp
    gamma_vec = 0.5 :* group_supp :/ tau

    // group_P: n x length(group_supp) matrix of probabilities
    group_P = J(n, cols(group_supp), .)
    for (i = 1; i <= n; i++) {
        log_probs_i = (Z[i] :* gamma_vec)'
        // Numerical stability: subtract max
        log_probs_i = log_probs_i :- max(log_probs_i)
        probs_i = exp(log_probs_i)
        probs_i = probs_i :/ sum(probs_i)
        group_P[i, .] = probs_i'
    }

    // Sample group for each individual
    G = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        U_draw = runiform(1, 1)
        cumprobs = runningsum(group_P[i, .])
        G[i] = group_supp[cols(group_supp)]  // default: last group
        for (j = 1; j <= cols(group_supp); j++) {
            if (U_draw <= cumprobs[j]) {
                G[i] = group_supp[j]
                break
            }
        }
    }

    // -----------------------------------------------------------------
    // Step 4: Individual effect
    // -----------------------------------------------------------------
    eta = rnormal(1, 1, G, J(n, 1, 1))

    // -----------------------------------------------------------------
    // Step 5: Heteroscedasticity setup
    // -----------------------------------------------------------------
    if (hc) {
        u_sd = 0.5 :+ normal(Z)
        v_sd = G :/ tau :+ normal(Z)
    }
    else {
        u_sd = J(n, 1, 1)
        v_sd = J(n, 1, 1)
    }

    // -----------------------------------------------------------------
    // Step 6: Data generation loop over time periods
    // -----------------------------------------------------------------
    data = J(n * tau, 4 + dimx, 0)

    for (t = 1; t <= tau; t++) {

        // Error term for untreated potential outcome
        u = rnormal(n, 1, 0, 1) :* u_sd

        // Untreated potential outcome
        delta_t = t
        beta_t0 = J(1, dimx, .)
        for (k = 1; k <= dimx; k++) {
            beta_t0[k] = t / k
        }

        if (dimx == 1) {
            Y_0 = J(n, 1, delta_t) :+ eta :+ X :* beta_t0 :+ u
        }
        else {
            Y_0 = J(n, 1, delta_t) :+ eta :+ X * beta_t0' :+ u
        }

        // Error term for treated potential outcome
        v = rnormal(n, 1, 0, 1) :* v_sd

        // Treatment effect term
        delta_e = J(n, 1, t) :- G :+ J(n, 1, 1)

        if (dgpy == 1) {
            mgt = (G :/ t) :* sin(pi() :* Z)
        }
        else {
            // DGPY == 2
            mgt = Z :* G :/ t
        }

        // Treated potential outcome
        Y_g = Y_0 :+ mgt :+ delta_e :+ v :- u

        // Observed outcome: switching equation
        treated = (G :!= 0) :& (G :<= t)
        Y = treated :* Y_g :+ (1 :- treated) :* Y_0

        // Fill data block for this period
        // Columns: id, period, Y, G, Z [, X1, X2, ...]
        data[((t-1)*n+1)::(t*n), 1] = (1::n)
        data[((t-1)*n+1)::(t*n), 2] = J(n, 1, t)
        data[((t-1)*n+1)::(t*n), 3] = Y
        data[((t-1)*n+1)::(t*n), 4] = G
        data[((t-1)*n+1)::(t*n), 5] = Z
        if (dimx > 1) {
            for (k = 2; k <= dimx; k++) {
                data[((t-1)*n+1)::(t*n), 4 + k] = X[., k]
            }
        }
    }

    // Sort by id, then period
    data = sort(data, (1, 2))

    return(data)
}

// -----------------------------------------------------------------------------
// _didhetero_simdata_to_stata()
//
// Bridge function: calls didhetero_simdata_core(), then writes the result
// matrix into Stata variables using st_addobs() and st_store().
// -----------------------------------------------------------------------------
void _didhetero_simdata_to_stata(
    real scalar n,
    real scalar tau,
    real scalar hc,
    real scalar dimx,
    real scalar dgpy,
    real scalar continuous)
{
    real matrix data
    real scalar nobs, k, idx

    // Generate data
    data = didhetero_simdata_core(n, tau, hc, dimx, dgpy, continuous)
    nobs = rows(data)

    // Add observations
    st_addobs(nobs)

    // Create variables and store data
    // Column order: id, period, Y, G, Z [, X1, X2, ...]
    idx = st_addvar("long",   "id")
    st_store(., idx, data[., 1])

    idx = st_addvar("long",   "period")
    st_store(., idx, data[., 2])

    idx = st_addvar("double", "Y")
    st_store(., idx, data[., 3])

    idx = st_addvar("long",   "G")
    st_store(., idx, data[., 4])

    idx = st_addvar("double", "Z")
    st_store(., idx, data[., 5])

    if (dimx > 1) {
        for (k = 2; k <= dimx; k++) {
            idx = st_addvar("double", "X" + strofreal(k - 1))
            st_store(., idx, data[., 4 + k])
        }
    }
}

end
