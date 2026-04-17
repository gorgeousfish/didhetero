mata:

// =============================================================================
// didhetero_intermediate.mata
// Intermediate variable construction for DR estimator
//
// Constructs R_{i,g,t}, Y_tilde, E_{i,g,t}, F_{i,g,t} for each (g,t) pair.
//
// Functions:
//   1. _didhetero_extract_gps()       - Extract GPS for a specific (g,[t]) pair
//   2. _didhetero_extract_or()        - Extract OR for a specific (g,t) pair
//   3. didhetero_intermediate_vars()  - Main: construct all intermediate vars
//
// Structs:
//   DidHeteroIntermediate - Container for intermediate variables
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//   Section 4.2.1 (intermediate variable construction)
// =============================================================================

// -----------------------------------------------------------------------------
// DidHeteroIntermediate
// Container for intermediate variables at a single (g,t) pair.
//
// All vectors are n x 1 (one entry per cross-sectional unit).
//
// Paper ref: Section 4.2.1, Eq. 14 (intermediate variable definitions)
// -----------------------------------------------------------------------------
struct DidHeteroIntermediate {
    real colvector R_g      // n x 1, R_{i,g,t} = p_hat * comp_ind / (1 - p_hat)
    real colvector Y_diff   // n x 1, Y_tilde_{i,g,t} = Y_t - Y_base - m_hat
    real colvector E_g_t    // n x 1, E_{i,g,t} = R_g * Y_diff
    real colvector F_g_t    // n x 1, F_{i,g,t} = G_ig * Y_diff
    real colvector G_ig     // n x 1, G_{i,g} indicator (1 if unit i in group g)
}

// -----------------------------------------------------------------------------
// _didhetero_extract_gps()
// Extract GPS estimates for a specific group (and time, if notyettreated).
//
// GPS matrix format:
//   nevertreated:  3 columns (id, g, est)
//   notyettreated: 4 columns (id, g, t, est)
//
// Args:
//   gps_mat       - GPS result matrix from didhetero_gps_estimate()
//   g1            - target group value
//   t1            - target time value (used only for notyettreated)
//   control_group - "nevertreated" or "notyettreated"
//   n             - expected number of rows in output
//
// Returns:
//   n x 1 vector of GPS estimates p_hat
//
// Paper ref: Section 4.2.1 (GPS extraction for DR estimand)
// -----------------------------------------------------------------------------
real colvector _didhetero_extract_gps(
    real matrix gps_mat,
    real scalar g1,
    real scalar t1,
    string scalar control_group,
    real scalar n)
{
    real colvector mask, p_hat
    real rowvector idx
    real scalar n_extracted

    if (control_group == "nevertreated") {
        // GPS matrix: (id, g, est) — 3 columns
        // Filter: column 2 == g1
        mask = (gps_mat[., 2] :== g1)
        idx = didhetero_selectindex(mask)

        if (cols(idx) == 0) {
            _error("GPS extract: no rows found for g=" + strofreal(g1))
        }

        p_hat = gps_mat[idx', 3]
    }
    else {
        // GPS matrix: (id, g, t, est) — 4 columns
        // Filter: column 2 == g1 AND column 3 == t1
        mask = (gps_mat[., 2] :== g1) :& (gps_mat[., 3] :== t1)
        idx = didhetero_selectindex(mask)

        if (cols(idx) == 0) {
            _error("GPS extract: no rows found for g=" + strofreal(g1) +
                   " t=" + strofreal(t1))
        }

        p_hat = gps_mat[idx', 4]
    }

    // Validate row count
    n_extracted = rows(p_hat)
    if (n_extracted != n) {
        _error("GPS extract: expected " + strofreal(n) + " rows but got " +
               strofreal(n_extracted) + " for g=" + strofreal(g1) +
               " t=" + strofreal(t1))
    }

    return(p_hat)
}

// -----------------------------------------------------------------------------
// _didhetero_extract_or()
// Extract OR predicted values for a specific (g,t) pair.
//
// OR matrix format: always 4 columns (id, g, t, est)
//
// Args:
//   or_mat - OR result matrix from didhetero_or_estimate()
//   g1     - target group value
//   t1     - target time value
//   n      - expected number of rows in output
//
// Returns:
//   n x 1 vector of OR predicted values m_hat
//
// Paper ref: Section 4.2.1 (OR extraction for DR estimand)
// -----------------------------------------------------------------------------
real colvector _didhetero_extract_or(
    real matrix or_mat,
    real scalar g1,
    real scalar t1,
    real scalar n)
{
    real colvector mask, m_hat
    real rowvector idx
    real scalar n_extracted

    // OR matrix: (id, g, t, est) — always 4 columns
    // Filter: column 2 == g1 AND column 3 == t1
    mask = (or_mat[., 2] :== g1) :& (or_mat[., 3] :== t1)
    idx = didhetero_selectindex(mask)

    if (cols(idx) == 0) {
        _error("OR extract: no rows found for g=" + strofreal(g1) +
               " t=" + strofreal(t1))
    }

    m_hat = or_mat[idx', 4]

    // Validate row count
    n_extracted = rows(m_hat)
    if (n_extracted != n) {
        _error("OR extract: expected " + strofreal(n) + " rows but got " +
               strofreal(n_extracted) + " for g=" + strofreal(g1) +
               " t=" + strofreal(t1))
    }

    return(m_hat)
}

// -----------------------------------------------------------------------------
// didhetero_intermediate_vars()
// Construct all intermediate variables for a single (g,t) pair.
//
// Implements the DR estimator building blocks:
//   R_{i,g,t}     = p_hat * comp_ind / (1 - p_hat)
//   Y_tilde       = Y_t - Y_{g-delta-1} - m_hat
//   E_{i,g,t}     = R_{i,g,t} * Y_tilde
//   F_{i,g,t}     = G_{i,g} * Y_tilde
//
// Args:
//   data          - DidHeteroData struct (panel data)
//   gps_mat       - GPS result matrix from Stage 1
//   or_mat        - OR result matrix from Stage 1
//   g1            - target group value
//   t1            - target time value
//   id_gt         - 1-based index of this (g,t) pair in gteval
//   control_group - "nevertreated" or "notyettreated"
//   anticipation  - anticipation parameter delta >= 0
//   G_g           - n x num_gteval matrix (modified in place: column id_gt set)
//
// Returns:
//   DidHeteroIntermediate scalar with all fields populated
//
// Paper ref: Section 4.2.1, Eq. 14 (full intermediate variable construction)
// -----------------------------------------------------------------------------
struct DidHeteroIntermediate scalar didhetero_intermediate_vars(
    struct DidHeteroData scalar data,
    real matrix gps_mat,
    real matrix or_mat,
    real scalar g1,
    real scalar t1,
    real scalar id_gt,
    string scalar control_group,
    real scalar anticipation,
    real matrix G_g)
{
    struct DidHeteroIntermediate scalar result
    real colvector p_hat, G_ig, comp_ind, R_g, G_ord
    real colvector Y_t, Y_base, m_hat, Y_diff
    real scalar col_t, col_base, threshold_ord, base_label, i

    // =========================================================================
    // Step 1: Extract GPS estimates for this (g,t) pair
    // Paper ref: Section 4.2.1, GPS estimation
    // =========================================================================
    p_hat = _didhetero_extract_gps(gps_mat, g1, t1, control_group, data.n)

    // =========================================================================
    // Step 2: Construct indicator variables
    // Paper ref: Section 4.2.1, indicator variables
    // =========================================================================

    // G_{i,g}: treatment group indicator
    G_ig = (data.G :== g1)

    if (control_group == "nevertreated") {
        // Comparison group: never-treated only
        comp_ind = (data.G :== 0)
    }
    else {
        // Comparison group: never-treated + units untreated at t + delta.
        threshold_ord = didhetero_period_ord(t1, data.t_vals) + anticipation
        G_ord = J(rows(data.G), 1, 0)
        for (i = 1; i <= rows(data.G); i++) {
            if (data.G[i] != 0) {
                G_ord[i] = didhetero_period_ord(data.G[i], data.t_vals)
            }
        }

        // comp_ind = 1*(G == 0) + 1*(G > threshold)
        // These groups are mutually exclusive so addition == logical OR numerically.
        comp_ind = (data.G :== 0) + (G_ord :> threshold_ord)
    }

    // =========================================================================
    // Step 3: Construct R_{i,g,t}
    // Paper ref: Section 4.2.1, Eq. 14
    // R_{i,g,t} = p_hat * comp_ind / (1 - p_hat)
    // =========================================================================
    R_g = (p_hat :* comp_ind) :/ (1 :- p_hat)

    // =========================================================================
    // Step 4: Store G_{i,g} into G_g matrix (passed by reference)
    // Store G_{i,g} column for later use in influence functions
    // =========================================================================
    G_g[., id_gt] = G_ig

    // =========================================================================
    // Step 5: Extract outcome variables Y_t and Y_{g-delta-1}
    // Paper ref: Section 4.2.1 (outcome differencing)
    // =========================================================================
    col_t = didhetero_period_col(t1, data.t_vals)
    base_label = didhetero_period_at(
        didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
        data.t_vals,
        "intermediate vars base period")
    col_base = didhetero_period_col(base_label, data.t_vals)

    Y_t    = data.Y_wide[., col_t]
    Y_base = data.Y_wide[., col_base]

    // =========================================================================
    // Step 6: Extract OR predicted values m_hat
    // Paper ref: Section 4.2.1, OR prediction
    // =========================================================================
    m_hat = _didhetero_extract_or(or_mat, g1, t1, data.n)

    // =========================================================================
    // Step 7: Construct Y_tilde = Y_t - Y_{g-delta-1} - m_hat
    // Y_tilde = Y_t - Y_{g-delta-1} - m_hat
    // =========================================================================
    Y_diff = Y_t - Y_base - m_hat

    // =========================================================================
    // Step 8: Construct E_{i,g,t} and F_{i,g,t}
    // Paper ref: Section 4.2.1, Eq. 14
    // E_{i,g,t} = R_{i,g,t} * Y_tilde
    // F_{i,g,t} = G_{i,g} * Y_tilde
    // =========================================================================

    // =========================================================================
    // Pack results into struct
    // =========================================================================
    result.R_g    = R_g
    result.Y_diff = Y_diff
    result.E_g_t  = R_g :* Y_diff       // E_{i,g,t} = R_{i,g,t} * Y_tilde
    result.F_g_t  = G_ig :* Y_diff       // F_{i,g,t} = G_{i,g} * Y_tilde
    result.G_ig   = G_ig

    return(result)
}

end
