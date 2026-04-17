mata:

// =============================================================================
// Intermediate variable construction for doubly robust estimation
//
// This module constructs the intermediate quantities required for computing
// CATT estimates: the weighting variable R, the adjusted outcome difference
// Y_tilde, and the products E and F.
//
// Functions:
//   _didhetero_extract_gps()       - Extract GPS estimates for a given (g,t)
//   _didhetero_extract_or()        - Extract outcome regression predictions
//   didhetero_intermediate_vars()  - Construct all intermediate variables
//
// Structure:
//   DidHeteroIntermediate - Container for intermediate variables
// =============================================================================

// -----------------------------------------------------------------------------
// DidHeteroIntermediate
// Container for intermediate variables at a single (g,t) pair.
// All vectors are n x 1, where n is the number of cross-sectional units.
// -----------------------------------------------------------------------------
struct DidHeteroIntermediate {
    real colvector R_g      // n x 1, weighting variable R_{i,g,t}
    real colvector Y_diff   // n x 1, adjusted outcome difference Y_tilde
    real colvector E_g_t    // n x 1, product E_{i,g,t} = R_{i,g,t} * Y_tilde
    real colvector F_g_t    // n x 1, product F_{i,g,t} = G_{i,g} * Y_tilde
    real colvector G_ig     // n x 1, group membership indicator G_{i,g}
}

// -----------------------------------------------------------------------------
// _didhetero_extract_gps()
// Extract generalized propensity score estimates for a specific (g,t) pair.
//
// The GPS matrix has distinct structures depending on the control group:
//   - nevertreated:  3 columns (id, g, estimate)
//   - notyettreated: 4 columns (id, g, t, estimate)
//
// Arguments:
//   gps_mat       - matrix containing GPS estimation results
//   g1            - target group identifier
//   t1            - target time period (relevant for not-yet-treated)
//   control_group - control group type ("nevertreated" or "notyettreated")
//   n             - expected number of observations
//
// Returns:
//   n x 1 vector of GPS estimates
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
        // GPS matrix with 3 columns: (id, g, estimate)
        // Select rows where the group column matches g1
        mask = (gps_mat[., 2] :== g1)
        idx = didhetero_selectindex(mask)

        if (cols(idx) == 0) {
            _error("GPS extract: no rows found for g=" + strofreal(g1))
        }

        p_hat = gps_mat[idx', 3]
    }
    else {
        // GPS matrix with 4 columns: (id, g, t, estimate)
        // Select rows where group equals g1 and time equals t1
        mask = (gps_mat[., 2] :== g1) :& (gps_mat[., 3] :== t1)
        idx = didhetero_selectindex(mask)

        if (cols(idx) == 0) {
            _error("GPS extract: no rows found for g=" + strofreal(g1) +
                   " t=" + strofreal(t1))
        }

        p_hat = gps_mat[idx', 4]
    }

    // Verify the number of extracted rows matches expected count
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
// Extract outcome regression predictions for a specific (g,t) pair.
//
// The OR matrix has 4 columns: (id, g, t, prediction).
//
// Arguments:
//   or_mat - matrix containing outcome regression predictions
//   g1     - target group identifier
//   t1     - target time period
//   n      - expected number of observations
//
// Returns:
//   n x 1 vector of predicted values
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

    // OR matrix with 4 columns: (id, g, t, prediction)
    // Select rows where group equals g1 and time equals t1
    mask = (or_mat[., 2] :== g1) :& (or_mat[., 3] :== t1)
    idx = didhetero_selectindex(mask)

    if (cols(idx) == 0) {
        _error("OR extract: no rows found for g=" + strofreal(g1) +
               " t=" + strofreal(t1))
    }

    m_hat = or_mat[idx', 4]

    // Verify the number of extracted rows matches expected count
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
// This function computes the components required for CATT estimation:
//   R_{i,g,t} = p_hat * comp_ind / (1 - p_hat)
//   Y_tilde   = Y_t - Y_{g-delta-1} - m_hat
//   E_{i,g,t} = R_{i,g,t} * Y_tilde
//   F_{i,g,t} = G_{i,g} * Y_tilde
//
// Arguments:
//   data          - panel data structure
//   gps_mat       - matrix of GPS estimates
//   or_mat        - matrix of outcome regression predictions
//   g1            - target group identifier
//   t1            - target time period
//   id_gt         - index of the (g,t) pair in the evaluation set
//   control_group - control group type ("nevertreated" or "notyettreated")
//   anticipation  - anticipation period (delta >= 0)
//   G_g           - matrix storing group indicators (modified in place)
//
// Returns:
//   structure containing all intermediate variables
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

    // Extract GPS estimates for the current (g,t) pair
    p_hat = _didhetero_extract_gps(gps_mat, g1, t1, control_group, data.n)

    // Construct group membership and comparison group indicators

    // G_{i,g}: indicator for membership in treatment group g
    G_ig = (data.G :== g1)

    if (control_group == "nevertreated") {
        // Comparison group consists of never-treated units only
        comp_ind = (data.G :== 0)
    }
    else {
        // Comparison group includes never-treated and not-yet-treated units
        threshold_ord = didhetero_period_ord(t1, data.t_vals) + anticipation
        G_ord = J(rows(data.G), 1, 0)
        for (i = 1; i <= rows(data.G); i++) {
            if (data.G[i] != 0) {
                G_ord[i] = didhetero_period_ord(data.G[i], data.t_vals)
            }
        }

        // Comparison indicator: never-treated or not-yet-treated by t + delta
        // The two conditions are mutually exclusive
        comp_ind = (data.G :== 0) + (G_ord :> threshold_ord)
    }

    // Compute the weighting variable R_{i,g,t}
    R_g = (p_hat :* comp_ind) :/ (1 :- p_hat)

    // Store group indicator in the G_g matrix for subsequent computations
    G_g[., id_gt] = G_ig

    // Extract outcome variables for periods t and g-delta-1
    col_t = didhetero_period_col(t1, data.t_vals)
    base_label = didhetero_period_at(
        didhetero_period_ord(g1, data.t_vals) - anticipation - 1,
        data.t_vals,
        "intermediate vars base period")
    col_base = didhetero_period_col(base_label, data.t_vals)

    Y_t    = data.Y_wide[., col_t]
    Y_base = data.Y_wide[., col_base]

    // Extract outcome regression predictions
    m_hat = _didhetero_extract_or(or_mat, g1, t1, data.n)

    // Compute the adjusted outcome difference Y_tilde
    Y_diff = Y_t - Y_base - m_hat

    // Assemble results into the output structure
    result.R_g    = R_g
    result.Y_diff = Y_diff
    result.E_g_t  = R_g :* Y_diff
    result.F_g_t  = G_ig :* Y_diff
    result.G_ig   = G_ig

    return(result)
}

end
