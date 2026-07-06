mata:

// =============================================================================
// didhetero_utils_domain.mata
// Period/group validation and (g,t) domain utility functions
//
// Provides:
//   - didhetero_period_ord()           // Time label -> ordinal position
//   - didhetero_period_at()            // Ordinal position -> time label
//   - didhetero_group_is_valid()       // Validate group under anticipation
//   - didhetero_time_is_valid()        // Validate time under anticipation
//   - didhetero_pair_in_domain()       // Check (g,t) pair admissibility
//   - didhetero_build_gteval()         // Build full (g,t) evaluation set
//   - didhetero_duplicate_gteval_pairs() // Detect duplicate pairs
//   - didhetero_user_gteval_in_domain()  // Check user pair in domain
//   - didhetero_validate_user_gteval()   // Validate user-specified gteval
//
// Requires:
//   - didhetero_types.mata      (DidHeteroData struct)
//   - didhetero_errors.mata     (_dh_error_input)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_period_ord()
// Return the 1-based ordinal position of a time label in the observed support.
//
// Args:
//   t      - time label
//   t_vals - sorted observed time support
//
// Returns:
//   1-based ordinal position, or missing if the label is not observed
// -----------------------------------------------------------------------------
real scalar didhetero_period_ord(real scalar t, real colvector t_vals)
{
    real scalar i

    for (i = 1; i <= rows(t_vals); i++) {
        if (t_vals[i] == t) return(i)
    }

    return(.)
}

// -----------------------------------------------------------------------------
// didhetero_period_at()
// Return the observed time label at a given 1-based ordinal position.
//
// Args:
//   ord     - 1-based ordinal position
//   t_vals  - sorted observed time support
//   context - caller label for error reporting
// -----------------------------------------------------------------------------
real scalar didhetero_period_at(real scalar ord, real colvector t_vals, string scalar context)
{
    if (ord < 1 | ord > rows(t_vals)) {
        _error(context + ": ordinal period index " + strofreal(ord) +
            " is outside the observed time support")
    }

    return(t_vals[ord])
}

// -----------------------------------------------------------------------------
// didhetero_group_is_valid()
// Check whether a treatment group has a valid unaffected base period under the
// current anticipation horizon.
// -----------------------------------------------------------------------------
real scalar didhetero_group_is_valid(real scalar g1, real scalar anticipation, real colvector t_vals)
{
    real scalar g_ord

    if (g1 == 0) return(0)

    g_ord = didhetero_period_ord(g1, t_vals)
    if (g_ord >= .) return(0)

    return(g_ord > anticipation + 1)
}

// -----------------------------------------------------------------------------
// didhetero_time_is_valid()
// Check whether an evaluation time is within the observed support and leaves
// enough forward periods for the anticipation horizon.
// -----------------------------------------------------------------------------
real scalar didhetero_time_is_valid(real scalar t1, real scalar anticipation, real colvector t_vals)
{
    real scalar t_ord

    t_ord = didhetero_period_ord(t1, t_vals)
    if (t_ord >= .) return(0)

    return(t_ord > 1 & t_ord <= rows(t_vals) - anticipation)
}

// -----------------------------------------------------------------------------
// didhetero_comp_threshold_ord()
// Ordinal period threshold defining the not-yet-treated comparison pool for a
// (g,t) pair.  Post-treatment pairs use the pool untreated through t + delta
// (R_{g,t} in eq. 2.6 of Imai, Qin, and Yanagi).  Pre-treatment pairs
// (t < g - delta) anchor the pool at t = g (R_{g,g} in the pre-treatment DR
// estimand of their Appendix F), which keeps cohorts treated within the
// pre-trend window out of the comparison pool.
// -----------------------------------------------------------------------------
real scalar didhetero_comp_threshold_ord(
    real scalar g1,
    real scalar t1,
    real scalar anticipation,
    real colvector t_vals)
{
    real scalar g_ord, t_ord

    g_ord = didhetero_period_ord(g1, t_vals)
    t_ord = didhetero_period_ord(t1, t_vals)

    if (t_ord < g_ord - anticipation) {
        return(g_ord + anticipation)
    }
    return(t_ord + anticipation)
}

// -----------------------------------------------------------------------------
// didhetero_pair_in_domain()
// Check whether a (g,t) pair is in the estimable domain when time labels are
// interpreted through their ordinal positions in the observed support.
// -----------------------------------------------------------------------------
real scalar didhetero_pair_in_domain(
    real scalar g1,
    real scalar t1,
    real scalar anticipation,
    string scalar control_group,
    real scalar pretrend,
    real scalar has_never,
    real scalar gbar,
    real colvector t_vals)
{
    real scalar g_ord, t_ord, gbar_ord

    if (!didhetero_group_is_valid(g1, anticipation, t_vals)) {
        return(0)
    }
    if (!didhetero_time_is_valid(t1, anticipation, t_vals)) {
        return(0)
    }

    g_ord = didhetero_period_ord(g1, t_vals)
    t_ord = didhetero_period_ord(t1, t_vals)

    if (pretrend != 0) {
        if (t_ord == g_ord - anticipation - 1) return(0)
        if (t_ord < g_ord - anticipation &
            control_group == "notyettreated" & !has_never) {
            gbar_ord = didhetero_period_ord(gbar, t_vals)
            if (gbar_ord >= .) return(0)
            return(gbar_ord > g_ord + anticipation)
        }
        return(1)
    }

    if (t_ord < g_ord - anticipation) {
        return(0)
    }

    if (control_group == "nevertreated") {
        return(has_never)
    }

    if (control_group == "notyettreated") {
        if (has_never) {
            return(1)
        }

        gbar_ord = didhetero_period_ord(gbar, t_vals)
        if (gbar_ord >= .) return(0)
        return(t_ord < gbar_ord - anticipation)
    }

    _error(198, "Invalid control_group: " + control_group)
    return(0)
}

// -----------------------------------------------------------------------------
// didhetero_build_gteval()
// Build geval, teval, gteval evaluation point matrices
//
// Determines which (g,t) combinations need CATT estimation.
// Modifies data struct in place and may override uniformall.
//
// Args:
//   data          - DidHeteroData struct (modified in place)
//   anticipation  - anticipation parameter delta >= 0
//   control_group - "nevertreated" or "notyettreated"
//   pretrend      - 0 or 1
//   uniformall    - 0 or 1 (may be overridden when num_gteval == 1)
//
// Side effects:
//   Fills data.supp_g, data.supp_t, data.period1, data.geval, data.teval,
//   data.gbar, data.gteval, data.num_gteval
//   May set uniformall = 0 when num_gteval == 1
// -----------------------------------------------------------------------------
void didhetero_build_gteval(
    struct DidHeteroData scalar data,
    real scalar anticipation,
    string scalar control_group,
    real scalar pretrend,
    real scalar uniformall)
{
    real colvector supp_g, supp_t, geval, teval
    real scalar period1, gbar, num_gteval, has_never
    real scalar i, j, g1, t1
    real matrix gteval
    
    // === Step 1: Basic quantities ===
    supp_g = sort(uniqrows(data.G), 1)
    supp_t = sort(uniqrows(data.t_vals), 1)
    period1 = supp_t[1]

    // Guard the support restriction even if callers bypass ado-level
    // validation. Treated cohorts must start strictly after the first period.
    if (sum((supp_g :!= 0) :& (supp_g :<= period1)) > 0) {
        _error(198,
            "group() must be 0 for never-treated units or a treatment time strictly after the first observed period")
    }
    
    has_never = sum(data.G :== 0) > 0

    // === Step 2: geval computation on ordinal time positions ===
    geval = J(0, 1, .)
    for (i = 1; i <= rows(supp_g); i++) {
        g1 = supp_g[i]
        if (didhetero_group_is_valid(g1, anticipation, supp_t)) {
            geval = geval \ g1
        }
    }
    
    if (rows(geval) == 0) {
        _error("No valid treatment groups found after applying anticipation filter")
    }
    
    // === Step 3: teval computation on ordinal time positions ===
    teval = J(0, 1, .)
    for (i = 1; i <= rows(supp_t); i++) {
        t1 = supp_t[i]
        if (didhetero_time_is_valid(t1, anticipation, supp_t)) {
            teval = teval \ t1
        }
    }
    
    if (rows(teval) == 0) {
        _error("No valid time periods found after applying anticipation filter")
    }
    
    // === Step 4: gbar computation ===
    if (!has_never) {
        gbar = max(supp_g)
    }
    else {
        gbar = .
    }
    
    // === Step 5: gteval construction ===
    gteval = J(0, 2, .)
    
    for (i = 1; i <= rows(geval); i++) {
        g1 = geval[i]
        for (j = 1; j <= rows(teval); j++) {
            t1 = teval[j]
            if (didhetero_pair_in_domain(
                    g1, t1, anticipation, control_group,
                    pretrend, has_never, gbar, supp_t)) {
                gteval = gteval \ (g1, t1)
            }
        }
    }
    
    if (rows(gteval) == 0) {
        _error("No valid (g,t) pairs found for estimation")
    }
    
    // === Step 6: num_gteval and uniformall override ===
    num_gteval = rows(gteval)
    
    if (num_gteval == 1) {
        uniformall = 0
    }
    
    // === Step 7: Store to struct ===
    data.supp_g = supp_g
    data.supp_t = supp_t
    data.period1 = period1
    data.geval = geval
    data.teval = teval
    data.gbar = gbar
    data.gteval = gteval
    data.num_gteval = num_gteval
}

// -----------------------------------------------------------------------------
// didhetero_duplicate_gteval_pairs()
// Return a space-delimited list of unique duplicate (g,t) pairs supplied by
// the user, preserving the order in which duplicate values are first repeated.
//
// Args:
//   gteval_user - K x 2 matrix of explicit user pairs
//
// Returns:
//   Empty string if there are no duplicates, otherwise "(g1,t1) (g2,t2) ..."
// -----------------------------------------------------------------------------
string scalar didhetero_duplicate_gteval_pairs(real matrix gteval_user)
{
    real scalar i, j, already_listed
    real matrix duplicate_pairs
    string scalar duplicate_labels, pair_label

    duplicate_pairs = J(0, 2, .)
    duplicate_labels = ""

    for (i = 1; i <= rows(gteval_user); i++) {
        for (j = 1; j < i; j++) {
            if (gteval_user[i, 1] == gteval_user[j, 1] &
                gteval_user[i, 2] == gteval_user[j, 2]) {
                already_listed = 0
                if (rows(duplicate_pairs) > 0) {
                    already_listed =
                        sum((duplicate_pairs[., 1] :== gteval_user[i, 1]) :&
                            (duplicate_pairs[., 2] :== gteval_user[i, 2])) > 0
                }
                if (!already_listed) {
                    duplicate_pairs = duplicate_pairs \ gteval_user[i, .]
                    pair_label = "(" +
                        strofreal(gteval_user[i, 1], "%9.0g") + "," +
                        strofreal(gteval_user[i, 2], "%9.0g") + ")"
                    if (strlen(duplicate_labels) > 0) {
                        duplicate_labels = duplicate_labels + " "
                    }
                    duplicate_labels = duplicate_labels + pair_label
                }
                break
            }
        }
    }

    return(duplicate_labels)
}

// -----------------------------------------------------------------------------
// didhetero_user_gteval_in_domain()
// Check whether an explicit (g,t) pair belongs to the same admissible domain
// used by the automatic gteval builder for the current sample and options.
//
// Args:
//   data          - DidHeteroData struct
//   g1            - treatment group value
//   t1            - evaluation time
//   anticipation  - anticipation parameter delta >= 0
//   control_group - "nevertreated" or "notyettreated"
//   pretrend      - 0 or 1
//
// Returns:
//   1 if the pair is admissible, 0 otherwise
// -----------------------------------------------------------------------------
real scalar didhetero_user_gteval_in_domain(
    struct DidHeteroData scalar data,
    real scalar g1,
    real scalar t1,
    real scalar anticipation,
    string scalar control_group,
    real scalar pretrend)
{
    real colvector supp_t
    real scalar gbar, has_never

    supp_t = sort(uniqrows(data.t_vals), 1)
    has_never = sum(data.G :== 0) > 0
    gbar = has_never ? . : max(sort(uniqrows(data.G), 1))

    return(didhetero_pair_in_domain(
        g1, t1, anticipation, control_group,
        pretrend, has_never, gbar, supp_t))
}

// -----------------------------------------------------------------------------
// didhetero_validate_user_gteval()
// Reject user-specified gteval() pairs that are outside the estimable domain
// implied by the sample support, control group, anticipation, and pretrend.
//
// Args:
//   data          - DidHeteroData struct
//   gteval_user   - K x 2 matrix of explicit user pairs
//   anticipation  - anticipation parameter delta >= 0
//   control_group - "nevertreated" or "notyettreated"
//   pretrend      - 0 or 1
//
// Side effects:
//   Raises r(198) with a user-facing message when any pair is invalid.
// -----------------------------------------------------------------------------
void didhetero_validate_user_gteval(
    struct DidHeteroData scalar data,
    real matrix gteval_user,
    real scalar anticipation,
    string scalar control_group,
    real scalar pretrend)
{
    real scalar i, g1, t1
    string scalar duplicate_pairs

    duplicate_pairs = didhetero_duplicate_gteval_pairs(gteval_user)
    if (strlen(duplicate_pairs) > 0) {
        _error(198,
            "gteval() contains duplicate (g,t) pairs: " + duplicate_pairs)
    }

    for (i = 1; i <= rows(gteval_user); i++) {
        g1 = gteval_user[i, 1]
        t1 = gteval_user[i, 2]

        if (!didhetero_user_gteval_in_domain(
                data, g1, t1, anticipation, control_group, pretrend)) {
            _error(198,
                "gteval() contains (g,t) outside the identification domain for control_group(" +
                control_group + "): (" + strofreal(g1) + ", " + strofreal(t1) + ")")
        }
    }
}

end
