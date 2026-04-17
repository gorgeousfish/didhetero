mata:

// =============================================================================
// didhetero_utils.mata
// Utility functions for didhetero-stata
//
// Functions:
//   1. didhetero_seq()            - Equally-spaced sequence generation
//   2. didhetero_trapz()          - Trapezoidal rule numerical integration
//   3. didhetero_quantile()       - Quantile calculation (Hyndman-Fan type 7)
//   4. didhetero_mammen_weights() - Mammen wild bootstrap weights
//   5. didhetero_analytical_cv()  - Analytical critical value
//   6. didhetero_parse_xformula_locals() - Parse xformula() into resolved terms
//   7. didhetero_unique_tokens()  - Remove duplicate xformula tokens, preserving order
//   8. _didhetero_intersect()     - Set intersection (sorted)
//   9. _didhetero_setdiff()       - Set difference (sorted)
//  10. didhetero_period_ord()     - Ordinal position of a time label in t_vals
//  11. didhetero_period_at()      - Time label at a given ordinal position
//  12. didhetero_group_is_valid() - Check whether a group has a valid base period
//  13. didhetero_time_is_valid()  - Check whether a time label is admissible
//  14. didhetero_pair_in_domain() - Check whether a (g,t) pair is estimable
//  15. didhetero_build_gteval()   - Build (g,t) evaluation point matrices
//  16. didhetero_duplicate_gteval_pairs() - Format duplicate user (g,t) pairs
//  17. didhetero_user_gteval_in_domain() - Check whether a user pair is estimable
//  18. didhetero_validate_user_gteval() - Reject explicit gteval() outside domain
//  19. didhetero_selectindex()    - Index selection (Stata <17 compat)
//  20. didhetero_gen_z_supp()     - Z support grid generation (100 points)
//  21. didhetero_init_arrays()    - Core estimation array initialization
//  22. didhetero_init_core_arrays() - Core estimation array initialization (delayed)
//  23. didhetero_prepare_data()   - Main data preparation (long->wide, struct)
//
// References:
//   Paper: Imai, Qin, Yanagi (2025)
//     - Section 4.2.4 Eq. 21-22 (analytical critical value)
//     - Section 4.2.4 (Mammen wild bootstrap)
//   Sections 4.2.1–4.2.4 (preparation, integration, critical values, bootstrap)
// =============================================================================

// -----------------------------------------------------------------------------
// didhetero_selectindex()
// Return indices where v is nonzero (compatible with Stata < 17)
//
// Args:
//   v - real colvector of values
//
// Returns:
//   real rowvector of 1-based indices where v != 0
//   Returns J(1, 0, .) (empty rowvector) if no nonzero elements
//
// Note: Replaces Stata 17+ selectindex() for backward compatibility
// -----------------------------------------------------------------------------
real rowvector didhetero_selectindex(real colvector v)
{
    real rowvector idx
    real scalar i, k, n
    
    n = rows(v)
    k = 0
    idx = J(1, n, .)
    
    for (i = 1; i <= n; i++) {
        if (v[i] != 0) {
            k++
            idx[k] = i
        }
    }
    
    if (k == 0) return(J(1, 0, .))
    return(idx[1..k])
}

// -----------------------------------------------------------------------------
// didhetero_trim()
// Trim leading and trailing whitespace.
// -----------------------------------------------------------------------------
string scalar didhetero_trim(string scalar s)
{
    return(strtrim(s))
}

// -----------------------------------------------------------------------------
// didhetero_has_top_level_char()
// Check whether a character appears outside parentheses.
// -----------------------------------------------------------------------------
real scalar didhetero_has_top_level_char(string scalar s, string scalar needle)
{
    real scalar i, depth
    string scalar ch

    depth = 0
    for (i = 1; i <= strlen(s); i++) {
        ch = substr(s, i, 1)
        if (ch == "(") depth++
        else if (ch == ")") {
            depth--
            if (depth < 0) {
                _error(198, "xformula() has mismatched parentheses")
            }
        }

        if (ch == needle & depth == 0) {
            return(1)
        }
    }

    if (depth != 0) {
        _error(198, "xformula() has mismatched parentheses")
    }

    return(0)
}

// -----------------------------------------------------------------------------
// didhetero_split_top()
// Split a string on a single-character delimiter outside parentheses.
// -----------------------------------------------------------------------------
string rowvector didhetero_split_top(string scalar s, string scalar delim)
{
    string rowvector parts
    string scalar token, ch
    real scalar i, depth

    parts = J(1, 0, "")
    token = ""
    depth = 0

    for (i = 1; i <= strlen(s); i++) {
        ch = substr(s, i, 1)

        if (ch == "(") depth++
        else if (ch == ")") {
            depth--
            if (depth < 0) {
                _error(198, "xformula() has mismatched parentheses")
            }
        }

        if (ch == delim & depth == 0) {
            parts = parts, didhetero_trim(token)
            token = ""
        }
        else {
            token = token + ch
        }
    }

    if (depth != 0) {
        _error(198, "xformula() has mismatched parentheses")
    }

    parts = parts, didhetero_trim(token)
    return(parts)
}

// -----------------------------------------------------------------------------
// didhetero_is_identifier()
// Check whether a token is a simple Stata-style variable name.
// -----------------------------------------------------------------------------
real scalar didhetero_is_identifier(string scalar s)
{
    real scalar i, code
    string scalar ch

    s = didhetero_trim(s)
    if (strlen(s) == 0) return(0)

    ch = substr(s, 1, 1)
    code = ascii(ch)
    if (!((code >= 65 & code <= 90) | (code >= 97 & code <= 122) | code == 95)) {
        return(0)
    }

    for (i = 2; i <= strlen(s); i++) {
        ch = substr(s, i, 1)
        code = ascii(ch)
        if (!((code >= 65 & code <= 90) | (code >= 97 & code <= 122) |
              (code >= 48 & code <= 57) | code == 95)) {
            return(0)
        }
    }

    return(1)
}

// -----------------------------------------------------------------------------
// didhetero_is_wrapped_I()
// Check whether a term is of the form I(expr).
// -----------------------------------------------------------------------------
real scalar didhetero_is_wrapped_I(string scalar term)
{
    real scalar i, depth
    string scalar ch

    term = didhetero_trim(term)
    if (strlen(term) < 4) return(0)
    if (substr(term, 1, 2) != "I(") return(0)
    if (substr(term, strlen(term), 1) != ")") return(0)

    depth = 0
    for (i = 2; i <= strlen(term); i++) {
        ch = substr(term, i, 1)
        if (ch == "(") depth++
        else if (ch == ")") {
            depth--
            if (depth == 0 & i < strlen(term)) {
                return(0)
            }
            if (depth < 0) {
                return(0)
            }
        }
    }

    return(depth == 0)
}

// -----------------------------------------------------------------------------
// didhetero_join_combo()
// Join selected factors with ":" to represent an interaction term.
// -----------------------------------------------------------------------------
string scalar didhetero_join_combo(string rowvector factors, real rowvector idx)
{
    real scalar j
    string scalar combo

    combo = didhetero_trim(factors[idx[1]])
    for (j = 2; j <= cols(idx); j++) {
        combo = combo + ":" + didhetero_trim(factors[idx[j]])
    }

    return(combo)
}

// -----------------------------------------------------------------------------
// didhetero_expand_star_term()
// Expand A*B*C into main effects plus interactions.
// -----------------------------------------------------------------------------
string rowvector didhetero_expand_star_term(string scalar term)
{
    string rowvector factors, expanded
    real rowvector idx
    real scalar i, j, k, n_factors

    factors = didhetero_split_top(term, "*")
    n_factors = cols(factors)
    if (n_factors == 1) {
        return(didhetero_trim(term))
    }

    expanded = J(1, 0, "")
    for (k = 1; k <= n_factors; k++) {
        idx = J(1, k, .)
        for (i = 1; i <= k; i++) {
            idx[i] = i
        }

        while (1) {
            expanded = expanded, didhetero_join_combo(factors, idx)

            for (i = k; i >= 1; i--) {
                if (idx[i] < n_factors - k + i) {
                    idx[i] = idx[i] + 1
                    for (j = i + 1; j <= k; j++) {
                        idx[j] = idx[j - 1] + 1
                    }
                    break
                }
            }

            if (i == 0) break
        }
    }

    return(expanded)
}

// -----------------------------------------------------------------------------
// didhetero_formula_term_expr()
// Convert a supported formula term into a numeric Stata expression.
// -----------------------------------------------------------------------------
string scalar didhetero_formula_term_expr(string scalar term)
{
    string rowvector parts
    string scalar expr
    real scalar i

    term = didhetero_trim(term)
    if (didhetero_is_wrapped_I(term)) {
        return(didhetero_trim(substr(term, 3, strlen(term) - 3)))
    }

    parts = didhetero_split_top(term, ":")
    if (cols(parts) > 1) {
        expr = "(" + didhetero_formula_term_expr(parts[1]) + ")"
        for (i = 2; i <= cols(parts); i++) {
            expr = expr + "*(" + didhetero_formula_term_expr(parts[i]) + ")"
        }
        return(expr)
    }

    if (didhetero_is_identifier(term)) {
        return(term)
    }

    _error(198,
        "unsupported xformula() term: " + term +
        " (supported: bare variables, I(expr), :, and *)")
    return("")
}

// -----------------------------------------------------------------------------
// didhetero_parse_xformula_locals()
// Parse xformula() into additive terms and expose them through Stata locals.
// -----------------------------------------------------------------------------
void didhetero_parse_xformula_locals()
{
    string scalar raw, body, term
    string rowvector additive_terms, expanded_terms, final_terms, minus_terms
    string rowvector kinds, payloads
    real scalar i, j, seen, legacy_mode, has_intercept

    raw = didhetero_trim(st_local("xformula"))
    st_local("dh_xformula_display", raw)
    st_local("dh_xformula_n", "0")
    st_local("dh_xformula_has_intercept", "1")
    has_intercept = 1

    if (raw == "") {
        return
    }

    legacy_mode = (substr(raw, 1, 1) != "~" &
        strpos(raw, "+") == 0 &
        strpos(raw, ":") == 0 &
        strpos(raw, "*") == 0 &
        strpos(raw, "I(") == 0 &
        cols(tokens(raw)) > 1)

    body = raw
    if (substr(body, 1, 1) == "~") {
        body = didhetero_trim(substr(body, 2, strlen(body) - 1))
    }

    if (!legacy_mode & didhetero_has_top_level_char(body, "-")) {
        minus_terms = didhetero_split_top(body, "-")
        body = didhetero_trim(minus_terms[1])
        for (i = 2; i <= cols(minus_terms); i++) {
            term = didhetero_trim(minus_terms[i])
            if (term == "1") {
                has_intercept = 0
            }
            else {
                _error(198,
                    "xformula() does not support term removal; only -1 is allowed to suppress the intercept")
            }
        }
    }

    if (legacy_mode) {
        additive_terms = tokens(raw)
    }
    else {
        additive_terms = didhetero_split_top(body, "+")
    }

    final_terms = J(1, 0, "")
    for (i = 1; i <= cols(additive_terms); i++) {
        term = didhetero_trim(additive_terms[i])
        if (term == "" | term == "1") continue
        if (term == "0" | term == "-1") {
            has_intercept = 0
            continue
        }

        expanded_terms = didhetero_expand_star_term(term)
        for (j = 1; j <= cols(expanded_terms); j++) {
            term = didhetero_trim(expanded_terms[j])
            if (term == "" | term == "1") continue

            seen = 0
            if (cols(final_terms) > 0) {
                seen = sum(final_terms :== term) > 0
            }
            if (!seen) {
                final_terms = final_terms, term
            }
        }
    }

    kinds = J(1, 0, "")
    payloads = J(1, 0, "")
    for (i = 1; i <= cols(final_terms); i++) {
        term = final_terms[i]
        if (didhetero_is_identifier(term)) {
            kinds = kinds, "var"
            payloads = payloads, term
        }
        else {
            kinds = kinds, "expr"
            payloads = payloads, didhetero_formula_term_expr(term)
        }
    }

    st_local("dh_xformula_n", strofreal(cols(final_terms), "%9.0g"))
    st_local("dh_xformula_has_intercept", strofreal(has_intercept, "%9.0g"))
    for (i = 1; i <= cols(final_terms); i++) {
        st_local(sprintf("dh_xformula_term%g", i), final_terms[i])
        st_local(sprintf("dh_xformula_kind%g", i), kinds[i])
        st_local(sprintf("dh_xformula_payload%g", i), payloads[i])
    }
}

// -----------------------------------------------------------------------------
// didhetero_unique_tokens()
// Remove duplicate space-delimited tokens while preserving first appearance.
// Repeated terms in a formula should not create extra design-matrix columns.
//
// Args:
//   xformula - raw xformula string from Stata option parsing
//
// Returns:
//   string rowvector of unique tokens in original order
// -----------------------------------------------------------------------------
string rowvector didhetero_unique_tokens(string scalar xformula)
{
    string rowvector raw_tokens, uniq_tokens
    real scalar i, j, n_raw, n_uniq, seen

    raw_tokens = tokens(xformula)
    n_raw = cols(raw_tokens)
    if (n_raw == 0) return(J(1, 0, ""))

    uniq_tokens = J(1, n_raw, "")
    n_uniq = 0

    for (i = 1; i <= n_raw; i++) {
        seen = 0
        for (j = 1; j <= n_uniq; j++) {
            if (uniq_tokens[j] == raw_tokens[i]) {
                seen = 1
                break
            }
        }
        if (!seen) {
            n_uniq++
            uniq_tokens[n_uniq] = raw_tokens[i]
        }
    }

    return(uniq_tokens[|1, 1 \ 1, n_uniq|])
}

// -----------------------------------------------------------------------------
// didhetero_seq()
// Generate equispaced sequence from `from` to `to` with `length` points
//
// Args:
//   from   - start value
//   to     - end value
//   length - number of points (must be >= 1)
//
// Returns:
//   real colvector of equispaced values
//
// Note: Uses Mata built-in rangen() for numerical consistency
// -----------------------------------------------------------------------------
real colvector didhetero_seq(real scalar from, real scalar to, real scalar length)
{
    // Input validation
    if (length < 1) _error("length must be at least 1")
    
    // Special case: single point
    if (length == 1) return(from)
    
    // Use Mata built-in rangen for numerical consistency
    return(rangen(from, to, length))
}

// -----------------------------------------------------------------------------
// didhetero_trapz()
// Trapezoidal rule numerical integration
//
// Args:
//   x - real colvector of grid points (abscissae)
//   y - real colvector of function values at grid points
//
// Returns:
//   real scalar approximation of integral of y over x
//
// Paper ref: Numerical integration for IMSE bandwidth selection
// -----------------------------------------------------------------------------
real scalar didhetero_trapz(real colvector x, real colvector y)
{
    real scalar n, integral, i
    
    // Input validation: x and y must have same length
    n = length(x)
    if (n != length(y)) _error("x and y must have same length")
    
    // Edge case: fewer than 2 points, no interval to integrate
    if (n < 2) return(0)
    
    // Trapezoidal rule: sum of (x[i+1]-x[i]) * (y[i]+y[i+1]) / 2
    integral = 0
    for (i = 1; i < n; i++) {
        integral = integral + (x[i+1] - x[i]) * (y[i] + y[i+1]) / 2
    }
    return(integral)
}

// -----------------------------------------------------------------------------
// didhetero_quantile()
// Quantile function using Hyndman & Fan (1996) type=7 algorithm
// (Hyndman & Fan type 7, the standard default in statistical software)
//
// Args:
//   x     - real colvector of data values (must be non-empty)
//   probs - real vector of probabilities in [0, 1]
//
// Returns:
//   real vector of quantile values corresponding to probs
//
// Note: Type 7 continuous index: h = (n-1)*p + 1, with linear interpolation
// -----------------------------------------------------------------------------
real vector didhetero_quantile(real colvector x, real vector probs)
{
    real colvector sorted_x
    real vector result
    real scalar n, i, h_idx, lo, hi
    
    // Input validation: x must be non-empty
    if (length(x) < 1) _error("x must be non-empty")
    
    // Input validation: probs must be in [0, 1]
    if (min(probs) < 0 | max(probs) > 1) _error("probs must be between 0 and 1")
    
    // Sort data ascending
    sorted_x = sort(x, 1)
    n = length(sorted_x)
    result = J(length(probs), 1, .)
    
    for (i = 1; i <= length(probs); i++) {
        // Linear interpolation quantile: continuous index h = (n-1)*p + 1
        h_idx = (n - 1) * probs[i] + 1
        lo = floor(h_idx)
        hi = ceil(h_idx)
        
        // Boundary protection
        if (lo < 1) lo = 1
        if (hi > n) hi = n
        
        if (lo == hi) {
            result[i] = sorted_x[lo]
        }
        else {
            // Linear interpolation
            result[i] = sorted_x[lo] + (h_idx - lo) * (sorted_x[hi] - sorted_x[lo])
        }
    }
    return(result)
}

// -----------------------------------------------------------------------------
// didhetero_mammen_weights()
// Generate Mammen (1993) two-point wild bootstrap weights
//
// Args:
//   n - number of weights to generate (must be >= 1)
//
// Returns:
//   real colvector of bootstrap weights
//
// Note: Paper Section 4.2.4: kappa = (sqrt(5)+1)/2 (golden ratio)
//       V* = 2-kappa w.p. kappa/sqrt(5), 1+kappa w.p. 1-kappa/sqrt(5)
//       Moment properties: E[V]=1, Var[V]=1, E[(V-1)^3]=1
//       Paper ref: Section 4.2.4 Eq. 23
// -----------------------------------------------------------------------------
real colvector didhetero_mammen_weights(real scalar n)
{
    real scalar kappa, p_low
    real colvector u, weights
    
    // Input validation
    if (n < 1) _error("n must be a positive integer")
    
    // Golden ratio
    // Paper Section 4.2.4: kappa = (sqrt(5) + 1) / 2
    kappa = (sqrt(5) + 1) / 2
    p_low = kappa / sqrt(5)
    
    // Two-point distribution via uniform threshold
    // V* = 2-kappa with prob kappa/sqrt(5), 1+kappa with prob 1-kappa/sqrt(5)
    u = runiform(n, 1)
    weights = (u :<= p_low) :* (2 - kappa) + (u :> p_low) :* (1 + kappa)
    
    return(weights)
}

// -----------------------------------------------------------------------------
// didhetero_analytical_cv()
// Analytical critical value from extreme value distribution
//
// Args:
//   a      - lower bound of the interval
//   b      - upper bound of the interval (must be > a)
//   h      - bandwidth (must be > 0)
//   lambda - kernel-dependent constant
//   alpha  - significance level in (0, 1)
//
// Returns:
//   real scalar critical value c_hat
//
// Note: Paper Eq. 21-22
//       a_n^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
//       c_hat = sqrt(a_n^2 - 2*log(log(1/sqrt(1-alpha))))
//       Paper ref: Section 4.2.4 Eq. 21-22
// -----------------------------------------------------------------------------
real scalar didhetero_analytical_cv(real scalar a, real scalar b,
                                     real scalar h, real scalar lambda,
                                     real scalar alpha)
{
    real scalar a_n_sq, c_hat, inner_log
    
    // Input validation
    if (h <= 0) _error("h must be positive")
    if (b <= a) _error("b must be greater than a")
    if (alpha <= 0 | alpha >= 1) _error("alpha must be between 0 and 1 (exclusive)")
    
    // Paper Eq. 21: a_n^2 = 2*log((b-a)/h) + 2*log(sqrt(lambda)/(2*pi))
    a_n_sq = 2 * log((b - a) / h) + 2 * log(sqrt(lambda) / (2 * c("pi")))
    
    // Paper Eq. 22: inner_log = log(log(1/sqrt(1-alpha)))
    inner_log = log(log(1 / sqrt(1 - alpha)))
    
    // Boundary protection: fallback to normal quantile when a_n^2 is too small
    if (a_n_sq - 2 * inner_log <= 0) {
        c_hat = invnormal(1 - alpha/2)
    }
    else {
        // Paper Eq. 22: c_hat = sqrt(a_n_sq - 2*inner_log)
        c_hat = sqrt(a_n_sq - 2 * inner_log)
    }
    
    return(c_hat)
}

// -----------------------------------------------------------------------------
// _didhetero_intersect()
// Return sorted intersection of two vectors
//
// Args:
//   a - real colvector
//   b - real colvector
//
// Returns:
//   real colvector of elements present in both a and b, sorted ascending
//   Returns J(0, 1, .) if intersection is empty
// -----------------------------------------------------------------------------
real colvector _didhetero_intersect(real colvector a, real colvector b)
{
    real colvector result
    real scalar i, n_a
    
    // Handle empty inputs
    n_a = rows(a)
    if (n_a == 0 | rows(b) == 0) return(J(0, 1, .))
    
    result = J(0, 1, .)
    for (i = 1; i <= n_a; i++) {
        if (any(b :== a[i])) {
            result = result \ a[i]
        }
    }
    
    if (rows(result) == 0) return(J(0, 1, .))
    return(sort(result, 1))
}

// -----------------------------------------------------------------------------
// _didhetero_setdiff()
// Return sorted set difference: elements in a but not in b
//
// Args:
//   a - real colvector
//   b - real colvector
//
// Returns:
//   real colvector of elements in a but not in b, sorted ascending
//   Returns J(0, 1, .) if result is empty
// -----------------------------------------------------------------------------
real colvector _didhetero_setdiff(real colvector a, real colvector b)
{
    real colvector result
    real scalar i, n_a
    
    // Handle empty inputs
    n_a = rows(a)
    if (n_a == 0) return(J(0, 1, .))
    if (rows(b) == 0) return(sort(a, 1))
    
    result = J(0, 1, .)
    for (i = 1; i <= n_a; i++) {
        if (!any(b :== a[i])) {
            result = result \ a[i]
        }
    }
    
    if (rows(result) == 0) return(J(0, 1, .))
    return(sort(result, 1))
}

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
        return(t_ord != g_ord - anticipation - 1)
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
//
// Paper ref: Section 2 (Setup), Section 4.2.1 (Estimation procedure)
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

    // Guard the paper's support restriction even if callers bypass ado-level
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
    
    // === Step 4: gbar computation (Paper Section 2) ===
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

// -----------------------------------------------------------------------------
// didhetero_gen_z_supp()
// Generate uniform support grid on [min(Z), max(Z)] with 100 points
//
// Args:
//   Z - real colvector of continuous treatment values
//
// Returns:
//   real colvector of 100 equally-spaced points from min(Z) to max(Z)
//
// Paper ref: Section 4.2.5 (Bandwidth selection)
// -----------------------------------------------------------------------------
real colvector didhetero_gen_z_supp(real colvector Z)
{
    return(rangen(min(Z), max(Z), 100))
}

// -----------------------------------------------------------------------------
// didhetero_init_arrays()
// Initialize the 6 core estimation arrays in DidHeteroData.
// Called after num_gteval is determined.
//
// Args:
//   d          - DidHeteroData struct (modified in-place)
//   n          - number of cross-sectional units
//   num_zeval  - number of Z evaluation points
//   num_gteval - number of (g,t) evaluation pairs
//
// Side effects:
//   Populates d.A_g_t, d.B_g_t, d.G_g, d.mu_G_g, d.mu_E_g_t, d.mu_F_g_t
//
// Paper ref: Section 4.2.1 (Three-step estimation procedure)
//       A_g_t, B_g_t are pointer arrays (1 x num_gteval), each n x num_zeval
//       G_g is n x num_gteval
//       mu_G_g, mu_E_g_t, mu_F_g_t are num_zeval x num_gteval
// -----------------------------------------------------------------------------
void didhetero_init_arrays(struct DidHeteroData scalar d,
                           real scalar n,
                           real scalar num_zeval,
                           real scalar num_gteval)
{
    real scalar k

    // A_g_t, B_g_t: pointer rowvector, each element is n x num_zeval matrix of missing
    d.A_g_t = J(1, num_gteval, NULL)
    d.B_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.A_g_t[k] = &(J(n, num_zeval, .))
        d.B_g_t[k] = &(J(n, num_zeval, .))
    }

    // G_g: n x num_gteval group indicator matrix
    d.G_g = J(n, num_gteval, .)

    // mu_G_g: num_zeval x num_gteval conditional group density
    d.mu_G_g = J(num_zeval, num_gteval, .)

    // mu_E_g_t, mu_F_g_t: num_zeval x num_gteval
    d.mu_E_g_t = J(num_zeval, num_gteval, .)
    d.mu_F_g_t = J(num_zeval, num_gteval, .)
}

// -----------------------------------------------------------------------------
// didhetero_get_uniformall()
// Return the effective uniformall flag stored in the global data struct.
// This is a test helper to verify propagation of the single (g,t) override.
// -----------------------------------------------------------------------------
real scalar didhetero_get_uniformall()
{
    external struct DidHeteroData scalar _dh_data
    return(_dh_data.uniformall)
}

// -----------------------------------------------------------------------------
// _didhetero_read_panel_id()
// Read a Stata panel id variable into a numeric column vector.
//
// Numeric ids are read directly. String ids are converted into stable
// consecutive numeric codes using the current row order. The caller guarantees
// the data are pre-sorted by (id, time), so equal ids are contiguous and this
// scan preserves panel blocks exactly.
// -----------------------------------------------------------------------------
real colvector _didhetero_read_panel_id(string scalar idvar)
{
    real scalar n_obs, i, id_code
    real colvector id_long
    string matrix id_str

    if (!st_isstrvar(idvar)) {
        return(st_data(., idvar))
    }

    id_str = st_sdata(., idvar)
    n_obs = rows(id_str)
    id_long = J(n_obs, 1, .)

    if (n_obs == 0) {
        return(id_long)
    }

    id_code = 1
    id_long[1] = id_code

    for (i = 2; i <= n_obs; i++) {
        if (id_str[i, 1] != id_str[i - 1, 1]) {
            id_code++
        }
        id_long[i] = id_code
    }

    return(id_long)
}

// -----------------------------------------------------------------------------
// didhetero_prepare_data()
// Main data preparation: read Stata data, reshape, and fill DidHeteroData struct
//
// Args:
//   depvar   - string scalar, name of the dependent variable
//   idvar    - string scalar, name of the panel id variable
//   timevar  - string scalar, name of the time variable
//   groupvar - string scalar, name of the group variable
//   zvar     - string scalar, name of the continuous treatment variable
//   zeval    - real colvector, evaluation points for Z
//   xformula - string scalar, space-separated covariate names (empty = none)
//
// Returns:
//   struct DidHeteroData scalar with all fields populated
//
// Paper ref: Section 2 (Setup), Section 4.2.1 (Estimation procedure)
//       Data must be pre-sorted by (id, time) via _didhetero_validate.ado
// -----------------------------------------------------------------------------
struct DidHeteroData scalar didhetero_prepare_data(
    string scalar depvar,
    string scalar idvar,
    string scalar timevar,
    string scalar groupvar,
    string scalar zvar,
    real colvector zeval,
    string scalar xformula)
{
    struct DidHeteroData scalar d
    real colvector Y_long, id_long, time_long, G_long, Z_long
    real colvector id_unique, t_vals, G, Z, id, Z_supp
    real colvector period1_idx, x_long
    real matrix Y_wide, X, X_sub
    real scalar n, T_num, period1, num_zeval, i, row_start, row_end, j, has_intercept
    string rowvector tok
    string scalar xformula_has_intercept

    // === Step 1: Read data from Stata ===
    // Data already sorted by _didhetero_validate.ado
    Y_long    = st_data(., depvar)
    id_long   = _didhetero_read_panel_id(idvar)
    time_long = st_data(., timevar)
    G_long    = st_data(., groupvar)
    Z_long    = st_data(., zvar)

    // === Step 2: Basic dimensions ===
    // Number of cross-sectional units
    id_unique = uniqrows(id_long)
    n = rows(id_unique)

    // Sorted unique time periods
    t_vals = sort(uniqrows(time_long), 1)
    T_num  = rows(t_vals)

    // First time period
    period1 = t_vals[1]

    // === Step 3: Long -> Wide Y_wide ===
    // Data is sorted by (id, time), so each block of T_num rows = one unit
    Y_wide = J(n, T_num, .)
    for (i = 1; i <= n; i++) {
        row_start = (i - 1) * T_num + 1
        row_end   = i * T_num
        Y_wide[i, .] = Y_long[row_start..row_end]'
    }

    // === Step 4: G vector (n x 1) ===
    // Extract one value per id (time-invariant, take first period)
    // Indices: 1, 1+T_num, 1+2*T_num, ..., 1+(n-1)*T_num
    // Note: range(1, n*T_num, T_num) fails when (n*T_num-1) is not divisible
    //       by T_num, because Mata's range() rounds up the point count.
    //       Use explicit vector arithmetic instead.
    period1_idx = (0::(n-1)) :* T_num :+ 1
    G = G_long[period1_idx]

    // === Step 5: Z vector (n x 1) ===
    Z = Z_long[period1_idx]

    // === Step 6: id vector (n x 1) ===
    id = id_unique

    // === Step 7: zeval sort and range validation ===
    zeval = sort(zeval, 1)
    num_zeval = rows(zeval)
    if (min(zeval) < min(Z) | max(zeval) > max(Z)) {
        _error("zeval values must be within the range of Z: [" +
               strofreal(min(Z)) + ", " + strofreal(max(Z)) + "]")
    }

    // === Step 8: Z_supp generation ===
    Z_supp = didhetero_gen_z_supp(Z)

    // === Step 9: Covariate matrix X ===
    // The default design matrix is [intercept, formula_vars].
    // Standard formula syntax can suppress the intercept via
    // `~ 0 + ...` or `~ ... - 1`.
    xformula_has_intercept = st_local("xformula_has_intercept")
    if (xformula_has_intercept == "") has_intercept = 1
    else {
        has_intercept = strtoreal(xformula_has_intercept)
        if (missing(has_intercept)) has_intercept = 1
    }

    if (xformula == "") {
        // No extra covariates: [intercept, Z] -> n x 2
        X = J(n, 1, 1), Z
    }
    else {
        // Collapse duplicate formula terms into unique columns before
        // extracting first-period covariates.
        tok = didhetero_unique_tokens(xformula)
        X_sub = J(n, cols(tok), .)
        for (j = 1; j <= cols(tok); j++) {
            x_long = st_data(., tok[j])
            X_sub[., j] = x_long[period1_idx]
        }
        // Do NOT prepend Z separately — it is included via xformula if needed.
        // Explicit no-intercept formulas keep only the requested covariate
        // columns.
        if (has_intercept) X = J(n, 1, 1), X_sub
        else X = X_sub
    }

    // === Step 10: Fill struct ===
    d.n         = n
    d.T_num     = T_num
    d.Y_wide    = Y_wide
    d.G         = G
    d.Z         = Z
    d.id        = id
    d.t_vals    = t_vals
    d.period1   = period1
    d.zeval     = zeval
    d.num_zeval = num_zeval
    d.Z_supp    = Z_supp
    d.X         = X

    return(d)
}

// -----------------------------------------------------------------------------
// didhetero_init_core_arrays()
// Initialize core estimation arrays after num_gteval is known.
// Called after gteval computation.
//
// Args:
//   d          - DidHeteroData struct (modified in place via pointer)
//   num_gteval - number of (g,t) evaluation pairs
//
// Side effects:
//   Populates d.A_g_t, d.B_g_t, d.G_g, d.mu_G_g, d.mu_E_g_t, d.mu_F_g_t
//
// Paper ref: Section 4.2.1 (Three-step estimation procedure)
// -----------------------------------------------------------------------------
void didhetero_init_core_arrays(struct DidHeteroData scalar d,
                                real scalar num_gteval)
{
    real scalar k

    // A_g_t: 1 x num_gteval pointer array, each -> n x num_zeval DR score matrix
    d.A_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.A_g_t[k] = &(J(d.n, d.num_zeval, .))
    }

    // B_g_t: 1 x num_gteval pointer array, each -> n x num_zeval influence function
    d.B_g_t = J(1, num_gteval, NULL)
    for (k = 1; k <= num_gteval; k++) {
        d.B_g_t[k] = &(J(d.n, d.num_zeval, .))
    }

    // G_g: n x num_gteval group indicator matrix
    d.G_g = J(d.n, num_gteval, .)

    // mu_G_g: num_zeval x num_gteval conditional group density
    d.mu_G_g = J(d.num_zeval, num_gteval, .)

    // mu_E_g_t, mu_F_g_t: num_zeval x num_gteval
    d.mu_E_g_t = J(d.num_zeval, num_gteval, .)
    d.mu_F_g_t = J(d.num_zeval, num_gteval, .)
}

// -----------------------------------------------------------------------------
// didhetero_init_from_ado()
// Wrapper called from catt_gt.ado to initialize data and kernel structs.
// Reads Stata locals set by catt_gt.ado, calls prepare_data and kernel_consts,
// and stores results as Mata external globals.
//
// This function exists because Stata's inline mata { } blocks inside ado files
// do NOT support the `external` keyword. By compiling this function in a .mata
// file, we can properly declare external structs.
//
// Side effects:
//   Creates external globals _dh_data (DidHeteroData) and _dh_kc (DidHeteroKernelConsts)
// -----------------------------------------------------------------------------
void didhetero_init_from_ado()
{
    external struct DidHeteroData scalar _dh_data
    external struct DidHeteroKernelConsts scalar _dh_kc
    
    string rowvector _zeval_tokens
    real rowvector _zeval_vec
    real colvector _zeval_col
    real scalar _i, _porder
    
    // Parse zeval from Stata local into Mata vector
    _zeval_tokens = tokens(st_local("zeval"))
    _zeval_vec = J(1, cols(_zeval_tokens), .)
    for (_i = 1; _i <= cols(_zeval_tokens); _i++) {
        _zeval_vec[_i] = strtoreal(_zeval_tokens[_i])
    }
    _zeval_col = _zeval_vec'
    
    // Call data preparation function
    _dh_data = didhetero_prepare_data(  
        st_local("depvar"),             
        st_local("id"),                 
        st_local("time"),               
        st_local("group"),              
        st_local("z"),                  
        _zeval_col,                     
        st_local("xformula")            
    )
    
    // Store configuration parameters into struct
    _dh_data.control_group = st_local("control")
    _dh_data.anticipation  = strtoreal(st_local("anticip"))
    _dh_data.porder        = strtoreal(st_local("porder"))
    _dh_data.kernel        = st_local("kernel")
    _dh_data.alp           = strtoreal(st_local("alp"))
    _dh_data.biters        = strtoreal(st_local("biters"))
    _dh_data.uniformall    = strtoreal(st_local("uniform"))
    
    // Initialize kernel constants
    _dh_kc = didhetero_kernel_consts(st_local("kernel"))
    
    // Copy key derived constants to DidHeteroData struct
    _dh_data.const_B1 = _dh_kc.const_B1
    _dh_data.const_B2 = _dh_kc.const_B2
    _dh_data.lambda   = _dh_kc.lambda
    
    // Select const_V based on porder (Paper Section 4.2.2)
    _porder = strtoreal(st_local("porder"))
    _dh_data.const_V = (_porder == 1) * _dh_kc.const_V1 + (_porder == 2) * _dh_kc.const_V2
}

// ---------------------------------------------------------------------------
// didhetero_init_param_results()
// Factory function for DidHeteroParamResults struct.
// ---------------------------------------------------------------------------
struct DidHeteroParamResults scalar didhetero_init_param_results()
{
    struct DidHeteroParamResults scalar r
    r.gps_mat  = J(0, 0, .)
    r.gps_coef = J(0, 0, .)
    r.or_mat   = J(0, 0, .)
    r.or_coef  = J(0, 0, .)
    r.ctrl_type = ""
    return(r)
}

// ---------------------------------------------------------------------------
// didhetero_init_stage1_results()
// Factory function for DidHeteroStage1Results struct.
// ---------------------------------------------------------------------------
struct DidHeteroStage1Results scalar didhetero_init_stage1_results()
{
    struct DidHeteroStage1Results scalar r
    r.gps_mat  = J(0, 0, .)
    r.gps_coef = J(0, 0, .)
    r.or_mat   = J(0, 0, .)
    r.or_coef  = J(0, 0, .)
    r.ctrl_type = ""
    r.Z_supp   = J(0, 1, .)
    r.kd0_Z    = J(0, 1, .)
    r.kd1_Z    = J(0, 1, .)
    return(r)
}

end
