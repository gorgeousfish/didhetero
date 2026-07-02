mata:

// =============================================================================
// didhetero_utils_formula.mata
// Formula parsing and index utility functions
//
// Provides:
//   - didhetero_selectindex()       // Return indices where v is nonzero
//   - didhetero_trim()              // Trim whitespace
//   - didhetero_has_top_level_char()// Check char outside parentheses
//   - didhetero_split_top()         // Split at top-level delimiter
//   - didhetero_is_identifier()     // Validate Stata identifier
//   - didhetero_is_wrapped_I()      // Detect I() wrapper
//   - didhetero_join_combo()        // Join factor combination
//   - didhetero_expand_star_term()  // Expand * interactions
//   - didhetero_formula_term_expr() // Convert formula term to expression
//   - didhetero_parse_xformula_locals() // Parse xformula to Stata locals
//   - didhetero_unique_tokens()     // Extract unique variable tokens
//
// Requires:
//   - (none -- pure string processing)
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

end
