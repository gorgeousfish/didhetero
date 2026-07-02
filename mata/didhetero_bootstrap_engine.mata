// =============================================================================
// didhetero_bootstrap_engine.mata
// Bootstrap Infrastructure Engine — Strategy Pattern Base Layer
//
// Provides:
//   - dh_boot_generate_weights()   : weight generation (Mammen / extensible)
//   - dh_boot_quantile()           : bootstrap quantile with missing handling
//   - dh_boot_compute_sup()        : supremum t-statistic extraction
//   - dh_boot_mammen_params()      : Mammen two-point distribution constants
//   - dh_boot_max_nonmissing()     : max ignoring missing values
//
// Requires:
//   - didhetero_types.mata (struct BootPrecomp)
//
// Design:
//   This engine module isolates bootstrap infrastructure (weight generation,
//   quantile computation, sup-t extraction) from application-specific logic
//   (CATT estimation, aggregation). Application modules call these primitives
//   rather than reimplementing bootstrap mechanics.
//
// Naming convention: dh_boot_ prefix for all public engine functions.
// =============================================================================

mata:

// =============================================================================
// Section 1: Mammen Weight Constants
// =============================================================================

// -----------------------------------------------------------------------------
// dh_boot_mammen_params()
//
// Return Mammen (1993) two-point distribution parameters.
// The distribution satisfies E[V*] = 1, Var[V*] = 1.
//
// Two-point values:
//   v1 = 2 - kappa = (3 - sqrt(5)) / 2  (low value)
//   v2 = 1 + kappa = (3 + sqrt(5)) / 2  (high value)
//   P(V* = v1) = kappa / sqrt(5)
//   P(V* = v2) = 1 - kappa / sqrt(5)
//
// where kappa = (sqrt(5) + 1) / 2 (golden ratio).
//
// Parameters (modified in place):
//   kappa  - golden ratio constant
//   prob_lo - probability of the low value v1
//   v_lo   - low value (2 - kappa)
//   v_hi   - high value (1 + kappa)
//
// Returns: void
// -----------------------------------------------------------------------------
void dh_boot_mammen_params(real scalar kappa,
                           real scalar prob_lo,
                           real scalar v_lo,
                           real scalar v_hi)
{
    kappa   = (sqrt(5) + 1) / 2
    prob_lo = kappa / sqrt(5)
    v_lo    = 2 - kappa
    v_hi    = 1 + kappa
}


// =============================================================================
// Section 2: Weight Generation
// =============================================================================

// -----------------------------------------------------------------------------
// dh_boot_generate_weights()
//
// Generate bootstrap weight matrix for multiplier bootstrap.
// Returns an n x B matrix where each column is an independent draw of n weights.
//
// Parameters:
//   n      - sample size (positive integer)
//   B      - number of bootstrap iterations (positive integer)
//   method - weight distribution: "mammen" (default) | "rademacher"
//
// Returns:
//   n x B matrix of bootstrap weights
//
// Notes:
//   - "mammen": Mammen (1993) two-point distribution, E[V*]=1, Var[V*]=1
//   - "rademacher": V* in {0, 2} with equal probability (E=1, Var=1)
//     (Reserved for future extension)
// -----------------------------------------------------------------------------
real matrix dh_boot_generate_weights(real scalar n, real scalar B,
                                     string scalar method)
{
    real scalar kappa, prob_lo, v_lo, v_hi
    real matrix u, weights

    // Input validation
    if (missing(n) | missing(B)) {
        _error("dh_boot_generate_weights: n and B must be non-missing")
    }
    if (n <= 0 | B <= 0) {
        _error("dh_boot_generate_weights: n and B must be positive")
    }

    // Truncate to integers
    n = floor(n)
    B = floor(B)

    if (method == "mammen" | method == "") {
        // Mammen two-point distribution
        dh_boot_mammen_params(kappa, prob_lo, v_lo, v_hi)

        // Generate n x B uniform random matrix
        u = runiform(n, B)

        // Vectorized two-point distribution
        weights = (u :< prob_lo) :* v_lo + (u :>= prob_lo) :* v_hi
    }
    else if (method == "rademacher") {
        // Rademacher: V* in {0, 2} with P=0.5 each => E=1, Var=1
        u = runiform(n, B)
        weights = (u :< 0.5) :* 0 + (u :>= 0.5) :* 2
    }
    else {
        _error("dh_boot_generate_weights: unknown method '" + method + "'")
    }

    return(weights)
}


// =============================================================================
// Section 3: Quantile Computation
// =============================================================================

// -----------------------------------------------------------------------------
// dh_boot_quantile()
//
// Compute quantile from bootstrap statistic distribution with missing handling.
// Uses linear interpolation between order statistics (Type 7 quantile).
//
// Parameters:
//   stats - B x 1 column vector of bootstrap statistics (may contain missing)
//   alpha - quantile level in (0, 1); e.g., 0.95 for 95th percentile
//
// Returns:
//   scalar quantile value; returns missing (.) if no valid observations
// -----------------------------------------------------------------------------
real scalar dh_boot_quantile(real colvector stats, real scalar alpha)
{
    real colvector x_valid
    real scalar B, h, lo, hi

    // Remove missing values
    x_valid = select(stats, stats :< .)
    if (rows(x_valid) == 0) return(.)

    _sort(x_valid, 1)
    B = rows(x_valid)

    if (B == 1) return(x_valid[1])

    // Linear interpolation (Type 7)
    h  = (B - 1) * alpha + 1
    lo = floor(h)
    hi = ceil(h)

    if (lo < 1) lo = 1
    if (hi > B) hi = B

    if (lo == hi) return(x_valid[lo])
    return(x_valid[lo] + (h - lo) * (x_valid[hi] - x_valid[lo]))
}


// =============================================================================
// Section 4: Supremum Operations
// =============================================================================

// -----------------------------------------------------------------------------
// dh_boot_max_nonmissing()
//
// Compute maximum of non-missing values in a vector.
// Returns missing (.) if all values are missing.
//
// Parameters:
//   x - real column vector (may contain missing values)
//
// Returns:
//   scalar maximum of non-missing elements
// -----------------------------------------------------------------------------
real scalar dh_boot_max_nonmissing(real colvector x)
{
    real colvector x_valid

    x_valid = select(x, x :< .)
    if (rows(x_valid) == 0) return(.)
    return(max(x_valid))
}


// -----------------------------------------------------------------------------
// dh_boot_compute_sup()
//
// Extract critical values from bootstrap supremum t-statistic distribution.
// Implements two strategies:
//   uniformall=1: Global sup across all columns, then quantile
//   uniformall=0: Per-column quantile
//
// Parameters:
//   T_star     - B x K matrix of bootstrap sup-t statistics
//                (B iterations x K groups/pairs)
//   uniformall - flag: 1 = uniform across all K, 0 = per-column
//   alpha      - significance level (e.g., 0.05)
//
// Returns:
//   K x 1 column vector of critical values
// -----------------------------------------------------------------------------
real colvector dh_boot_compute_sup(real matrix T_star,
                                   real scalar uniformall,
                                   real scalar alpha)
{
    real scalar B, K, b, k, c_val
    real colvector c_check, mb_sup_all

    B = rows(T_star)
    K = cols(T_star)
    c_check = J(K, 1, .)

    if (uniformall) {
        // Global: max across all K columns per iteration, then quantile
        mb_sup_all = J(B, 1, .)
        for (b = 1; b <= B; b++) {
            mb_sup_all[b] = dh_boot_max_nonmissing(T_star[b, .]')
        }
        c_val = dh_boot_quantile(mb_sup_all, 1 - alpha)
        c_check = J(K, 1, c_val)
    }
    else {
        // Per-column: quantile of each column independently
        for (k = 1; k <= K; k++) {
            c_check[k] = dh_boot_quantile(T_star[., k], 1 - alpha)
        }
    }

    return(c_check)
}


end
