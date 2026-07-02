clear all
set more off
mata: mata clear
do "mata/didhetero_types.mata"
do "mata/didhetero_errors.mata"
do "mata/didhetero_kernel.mata"
do "mata/didhetero_utils_formula.mata"
do "mata/didhetero_utils_numerical.mata"
do "mata/didhetero_utils_domain.mata"
do "mata/didhetero_utils_init.mata"
do "mata/didhetero_lpr.mata"
do "mata/didhetero_bwselect_lp.mata"

// Test: explicit brute force vs fast for n=500 and n=1000
mata:
// Inline brute force (no dispatch)
real colvector _brute_ref(real colvector Y, real colvector X, real scalar J_nn)
{
    real scalar n, i, j, count, ss, Ybar_nn
    real colvector sigma2, dist_i, Y_nn
    real matrix sort_order
    n = rows(X)
    sigma2 = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        dist_i = abs(X :- X[i])
        dist_i[i] = .
        sort_order = order(dist_i, 1)
        Y_nn = J(J_nn, 1, .)
        count = 0
        for (j = 1; j <= n; j++) {
            if (dist_i[sort_order[j]] < .) {
                count++
                Y_nn[count] = Y[sort_order[j]]
                if (count >= J_nn) break
            }
        }
        Ybar_nn = mean(Y_nn[1..count])
        ss = 0
        for (j = 1; j <= count; j++) ss = ss + (Y_nn[j] - Ybar_nn)^2
        sigma2[i] = ss / count
    }
    return(sigma2)
}
end

// n=500, J=3
mata:
X = rnormal(500, 1, 0, 1)
Y = rnormal(500, 1, 0, 1)
rb = _brute_ref(Y, X, 3)
rf = _didhetero_nn_variance_fast(Y, X, 3)
st_numscalar("d500", max(abs(rb - rf)))
end
scalar list d500
assert d500 <= 1e-12

// n=1000, J=5
mata:
X = rnormal(1000, 1, 0, 1)
Y = rnormal(1000, 1, 0, 1)
rb = _brute_ref(Y, X, 5)
rf = _didhetero_nn_variance_fast(Y, X, 5)
st_numscalar("d1000", max(abs(rb - rf)))
end
scalar list d1000
assert d1000 <= 1e-12

// edge: n=5, J=3
mata:
X = (1\3\5\7\9)
Y = (10\20\30\40\50)
rb = _brute_ref(Y, X, 3)
rf = _didhetero_nn_variance_fast(Y, X, 3)
st_numscalar("d5", max(abs(rb - rf)))
end
scalar list d5
assert d5 <= 1e-12

// edge: n=50, J=49
mata:
X = rnormal(50, 1, 0, 1)
Y = rnormal(50, 1, 0, 1)
rb = _brute_ref(Y, X, 49)
rf = _didhetero_nn_variance_fast(Y, X, 49)
st_numscalar("d49", max(abs(rb - rf)))
end
scalar list d49
assert d49 <= 1e-12

display as result "ALL TESTS PASSED"
