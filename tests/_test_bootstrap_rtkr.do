// Quick verification: bootstrap CI with precomputed RtKR optimization
// Confirms compilation + reasonable CI output
clear all
set more off

// Load package
quietly adopath ++ "/Users/cxy/Desktop/didhetero/didhetero-main/ado"
quietly mata: mata mlib index

// Generate simulated data (guaranteed well-formed)
display _n "Generating simulated data..."
didhetero_simdata, n(300) tau(4) seed(42) clear

// Run CATT with bootstrap — manual bandwidth
display _n "Running CATT with bootstrap (biters=50, seed=12345, manual bw)..."
catt_gt Y, group(G) time(period) id(id) z(Z) ///
    zeval(-0.5 0.0 0.5) ///
    gteval(2 2 3 3) ///
    porder(1) kernel("gau") bwselect("manual") bw(0.5) ///
    bstrap("true") biters(50) seed(12345) ///
    uniformall(false) ///
    control_group("notyettreated")

// Check results
matrix R = e(results)
display _n "=== Bootstrap CI Verification ==="
matrix list R, format(%9.5f)

// Check first row
display _n "Row 1 checks:"
display "  Point estimate: " R[1,4]
display "  SE:             " R[1,5]
display "  CI2 lower:      " R[1,8]
display "  CI2 upper:      " R[1,9]

// Verify CI2 (bootstrap) is not missing for at least one row
local any_ci2 = 0
forvalues i = 1/`=rowsof(R)' {
    if R[`i',8] != . & R[`i',9] != . {
        local any_ci2 = 1
    }
}
if `any_ci2' == 0 {
    display as error "FAIL: All Bootstrap CI values are missing!"
    exit 198
}

// Verify CI2 contains the point estimate (for non-missing rows)
forvalues i = 1/`=rowsof(R)' {
    if R[`i',8] != . & R[`i',9] != . {
        if R[`i',8] > R[`i',4] | R[`i',9] < R[`i',4] {
            display as error "FAIL: Row `i' - Bootstrap CI does not contain point estimate!"
            exit 198
        }
    }
}

display _n as result "PASS: Bootstrap precomputed RtKR optimization verified."
