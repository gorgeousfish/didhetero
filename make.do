// =============================================================================
// make.do — Build script for didhetero-stata
//
// Usage:
//   do make.do          Compile Mata code and create ldidhetero.mlib
//
// Prerequisites:
//   Working directory must be set to didhetero-stata/
// =============================================================================

// Clear memory and set working directory
clear all
set more off

// Determine package root from do-file location
local PKGROOT ""
local DOFILE `"`c(filename)'"'
if strpos(`"`DOFILE'"', "/make.do") {
    local PKGROOT : subinstr local DOFILE "/make.do" "", all
}

if `"`PKGROOT'"' == "" {
    local _dh_env_root : env DIDHETERO_STATA_ROOT
    if `"`_dh_env_root'"' != "" {
        capture confirm file "`_dh_env_root'/make.do"
        if !_rc {
            local PKGROOT `"`_dh_env_root'"'
        }
    }
}

if `"`PKGROOT'"' == "" {
    capture confirm file "make.do"
    if !_rc {
        local PKGROOT `"`c(pwd)'"'
    }
    else {
        capture confirm file "didhetero-stata/make.do"
        if !_rc {
            local PKGROOT `"`c(pwd)'/didhetero-stata"'
        }
    }
}

if `"`PKGROOT'"' == "" {
    di as error "Cannot locate didhetero-stata root. Set DIDHETERO_STATA_ROOT or invoke make.do from the package/repo root."
    exit 601
}

quietly cd "`PKGROOT'"

capture confirm file "make.do"
if _rc {
    di as error "Resolved package root does not contain make.do: `PKGROOT'"
    exit 601
}

// Load test environment if available
capture noisily do "tests/_dh_test_env.do"
if _rc {
    capture noisily do "didhetero-stata/tests/_dh_test_env.do"
}
if _rc {
    display as text "(tests/_dh_test_env.do not found — skipping QA env setup)"
    local _rc 0
}

display _n "{hline 60}"
display "  didhetero-stata build system"
display "{hline 60}"

// =========================================================================
// Stage 1: Clear Mata memory
// =========================================================================
display _n "Stage 1: Clearing Mata memory..."
mata: mata clear

// =========================================================================
// Stage 2: Clean old .mlib files
// =========================================================================
display "Stage 2: Cleaning old .mlib files..."
capture erase "mata/ldidhetero.mlib"
capture erase "ado/ldidhetero.mlib"

// =========================================================================
// Stage 3: Compile .mata files in dependency order
// =========================================================================
display _n "Stage 3: Compiling Mata source files..."

// Data structure definitions
display as text "  Compiling: mata/didhetero_types.mata"
capture noisily do "mata/didhetero_types.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_types.mata (rc = " _rc ")"
    exit _rc
}

// Kernel weight computation
display as text "  Compiling: mata/didhetero_kernel.mata"
capture noisily do "mata/didhetero_kernel.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_kernel.mata (rc = " _rc ")"
    exit _rc
}

// Utility functions
display as text "  Compiling: mata/didhetero_utils.mata"
capture noisily do "mata/didhetero_utils.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_utils.mata (rc = " _rc ")"
    exit _rc
}

// Local polynomial regression
display as text "  Compiling: mata/didhetero_lpr.mata"
capture noisily do "mata/didhetero_lpr.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_lpr.mata (rc = " _rc ")"
    exit _rc
}

// Bandwidth selection
display as text "  Compiling: mata/didhetero_bwselect.mata"
capture noisily do "mata/didhetero_bwselect.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_bwselect.mata (rc = " _rc ")"
    exit _rc
}

// Kernel density estimation
display as text "  Compiling: mata/didhetero_kde.mata"
capture noisily do "mata/didhetero_kde.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_kde.mata (rc = " _rc ")"
    exit _rc
}

// Estimation modules

// Generalized propensity score estimation
display as text "  Compiling: mata/didhetero_gps.mata"
capture noisily do "mata/didhetero_gps.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_gps.mata (rc = " _rc ")"
    exit _rc
}

// Outcome regression estimation
display as text "  Compiling: mata/didhetero_or.mata"
capture noisily do "mata/didhetero_or.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_or.mata (rc = " _rc ")"
    exit _rc
}

// Stage 1 estimation dispatch
display as text "  Compiling: mata/didhetero_stage1.mata"
capture noisily do "mata/didhetero_stage1.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_stage1.mata (rc = " _rc ")"
    exit _rc
}

// Intermediate variable construction
display as text "  Compiling: mata/didhetero_intermediate.mata"
capture noisily do "mata/didhetero_intermediate.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_intermediate.mata (rc = " _rc ")"
    exit _rc
}

// CATT estimation core
display as text "  Compiling: mata/didhetero_catt_core.mata"
capture noisily do "mata/didhetero_catt_core.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_catt_core.mata (rc = " _rc ")"
    exit _rc
}
// Stage 2 and 3 estimation
display as text "  Compiling: mata/didhetero_stage23.mata"
capture noisily do "mata/didhetero_stage23.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_stage23.mata (rc = " _rc ")"
    exit _rc
}
// Standard error and analytical uniform confidence band
display as text "  Compiling: mata/didhetero_se.mata"
capture noisily do "mata/didhetero_se.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_se.mata (rc = " _rc ")"
    exit _rc
}

// Bootstrap uniform confidence band
display as text "  Compiling: mata/didhetero_bootstrap.mata"
capture noisily do "mata/didhetero_bootstrap.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_bootstrap.mata (rc = " _rc ")"
    exit _rc
}

// Bootstrap weight generation
display as text "  Compiling: mata/didhetero_boot.mata"
capture noisily do "mata/didhetero_boot.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_boot.mata (rc = " _rc ")"
    exit _rc
}

// Bootstrap optimization
display as text "  Compiling: mata/didhetero_bootstrap_opt.mata"
capture noisily do "mata/didhetero_bootstrap_opt.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_bootstrap_opt.mata (rc = " _rc ")"
    exit _rc
}


// Estimation orchestration
display as text "  Compiling: mata/didhetero_run.mata"
capture noisily do "mata/didhetero_run.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_run.mata (rc = " _rc ")"
    exit _rc
}

// Aggregation of treatment effects
display as text "  Compiling: mata/didhetero_aggte.mata"
capture noisily do "mata/didhetero_aggte.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_aggte.mata (rc = " _rc ")"
    exit _rc
}

// Simulated data generation
display as text "  Compiling: mata/didhetero_simdata.mata"
capture noisily do "mata/didhetero_simdata.mata"
if _rc {
    display as error "ERROR: Failed to compile didhetero_simdata.mata (rc = " _rc ")"
    exit _rc
}

// =========================================================================
// Stage 4: Create .mlib library
// =========================================================================
display _n "Stage 4: Creating ldidhetero.mlib..."
capture noisily mata: mata mlib create ldidhetero, dir("mata") replace
if _rc {
    display as error "ERROR: Failed to create ldidhetero.mlib (rc = " _rc ")"
    exit _rc
}

capture noisily mata: mata mlib add ldidhetero *(), dir("mata") complete
if _rc {
    display as error "ERROR: Failed to add functions to ldidhetero.mlib (rc = " _rc ")"
    exit _rc
}

// =========================================================================
// Stage 5: Copy .mlib to ado/ directory
// =========================================================================
display "Stage 5: Copying .mlib to ado/ directory..."
capture noisily copy "mata/ldidhetero.mlib" "ado/ldidhetero.mlib", replace
if _rc {
    display as error "ERROR: Failed to copy ldidhetero.mlib to ado/ (rc = " _rc ")"
    exit _rc
}

// =========================================================================
// Stage 6: Rebuild Mata library index
// =========================================================================
display "Stage 6: Rebuilding Mata library index..."
mata: mata mlib index

// =========================================================================
// Build complete
// =========================================================================
display _n as text "{hline 50}"
display as result "Build successful!"
display as text "  Library: mata/ldidhetero.mlib"
display as text "  Copy:    ado/ldidhetero.mlib"
display as text "{hline 50}"

// =========================================================================
// Stage 7: Test entry point removed
// =========================================================================
if "`0'" == "test" {
    display as error _n "ERROR: make.do no longer accepts the test subcommand."
    display as text "Build and test entry points are now decoupled."
    display as text "Use an explicit runner instead, for example:"
    display as text "  do tests/_run_tests.do"
    display as text "  do tests/run_test_gps.do"
    exit 198
}
