*! Backend loader - ensures Mata functions are available
*! Guarantees silent loading for end users (no compilation output)

program define _dh_ensure_backend
    version 16.0

    // ─── Step 1: Check if key functions are already in memory ───────────
    capture mata: mata which didhetero_run_from_ado()
    if _rc == 0 {
        // Functions already available (mlib loaded or previously compiled)
        exit
    }

    // ─── Step 2: Try rebuilding mlib index ─────────────────────────────
    capture quietly mata: mata mlib index
    capture mata: mata which didhetero_run_from_ado()
    if _rc == 0 {
        exit
    }

    // ─── Step 3: Try to directly locate and load the mlib file ─────────
    capture findfile ldidhetero.mlib
    if _rc == 0 {
        // File exists but wasn't indexed; force load
        quietly mata: mata mlib add `"`r(fn)'"'
        capture mata: mata which didhetero_run_from_ado()
        if _rc == 0 {
            exit
        }
    }

    // ─── Step 4: Determine dev vs installed mode ───────────────────────
    // Locate package ado directory
    local ado_file ""
    foreach probe in _dh_ensure_backend.ado didhetero.ado catt_gt.ado aggte_gt.ado didhetero_simdata.ado {
        capture quietly findfile `probe'
        if !_rc {
            local ado_file `"`r(fn)'"'
            continue, break
        }
    }

    // Resolve mata source directory (exists only in dev workspace)
    local mata_dir ""
    if `"`ado_file'"' != "" {
        local mata_dir = subinstr(`"`ado_file'"', "/ado/_dh_ensure_backend.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'", "/ado/didhetero.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'", "/ado/catt_gt.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'", "/ado/aggte_gt.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'", "/ado/didhetero_simdata.ado", "/mata", 1)
        capture confirm file "`mata_dir'/didhetero_types.mata"
        if _rc {
            local mata_dir ""
        }
    }

    // ─── Step 5: No mlib found, no dev source → error ──────────────────
    if `"`mata_dir'"' == "" {
        di as error "didhetero: Mata library (ldidhetero.mlib) not found."
        di as error "Please reinstall the package:"
        di as error `"    net install didhetero, from("https://raw.githubusercontent.com/xxx/didhetero/main") replace"'
        exit 601
    }

    // ─── Step 6: Dev mode — silently compile from source ───────────────
    // Only developers who have the mata/ source tree will reach here.
    di as text "(didhetero: compiling Mata source files...)"

    // Clear stale functions
    capture mata: mata drop didhetero_*()
    capture mata: mata drop _didhetero_*()
    capture mata: mata drop _gteeval_*()
    capture mata: mata drop _aggte_*()
    capture mata: mata drop _dh_*()
    capture mata: mata drop dh_boot_*()
    capture mata: mata drop DH_ERR_*()
    capture mata: mata drop dh_throw_error()
    // Clear stale struct definitions
    capture mata: mata drop DidHeteroData()
    capture mata: mata drop DidHeteroParamResults()
    capture mata: mata drop DidHeteroStage1Results()
    capture mata: mata drop DidHeteroKernelConsts()
    capture mata: mata drop BootPrecomp()
    capture mata: mata drop DidHeteroEstResult()
    capture mata: mata drop DidHeteroAggteResult()
    capture mata: mata drop DidHeteroCattResult()
    capture mata: mata drop DidHeteroIntermediate()
    capture mata: mata drop AggtResult()

    local mata_files ///
        didhetero_types.mata ///
        didhetero_errors.mata ///
        didhetero_kernel.mata ///
        didhetero_utils_formula.mata ///
        didhetero_utils_numerical.mata ///
        didhetero_utils_domain.mata ///
        didhetero_utils_init.mata ///
        didhetero_lpr.mata ///
        didhetero_bwselect_lp.mata ///
        didhetero_bwselect_kde.mata ///
        didhetero_kde.mata ///
        didhetero_bwselect_lpdensity.mata ///
        didhetero_gps.mata ///
        didhetero_or.mata ///
        didhetero_stage1.mata ///
        didhetero_intermediate.mata ///
        didhetero_catt_core.mata ///
        didhetero_stage23.mata ///
        didhetero_se.mata ///
        didhetero_bootstrap_engine.mata ///
        didhetero_bootstrap_unified.mata ///
        didhetero_bootstrap.mata ///
        didhetero_boot.mata ///
        didhetero_bootstrap_opt.mata ///
        didhetero_run.mata ///
        didhetero_aggte.mata ///
        didhetero_simdata.mata

    // All compilation wrapped in quietly — no output to user
    quietly {
        foreach f of local mata_files {
            local mata_file "`mata_dir'/`f'"
            capture noisily do "`mata_file'"
            if _rc {
                noisily di as error "didhetero: Mata backend failed to compile `f'"
                exit 601
            }
        }
    }

    di as text "(didhetero: compilation complete)"
end
