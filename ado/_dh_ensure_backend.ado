*! Backend loader - ensures Mata functions are available

program define _dh_ensure_backend
    version 16.0

    // Locate package directory to determine dev vs installed mode
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
        local mata_dir = subinstr(`"`mata_dir'"', "/ado/didhetero.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'"', "/ado/catt_gt.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'"', "/ado/aggte_gt.ado", "/mata", 1)
        local mata_dir = subinstr(`"`mata_dir'"', "/ado/didhetero_simdata.ado", "/mata", 1)
        capture confirm file "`mata_dir'/didhetero_types.mata"
        if _rc {
            local mata_dir ""
        }
    }

    // If mata source directory exists (dev workspace): always compile from source
    // to guarantee code matches the latest edits, bypassing any stale .mlib
    if `"`mata_dir'"' != "" {
        // Clear stale functions (by name pattern)
        capture mata: mata drop didhetero_*()
        capture mata: mata drop _didhetero_*()
        capture mata: mata drop _gteeval_*()
        capture mata: mata drop _aggte_*()
        capture mata: mata drop _dh_*()
        capture mata: mata drop dh_boot_*()
        capture mata: mata drop DH_ERR_*()
        capture mata: mata drop dh_throw_error()
        // Clear stale struct definitions (auto-generated constructors)
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

        foreach f of local mata_files {
            local mata_file "`mata_dir'/`f'"
            capture noisily do "`mata_file'"
            if _rc {
                di as error "didhetero Mata backend failed to load `f'"
                exit 601
            }
        }
        exit
    }

    // Installed mode: rely on pre-compiled .mlib
    capture quietly mata: mata which didhetero_init_from_ado()
    local has_init = (_rc == 0)
    capture quietly mata: mata which didhetero_parse_xformula_locals()
    local has_xformula = (_rc == 0)
    capture quietly mata: mata which didhetero_polynomial()
    local has_lpr = (_rc == 0)
    capture quietly mata: mata which didhetero_run_from_ado()
    local has_run = (_rc == 0)

    if `has_init' & `has_xformula' & `has_lpr' & `has_run' {
        exit
    }

    // Try rebuilding mlib index
    capture quietly mata: mata mlib index

    capture quietly mata: mata which didhetero_init_from_ado()
    local has_init = (_rc == 0)
    capture quietly mata: mata which didhetero_polynomial()
    local has_lpr = (_rc == 0)
    capture quietly mata: mata which didhetero_run_from_ado()
    local has_run = (_rc == 0)

    if !(`has_init' & `has_lpr' & `has_run') {
        di as error "didhetero Mata backend failed to initialize"
        di as error "Please reinstall the package or run make.do to rebuild"
        exit 3499
    }

end
