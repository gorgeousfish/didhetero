*! Backend loader - ensures Mata functions are available

program define _dh_ensure_backend
    version 16.0

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

    // Rebuild mlib index without clearing Mata session
    capture quietly mata: mata mlib index

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

    // Fall back to source loading from adjacent mata/ directory
    local ado_file ""
    foreach probe in didhetero.ado catt_gt.ado aggte_gt.ado didhetero_simdata.ado _dh_ensure_backend.ado {
        capture quietly findfile `probe'
        if !_rc {
            local ado_file `"`r(fn)'"'
            continue, break
        }
    }

    if `"`ado_file'"' == "" {
        di as error "didhetero Mata backend loader could not locate package ado files"
        exit 601
    }

    local mata_dir = subinstr(`"`ado_file'"', "/ado/_dh_ensure_backend.ado", "/mata", 1)
    local mata_dir = subinstr(`"`mata_dir'"', "/ado/didhetero.ado", "/mata", 1)
    local mata_dir = subinstr(`"`mata_dir'"', "/ado/catt_gt.ado", "/mata", 1)
    local mata_dir = subinstr(`"`mata_dir'"', "/ado/aggte_gt.ado", "/mata", 1)
    local mata_dir = subinstr(`"`mata_dir'"', "/ado/didhetero_simdata.ado", "/mata", 1)

    local mata_files ///
        didhetero_types.mata ///
        didhetero_kernel.mata ///
        didhetero_utils.mata ///
        didhetero_lpr.mata ///
        didhetero_bwselect.mata ///
        didhetero_kde.mata ///
        didhetero_gps.mata ///
        didhetero_or.mata ///
        didhetero_stage1.mata ///
        didhetero_intermediate.mata ///
        didhetero_catt_core.mata ///
        didhetero_stage23.mata ///
        didhetero_se.mata ///
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
            di as error "didhetero Mata backend failed to load"
            exit 601
        }
    }

    capture quietly mata: mata which didhetero_init_from_ado()
    local has_init = (_rc == 0)
    capture quietly mata: mata which didhetero_parse_xformula_locals()
    local has_xformula = (_rc == 0)
    capture quietly mata: mata which didhetero_polynomial()
    local has_lpr = (_rc == 0)
    capture quietly mata: mata which didhetero_run_from_ado()
    local has_run = (_rc == 0)

    if !(`has_init' & `has_xformula' & `has_lpr' & `has_run') {
        di as error "didhetero Mata backend failed to initialize"
        exit 3499
    }

end
