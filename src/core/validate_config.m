function validate_config(cfg)
%VALIDATE_CONFIG Fail early on nonphysical or unsupported inputs.

    assert(cfg.structure.Nz >= 50, 'Nz must be at least 50.');
    assert(cfg.structure.n_states >= 2, 'At least two states are required.');
    assert(cfg.structure.Lz_nm > 0, 'Lz must be positive.');
    assert(cfg.structure.U0_meV > 0, 'U0 must be positive.');
    assert(cfg.temperature_K > 0, 'Temperature must be positive.');
    assert(cfg.doping.Nd_sheet_cm2 >= 0, 'Sheet density cannot be negative.');
    assert(cfg.oap.Nqz >= 8 && cfg.oap.Nqperp >= 8, ...
        'q grids are too small.');
    assert(all(diff(cfg.oap.photon_energy_meV) > 0), ...
        'Photon-energy grid must be strictly increasing.');
    assert(size(cfg.oap.transitions, 2) == 2, ...
        'cfg.oap.transitions must have two columns [initial final].');
    assert(all(cfg.oap.transitions(:) >= 1), ...
        'Transition indices are MATLAB 1-based and must be >= 1.');
    assert(max(cfg.oap.transitions(:)) <= cfg.structure.n_states, ...
        'A transition index exceeds structure.n_states.');
    assert(all(cfg.oap.photon_orders >= 1), ...
        'Photon orders must be positive integers.');
end
