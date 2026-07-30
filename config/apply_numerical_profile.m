function cfg = apply_numerical_profile(cfg, profile)
%APPLY_NUMERICAL_PROFILE Set grid sizes without changing physical inputs.

    switch lower(profile)
        case 'quick'
            cfg.structure.Nz = 500;
            cfg.structure.n_states = min(cfg.structure.n_states, 4);
            cfg.sp.max_iter = 160;
            cfg.oap.Nqz = 38;
            cfg.oap.Nqperp = 38;
            cfg.oap.photon_energy_meV = 2.0:0.25:160.0;
        case 'standard'
            cfg.structure.Nz = 1200;
            cfg.sp.max_iter = 300;
            cfg.oap.Nqz = 90;
            cfg.oap.Nqperp = 90;
            cfg.oap.photon_energy_meV = 2.0:0.10:160.0;
        case 'high_accuracy'
            cfg.structure.Nz = 2000;
            cfg.sp.max_iter = 500;
            cfg.oap.Nqz = 150;
            cfg.oap.Nqperp = 150;
            cfg.oap.photon_energy_meV = 2.0:0.04:180.0;
        otherwise
            error('Unknown numerical profile: %s', profile);
    end
end
