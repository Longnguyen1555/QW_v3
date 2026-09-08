function [a0, Omega] = analytical_laser_displacement(Ephot, cfg)
%ANALYTICAL_LASER_DISPLACEMENT Laser displacement used by Eq. (5)/(7).

    c = cfg.constants;
    E0 = cfg.laser.E0_kVcm * 1.0e5;
    mstar = cfg.material.mstar_rel * c.m0;
    Omega = Ephot / c.hbar;

    switch lower(cfg.laser.a0_mode)
        case 'constant'
            a0 = cfg.laser.a0_nm * c.nm .* ones(size(Ephot));
        case 'dynamic'
            if any(~isfinite(Ephot(:)) | Ephot(:) <= 0)
                error('QW:AnalyticalPhotonEnergy', ...
                    'Dynamic a0 requires finite, strictly positive photon energies.');
            end
            a0 = c.e * E0 ./ (mstar * Omega.^2);
        otherwise
            error('QW:AnalyticalA0Mode', ...
                'Unknown laser.a0_mode: %s', cfg.laser.a0_mode);
    end

    if any(~isfinite(a0(:)))
        error('QW:AnalyticalA0Nonfinite', ...
            'Analytical-series laser displacement a0 is non-finite.');
    end
end
