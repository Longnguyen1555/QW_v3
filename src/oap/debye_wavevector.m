function [qd, ne3d] = debye_wavevector(sp, cfg)
%DEBYE_WAVEVECTOR qd^2 = ne e^2/(eps_s eps0 kBT).

    c = cfg.constants;

    switch lower(cfg.oap.screening_density_mode)
        case 'from_sheet'
            thickness = cfg.oap.active_thickness_nm * c.nm;
            if thickness <= 0
                thickness = max(cfg.structure.Lz_nm*c.nm, sp.dz_m);
            end
            ne3d = sp.Nd_sheet_m2 / thickness;
        case 'fixed'
            ne3d = cfg.oap.fixed_ne_cm3 * 1.0e6;
        otherwise
            error('Unknown screening_density_mode.');
    end

    qd = sqrt(max(ne3d,0) * c.e^2 / ...
        (cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
end
