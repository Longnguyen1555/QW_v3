function populations = electron_populations_for_oap(sp, cfg)
%ELECTRON_POPULATIONS_FOR_OAP Initial-state sheet populations.

    c = cfg.constants;
    E = sp.E_J(:);

    switch lower(cfg.oap.population_mode)
        case 'maxwell_boltzmann'
            weights = exp(-(E-min(E))/(c.kB*cfg.temperature_K));
            populations = sp.Nd_sheet_m2 .* weights ./ sum(weights);

        case 'fermi_dirac'
            populations = sp.subband_sheet_m2(:);

        case 'normalized'
            weights = exp(-(E-min(E))/(c.kB*cfg.temperature_K));
            populations = weights ./ sum(weights);

        otherwise
            error('Unsupported population_mode: %s', ...
                cfg.oap.population_mode);
    end
end
