function spectrum = normalize_spectrum_for_plot(spectrum, cfg)
%NORMALIZE_SPECTRUM_FOR_PLOT Keep raw results and add plotting arrays.

    scale = 1.0;
    if strcmpi(cfg.oap.plot_normalization, 'global_max')
        scale = max(spectrum.total_raw);
        if ~isfinite(scale) || scale <= 0
            scale = 1.0;
        end
    elseif ~strcmpi(cfg.oap.plot_normalization, 'none')
        error('Unknown plot_normalization option.');
    end

    spectrum.plot_scale = scale;
    spectrum.total_plot = spectrum.total_raw / scale;

    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        spectrum.mechanism.(mk).total_plot = ...
            spectrum.mechanism.(mk).total_raw / scale;
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            spectrum.mechanism.(mk).(ok).emission_plot = ...
                spectrum.mechanism.(mk).(ok).emission_raw/scale;
            spectrum.mechanism.(mk).(ok).absorption_plot = ...
                spectrum.mechanism.(mk).(ok).absorption_raw/scale;
            spectrum.mechanism.(mk).(ok).total_plot = ...
                spectrum.mechanism.(mk).(ok).total_raw/scale;
        end
    end
end
