function metrics = extract_spectrum_metrics(spectrum, cfg)
%EXTRACT_SPECTRUM_METRICS Peak position, FWHM and HWHM for each component.

    E = spectrum.energy_meV;
    metrics = struct();
    metrics.total = profile_fwhm(E, spectrum.total_raw);

    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        metrics.(mk).total = profile_fwhm(E, ...
            spectrum.mechanism.(mk).total_raw);
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d', order);
            metrics.(mk).(ok).emission = profile_fwhm(E, ...
                spectrum.mechanism.(mk).(ok).emission_raw);
            metrics.(mk).(ok).absorption = profile_fwhm(E, ...
                spectrum.mechanism.(mk).(ok).absorption_raw);
            metrics.(mk).(ok).total = profile_fwhm(E, ...
                spectrum.mechanism.(mk).(ok).total_raw);
        end
    end
end
