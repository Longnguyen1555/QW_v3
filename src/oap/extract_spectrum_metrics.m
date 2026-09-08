function metrics = extract_spectrum_metrics(spectrum, cfg)
%EXTRACT_SPECTRUM_METRICS Peak position, FWHM and HWHM for each component.

    E = spectrum.energy_meV;
    use_binned_raw = isfield(spectrum, 'meta') && ...
        isfield(spectrum.meta, 'linewidth_source') && ...
        strcmpi(spectrum.meta.linewidth_source, 'binned_raw');
    metrics = struct();
    if use_binned_raw
        metrics.total = profile_fwhm(E, spectrum.total_binned_raw);
    else
        metrics.total = profile_fwhm(E, spectrum.total_raw);
    end

    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        metrics.(mk).total = profile_fwhm(E, channel_for_linewidth( ...
            spectrum.mechanism.(mk), 'total', use_binned_raw));
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d', order);
            component = spectrum.mechanism.(mk).(ok);
            metrics.(mk).(ok).emission = profile_fwhm(E, ...
                channel_for_linewidth(component, 'emission', use_binned_raw));
            metrics.(mk).(ok).absorption = profile_fwhm(E, ...
                channel_for_linewidth(component, 'absorption', use_binned_raw));
            metrics.(mk).(ok).total = profile_fwhm(E, ...
                channel_for_linewidth(component, 'total', use_binned_raw));
        end
    end
end

function y = channel_for_linewidth(channel, prefix, use_binned_raw)
    if use_binned_raw
        y = channel.([prefix '_binned_raw']);
    else
        y = channel.([prefix '_raw']);
    end
end
