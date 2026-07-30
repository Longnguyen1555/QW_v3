function linewidth = collect_linewidth_metrics(cases, values, cfg)
%COLLECT_LINEWIDTH_METRICS Build arrays for FWHM/HWHM sweep plots.

    linewidth.values = values(:).';
    mech_names = {'optical','piezoelectric'};

    for im = 1:numel(mech_names)
        mk = mech_names{im};
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            linewidth.(mk).(ok).emission_fwhm = nan(size(values));
            linewidth.(mk).(ok).absorption_fwhm = nan(size(values));
            linewidth.(mk).(ok).emission_hwhm = nan(size(values));
            linewidth.(mk).(ok).absorption_hwhm = nan(size(values));

            for i = 1:numel(cases)
                a = cases{i}.metrics.(mk).(ok).emission;
                b = cases{i}.metrics.(mk).(ok).absorption;
                linewidth.(mk).(ok).emission_fwhm(i) = a.fwhm;
                linewidth.(mk).(ok).absorption_fwhm(i) = b.fwhm;
                linewidth.(mk).(ok).emission_hwhm(i) = a.hwhm;
                linewidth.(mk).(ok).absorption_hwhm(i) = b.hwhm;
            end
        end
    end
end
