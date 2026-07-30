function save_case_outputs(case_result, cfg, tag, figures)
%SAVE_CASE_OUTPUTS Save MAT, CSV and figures for one case.

    ensure_directory(cfg.output.directory);
    safe_tag = regexprep(tag, '[^A-Za-z0-9_-]', '_');

    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, [safe_tag '.mat']), ...
             'case_result', '-v7.3');
    end

    if cfg.output.save_csv
        E = case_result.spectrum.energy_meV(:);
        total = case_result.spectrum.total_raw(:);
        optical = case_result.spectrum.mechanism.optical.total_raw(:);
        piezo = case_result.spectrum.mechanism.piezoelectric.total_raw(:);
        T = table(E,total,optical,piezo, ...
            'VariableNames',{'PhotonEnergy_meV','TotalRaw', ...
                             'OpticalRaw','PiezoelectricRaw'});
        writetable(T, fullfile(cfg.output.directory, ...
                   [safe_tag '_spectrum.csv']));

        write_metrics_csv(case_result.metrics, cfg, ...
            fullfile(cfg.output.directory,[safe_tag '_metrics.csv']));
    end

    if cfg.output.save_figures
        names = {'structure','spectrum','contributions'};
        for i = 1:min(numel(figures),numel(names))
            if ishghandle(figures{i})
                export_figure_compat(figures{i}, ...
                    fullfile(cfg.output.directory, ...
                    [safe_tag '_' names{i} '.' cfg.output.figure_format]), ...
                    cfg.output.figure_dpi);
            end
        end
    end
end
