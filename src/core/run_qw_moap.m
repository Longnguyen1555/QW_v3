function results = run_qw_moap(cfg)
%RUN_QW_MOAP Top-level dispatcher.

    validate_config(cfg);
    ensure_directory(cfg.output.directory);

    task = lower(cfg.run.task);
    results = struct();
    results.config = cfg;

    switch task
        case 'single'
            results.single = run_single_case(cfg, 'single');

        case 'sweep_spectra'
            results.sweep = run_parameter_sweep(cfg, false);

        case 'linewidth_sweep'
            results.linewidth = run_parameter_sweep(cfg, true);

        case 'all_demo'
            results.single = run_single_case(cfg, 'single');

            % Run the expensive SP/OAP sweep only once, then reuse all cases
            % for the linewidth curves.
            results.sweep = run_parameter_sweep(cfg, false);
            lw = collect_linewidth_metrics(results.sweep.cases, ...
                                           results.sweep.values, cfg);
            results.linewidth = struct();
            results.linewidth.parameter = cfg.sweep.parameter;
            results.linewidth.values = cfg.sweep.values;
            results.linewidth.linewidth = lw;
            fig = plot_linewidth_sweep(lw, cfg);
            save_sweep_outputs(results.linewidth, cfg, ...
                               'linewidth_sweep', {fig});

        otherwise
            error('Unsupported cfg.run.task: %s', cfg.run.task);
    end

    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'project_results.mat'), ...
             'results', 'cfg', '-v7.3');
    end
end
