function results = run_qw_moap(cfg)
%RUN_QW_MOAP Top-level dispatcher.

    validate_config(cfg);

    ensure_directory(cfg.output.directory);


    task = lower(cfg.run.task);


    results = struct();

    results.config = cfg;


    switch task

        % ============================================================
        % SINGLE CASE
        % ============================================================
        case 'single'

            results.single = run_single_case( ...
                cfg, ...
                'single');


        % ============================================================
        % SPECTRUM SWEEP
        % ============================================================
        case 'sweep_spectra'

            results.sweep = run_parameter_sweep( ...
                cfg, ...
                false);


        % ============================================================
        % LINEWIDTH SWEEP
        % ============================================================
        case 'linewidth_sweep'

            results.linewidth = run_parameter_sweep( ...
                cfg, ...
                true);


        % ============================================================
        % ALL DEMO
        % ============================================================
        case 'all_demo'

            % --------------------------------------------------------
            % Single case
            % --------------------------------------------------------
            results.single = run_single_case( ...
                cfg, ...
                'single');


            % --------------------------------------------------------
            % Expensive SP/OAP sweep:
            % calculate only once
            % --------------------------------------------------------
            results.sweep = run_parameter_sweep( ...
                cfg, ...
                false);


            % --------------------------------------------------------
            % Reuse the same sweep cases for linewidth
            % --------------------------------------------------------
            lw = collect_linewidth_metrics( ...
                results.sweep.cases, ...
                results.sweep.values, ...
                cfg);


            results.linewidth = struct();

            results.linewidth.parameter = ...
                cfg.sweep.parameter;

            results.linewidth.values = ...
                cfg.sweep.values;

            results.linewidth.linewidth = lw;


            % --------------------------------------------------------
            % Separate linewidth figures
            % --------------------------------------------------------
            figures = {};
            figure_names = {};


            mechanisms = cfg.oap.mechanisms;

            if ischar(mechanisms)
                mechanisms = {mechanisms};
            end


            for im = 1:numel(mechanisms)

                mk = lower(mechanisms{im});


                if ~ismember( ...
                        mk, ...
                        {'optical', 'piezoelectric'})

                    continue;

                end


                fig = plot_linewidth_sweep( ...
                    lw, ...
                    cfg, ...
                    mk);


                figures{end+1} = fig; %#ok<AGROW>

                figure_names{end+1} = ...
                    ['linewidth_sweep_' mk]; %#ok<AGROW>

            end


            save_sweep_outputs( ...
                results.linewidth, ...
                cfg, ...
                'linewidth_sweep', ...
                figures, ...
                figure_names);


        otherwise

            error( ...
                'Unsupported cfg.run.task: %s', ...
                cfg.run.task);

    end


    % ================================================================
    % Save global project result
    % ================================================================
    if cfg.output.save_mat

        save( ...
            fullfile( ...
                cfg.output.directory, ...
                'project_results.mat'), ...
            'results', ...
            'cfg', ...
            '-v7.3');

    end

end