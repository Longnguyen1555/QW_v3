function sweep_result = run_parameter_sweep(cfg, linewidth_only)
%RUN_PARAMETER_SWEEP
% Recompute SP and MOAP while varying one parameter.
%
% The optical and piezoelectric mechanisms are plotted
% in completely separate figures.

    values = cfg.sweep.values(:).';

    n_cases = numel(values);

    cases = cell(1, n_cases);
    labels = cell(1, n_cases);


    % ================================================================
    % Compute all sweep cases
    % ================================================================
    for i = 1:n_cases

        cfg_i = set_sweep_parameter( ...
            cfg, ...
            cfg.sweep.parameter, ...
            values(i));


        % Do not save every intermediate sweep case separately
        cfg_i.output.save_figures = false;
        cfg_i.output.save_csv = false;
        cfg_i.output.save_mat = false;


        tag = sprintf( ...
            '%s_%g', ...
            cfg.sweep.parameter, ...
            values(i));


        cases{i} = run_single_case_no_plots(cfg_i);


        labels{i} = sweep_label( ...
            cfg.sweep.parameter, ...
            values(i));


        if cfg.run.verbose

            fprintf( ...
                '[SWEEP] %s: %d/%d complete\n', ...
                cfg.sweep.parameter, ...
                i, ...
                n_cases);

        end

    end


    % ================================================================
    % Store numerical sweep result
    % ================================================================
    sweep_result = struct();

    sweep_result.parameter = cfg.sweep.parameter;
    sweep_result.values = values;
    sweep_result.labels = labels;
    sweep_result.cases = cases;


    % ================================================================
    % Determine requested mechanisms
    % ================================================================
    mechanisms = cfg.oap.mechanisms;

    if ischar(mechanisms)
        mechanisms = {mechanisms};
    end


    % ================================================================
    % LINEWIDTH SWEEP
    % ================================================================
    if linewidth_only

        linewidth = collect_linewidth_metrics( ...
            cases, ...
            values, ...
            cfg);


        sweep_result.linewidth = linewidth;


        figures = {};
        figure_names = {};


        for im = 1:numel(mechanisms)

            mk = lower(mechanisms{im});


            if ~ismember( ...
                    mk, ...
                    {'optical', 'piezoelectric'})

                continue;

            end


            fig = plot_linewidth_sweep( ...
                linewidth, ...
                cfg, ...
                mk);


            figures{end+1} = fig; %#ok<AGROW>

            figure_names{end+1} = ...
                ['linewidth_sweep_' mk]; %#ok<AGROW>

        end


        save_sweep_outputs( ...
            sweep_result, ...
            cfg, ...
            'linewidth_sweep', ...
            figures, ...
            figure_names);


    % ================================================================
    % SPECTRUM SWEEP
    % ================================================================
    else

        figures = {};
        figure_names = {};


        for im = 1:numel(mechanisms)

            mk = lower(mechanisms{im});


            if ~ismember( ...
                    mk, ...
                    {'optical', 'piezoelectric'})

                continue;

            end


            fig = plot_sweep_spectra( ...
                cases, ...
                labels, ...
                cfg, ...
                mk);


            figures{end+1} = fig; %#ok<AGROW>

            figure_names{end+1} = ...
                ['sweep_spectra_' mk]; %#ok<AGROW>

        end


        save_sweep_outputs( ...
            sweep_result, ...
            cfg, ...
            'sweep_spectra', ...
            figures, ...
            figure_names);

    end

end


% ========================================================================
% Local helper
% ========================================================================
function out = run_single_case_no_plots(cfg)

    sp = solve_schrodinger_poisson(cfg);

    td = compute_transition_data( ...
        sp, ...
        cfg);

    spectrum = compute_moap_spectrum( ...
        sp, ...
        td, ...
        cfg);

    metrics = extract_spectrum_metrics( ...
        spectrum, ...
        cfg);


    out = struct( ...
        'cfg', cfg, ...
        'sp', sp, ...
        'transitions', td, ...
        'spectrum', spectrum, ...
        'metrics', metrics);

end