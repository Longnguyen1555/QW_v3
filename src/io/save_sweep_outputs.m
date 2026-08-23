function save_sweep_outputs( ...
    sweep_result, cfg, tag, figures, figure_names)
%SAVE_SWEEP_OUTPUTS
% Save sweep MAT, CSV and named figures.

    ensure_directory(cfg.output.directory);


    % ================================================================
    % Optional figure names
    % ================================================================
    if nargin < 5 || isempty(figure_names)

        figure_names = cell(1, numel(figures));

        for i = 1:numel(figures)

            figure_names{i} = sprintf( ...
                '%s_%d', ...
                tag, ...
                i);

        end

    end


    % ================================================================
    % MAT
    % ================================================================
    if cfg.output.save_mat

        save( ...
            fullfile( ...
                cfg.output.directory, ...
                [tag '.mat']), ...
            'sweep_result', ...
            '-v7.3');

    end


    % ================================================================
    % CSV
    % ================================================================
    if cfg.output.save_csv

        if isfield(sweep_result, 'cases')

            write_spectrum_sweep_csv( ...
                sweep_result, ...
                cfg, ...
                tag);

        end


        if isfield(sweep_result, 'linewidth')

            write_linewidth_sweep_csv( ...
                sweep_result, ...
                cfg, ...
                tag);

        end

    end


    % ================================================================
    % Figures
    % ================================================================
    if cfg.output.save_figures

        nfig = min( ...
            numel(figures), ...
            numel(figure_names));


        for i = 1:nfig

            if ishghandle(figures{i})

                filename = fullfile( ...
                    cfg.output.directory, ...
                    [ ...
                        figure_names{i}, ...
                        '.', ...
                        cfg.output.figure_format ...
                    ]);


                export_figure_compat( ...
                    figures{i}, ...
                    filename, ...
                    cfg.output.figure_dpi);

            end

        end

    end

end


% ========================================================================
% Spectrum sweep CSV
% ========================================================================
function write_spectrum_sweep_csv(sweep_result, cfg, tag)

    E = sweep_result.cases{1}.spectrum.energy_meV(:);

    data = E;

    names = {'PhotonEnergy_meV'};


    % ================================================================
    % Total spectrum
    % ================================================================
    for i = 1:numel(sweep_result.cases)

        data(:,end+1) = ...
            sweep_result.cases{i}.spectrum.total_raw(:); %#ok<AGROW>


        names{end+1} = matlab.lang.makeValidName( ...
            sprintf( ...
                'Total_%s_%g', ...
                cfg.sweep.parameter, ...
                sweep_result.values(i))); %#ok<AGROW>

    end


    % ================================================================
    % Optical spectrum
    % ================================================================
    for i = 1:numel(sweep_result.cases)

        data(:,end+1) = ...
            sweep_result.cases{i}.spectrum. ...
            mechanism.optical.total_raw(:); %#ok<AGROW>


        names{end+1} = matlab.lang.makeValidName( ...
            sprintf( ...
                'Optical_%s_%g', ...
                cfg.sweep.parameter, ...
                sweep_result.values(i))); %#ok<AGROW>

    end


    % ================================================================
    % Piezoelectric spectrum
    % ================================================================
    for i = 1:numel(sweep_result.cases)

        data(:,end+1) = ...
            sweep_result.cases{i}.spectrum. ...
            mechanism.piezoelectric.total_raw(:); %#ok<AGROW>


        names{end+1} = matlab.lang.makeValidName( ...
            sprintf( ...
                'Piezoelectric_%s_%g', ...
                cfg.sweep.parameter, ...
                sweep_result.values(i))); %#ok<AGROW>

    end


    T = array2table( ...
        data, ...
        'VariableNames', names);


    writetable( ...
        T, ...
        fullfile( ...
            cfg.output.directory, ...
            [tag '.csv']));

end


% ========================================================================
% Linewidth sweep CSV
% ========================================================================
function write_linewidth_sweep_csv(sweep_result, cfg, tag)

    lw = sweep_result.linewidth;


    parameter_value = [];
    mechanism = {};
    order_col = [];
    process = {};
    fwhm = [];
    hwhm = [];


    mechs = { ...
        'optical', ...
        'piezoelectric'};


    for i = 1:numel(lw.values)

        for im = 1:numel(mechs)

            mk = mechs{im};


            for order = cfg.oap.photon_orders

                ok = sprintf( ...
                    'order_%d', ...
                    order);


                % ----------------------------------------------------
                % Emission
                % ----------------------------------------------------
                parameter_value(end+1,1) = ...
                    lw.values(i); %#ok<AGROW>

                mechanism{end+1,1} = ...
                    mk; %#ok<AGROW>

                order_col(end+1,1) = ...
                    order; %#ok<AGROW>

                process{end+1,1} = ...
                    'emission'; %#ok<AGROW>

                fwhm(end+1,1) = ...
                    lw.(mk).(ok).emission_fwhm(i); %#ok<AGROW>

                hwhm(end+1,1) = ...
                    lw.(mk).(ok).emission_hwhm(i); %#ok<AGROW>


                % ----------------------------------------------------
                % Absorption
                % ----------------------------------------------------
                parameter_value(end+1,1) = ...
                    lw.values(i); %#ok<AGROW>

                mechanism{end+1,1} = ...
                    mk; %#ok<AGROW>

                order_col(end+1,1) = ...
                    order; %#ok<AGROW>

                process{end+1,1} = ...
                    'absorption'; %#ok<AGROW>

                fwhm(end+1,1) = ...
                    lw.(mk).(ok).absorption_fwhm(i); %#ok<AGROW>

                hwhm(end+1,1) = ...
                    lw.(mk).(ok).absorption_hwhm(i); %#ok<AGROW>

            end

        end

    end


    T = table( ...
        parameter_value, ...
        mechanism, ...
        order_col, ...
        process, ...
        fwhm, ...
        hwhm, ...
        'VariableNames', { ...
            cfg.sweep.parameter, ...
            'Mechanism', ...
            'PhotonOrder', ...
            'Process', ...
            'FWHM_meV', ...
            'HWHM_meV' ...
        });


    writetable( ...
        T, ...
        fullfile( ...
            cfg.output.directory, ...
            [tag '.csv']));

end