function fig = plot_sweep_spectra(cases, labels, cfg, mechanism)
%PLOT_SWEEP_SPECTRA
% Publication-style overlay of MOAP spectra for a parameter sweep.
%
% Visual convention:
%
%   Color + marker -> sweep-parameter value
%   Line style     -> solid for every case
%
% Only three representative sweep cases are displayed:
%
%   start -> first sweep value
%   mid   -> middle sweep value
%   end   -> last sweep value
%
% Example:
%
%   cfg.sweep.values = 2.0:1.0:10.0;
%
% Full calculated sweep:
%
%   2 3 4 5 6 7 8 9 10
%
% Displayed spectra:
%
%   2 6 10
%
% The complete sweep is still calculated outside this plotting function.
% This function only reduces the number of curves displayed in the
% MOAP sweep figure.
%
% Main plotting rules:
%
%   1) displayed sweep cases are plotted on the SAME vertical scale,
%   2) photon energy is written as hbar*Omega,
%   3) sparse numerical markers are placed on each curve,
%   4) minor ticks are enabled,
%   5) the axes box is closed,
%   6) unused photon-energy range is removed,
%   7) no title is placed inside the plotting area.


    %% ====================================================================
    %  Input validation
    % =====================================================================
    if nargin < 4

        error('plot_sweep_spectra requires a mechanism.');

    end


    mk = lower(mechanism);


    switch mk

        case 'optical'
            % Electron-LO phonon spectrum.

        case 'piezoelectric'
            % Electron-piezoelectric phonon spectrum.

        otherwise

            error( ...
                'Unsupported mechanism: %s', ...
                mechanism);

    end


    if isempty(cases)

        error('No sweep cases are available for plotting.');

    end


    n_cases_all = numel(cases);


    %% ====================================================================
    %  Select start - mid - end
    %
    %  The full sweep has already been calculated by run_parameter_sweep.
    %  Only three representative cases are selected here for visualization.
    %
    %  For an odd number of cases:
    %
    %       2:1:10
    %       -> 9 cases
    %       -> indices [1 5 9]
    %       -> values  [2 6 10]
    %
    %  For an even number of cases, the lower of the two central grid
    %  points is selected.
    %
    %  unique(...,'stable') also handles sweeps containing only one or
    %  two cases.
    % =====================================================================
    idx_start = 1;

    idx_mid = ...
        floor((n_cases_all + 1) / 2);

    idx_end = ...
        n_cases_all;


    idx_plot = ...
        unique( ...
            [idx_start, idx_mid, idx_end], ...
            'stable');


    n_cases = ...
        numel(idx_plot);


    %% ====================================================================
    %  Verify spectra and collect raw data
    %
    %  IMPORTANT:
    %
    %  For a sweep figure, relative peak heights between different cases
    %  are physically meaningful.
    %
    %  Therefore this function reconstructs the plotted curves from
    %  total_raw and, if requested, applies ONE common normalization
    %  factor to the three displayed sweep cases.
    % =====================================================================
    Ecell = cell(1, n_cases);

    Praw = cell(1, n_cases);


    common_max = 0;


    for i = 1:n_cases

        case_idx = ...
            idx_plot(i);


        spectrum = ...
            cases{case_idx}.spectrum;


        if ~isfield(spectrum.mechanism, mk)

            error( ...
                'Mechanism "%s" not found in sweep case %d.', ...
                mk, ...
                case_idx);

        end


        E_i = ...
            spectrum.energy_meV(:);

        P_i = ...
            spectrum.mechanism.(mk).total_raw(:);


        if numel(E_i) ~= numel(P_i)

            error( ...
                ['Photon-energy and MOAP array sizes differ ' ...
                 'in case %d.'], ...
                case_idx);

        end


        Ecell{i} = ...
            E_i;

        Praw{i} = ...
            P_i;


        finite_P = ...
            P_i(isfinite(P_i));


        if ~isempty(finite_P)

            Pi_max = ...
                max(abs(finite_P));


            if isfinite(Pi_max)

                common_max = ...
                    max( ...
                        common_max, ...
                        Pi_max);

            end

        end

    end


    %% ====================================================================
    %  Common normalization
    %
    %  One scale for start, mid and end.
    % =====================================================================
    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        common_scale = ...
            common_max;


        if ~isfinite(common_scale) || ...
                common_scale <= 0

            common_scale = ...
                1.0;

        end


    elseif strcmpi(cfg.oap.plot_normalization, 'none')

        common_scale = ...
            1.0;


    else

        error( ...
            'Unknown cfg.oap.plot_normalization: %s', ...
            cfg.oap.plot_normalization);

    end


    %% ====================================================================
    %  Figure
    % =====================================================================
    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 710 520]);


    ax = ...
        axes(fig);


    hold(ax, 'on');


    %% ====================================================================
    %  Sweep colors
    %
    %  For the three representative curves:
    %
    %       start -> magenta
    %       mid   -> blue
    %       end   -> dark red
    % =====================================================================
    sweep_colors = [ ...
        0.850 0.000 0.850; ...   % magenta
        0.050 0.150 0.650; ...   % blue
        0.600 0.050 0.000; ...   % dark red
        0.050 0.500 0.150; ...   % green
        0.100 0.100 0.100; ...   % black
        0.850 0.450 0.000; ...   % orange
        0.100 0.600 0.650  ...   % cyan
    ];


    sweep_markers = { ...
        'o', ...
        's', ...
        '^', ...
        'd', ...
        'v', ...
        '>', ...
        '<' ...
    };


    %% ====================================================================
    %  Storage
    % =====================================================================
    legend_handles = ...
        gobjects(0);

    legend_labels = ...
        {};

    x_active_all = ...
        [];

    y_all = ...
        [];


    %% ====================================================================
    %  Plot start, mid and end
    % =====================================================================
    for i = 1:n_cases

        case_idx = ...
            idx_plot(i);


        E = ...
            Ecell{i};

        P = ...
            Praw{i} ./ common_scale;


        %% ----------------------------------------------------------------
        %  Style
        % -----------------------------------------------------------------
        color_idx = ...
            mod(i - 1, size(sweep_colors, 1)) + 1;


        marker_idx = ...
            mod(i - 1, numel(sweep_markers)) + 1;


        curve_color = ...
            sweep_colors(color_idx, :);


        curve_marker = ...
            sweep_markers{marker_idx};


        %% ----------------------------------------------------------------
        %  Main spectrum
        %
        %  All lines remain solid, as in the original sweep plot.
        % -----------------------------------------------------------------
        h = plot( ...
            ax, ...
            E, ...
            P, ...
            '-', ...
            'Color', curve_color, ...
            'LineWidth', 1.65);


        legend_handles(end+1) = ...
            h; %#ok<AGROW>


        %% ----------------------------------------------------------------
        %  Legend label
        %
        %  IMPORTANT:
        %
        %  case_idx is the index in the FULL sweep.
        %
        %  Therefore:
        %
        %       cfg.sweep.values(case_idx)
        %
        %  must be used instead of cfg.sweep.values(i).
        % -----------------------------------------------------------------
        sweep_value = ...
            NaN;


        if isfield(cfg, 'sweep') && ...
                isfield(cfg.sweep, 'values') && ...
                numel(cfg.sweep.values) >= case_idx

            sweep_value = ...
                cfg.sweep.values(case_idx);

        end


        if isfinite(sweep_value)

            legend_labels{end+1} = ...
                sweep_legend_label( ...
                    cfg.sweep.parameter, ...
                    sweep_value); %#ok<AGROW>


        elseif nargin >= 2 && ...
                numel(labels) >= case_idx

            % Compatibility fallback.
            legend_labels{end+1} = ...
                labels{case_idx}; %#ok<AGROW>


        else

            legend_labels{end+1} = ...
                sprintf( ...
                    'Case %d', ...
                    case_idx); %#ok<AGROW>

        end


        %% ----------------------------------------------------------------
        %  Determine active spectral region for THIS case
        %
        %  Threshold is relative to the maximum of each case separately.
        %
        %  This ensures that a weak spectrum is not discarded merely
        %  because another displayed sweep case has a much larger
        %  absorption peak.
        % -----------------------------------------------------------------
        P_abs = ...
            abs(P);


        finite_mask = ...
            isfinite(E) & ...
            isfinite(P);


        Pfinite = ...
            P_abs(finite_mask);


        if ~isempty(Pfinite)

            Pmax_case = ...
                max(Pfinite);

        else

            Pmax_case = ...
                0;

        end


        if isfinite(Pmax_case) && ...
                Pmax_case > 0

            active_threshold = ...
                1.0e-4 * Pmax_case;


            idx_active = ...
                find( ...
                    finite_mask & ...
                    P_abs >= active_threshold);


        else

            idx_active = ...
                [];

        end


        %% ----------------------------------------------------------------
        %  Sparse numerical markers
        %
        %  Spectra commonly contain thousands of photon-energy points.
        %  Markers are therefore placed only at a small number of points
        %  in the physically active spectral interval.
        % -----------------------------------------------------------------
        if numel(idx_active) >= 2

            n_marker = ...
                min(24, numel(idx_active));


            marker_pos = ...
                unique( ...
                    round( ...
                        linspace( ...
                            1, ...
                            numel(idx_active), ...
                            n_marker)));


            idx_marker = ...
                idx_active(marker_pos);


            plot( ...
                ax, ...
                E(idx_marker), ...
                P(idx_marker), ...
                curve_marker, ...
                'LineStyle', 'none', ...
                'Color', curve_color, ...
                'MarkerEdgeColor', curve_color, ...
                'MarkerFaceColor', curve_color, ...
                'MarkerSize', 3.1, ...
                'HandleVisibility', 'off');


            x_active_all = [ ...
                x_active_all; ...
                E(idx_active) ...
            ]; %#ok<AGROW>


        else

            % If no active region can be identified, retain the complete
            % finite energy interval for safe axis determination.
            x_active_all = [ ...
                x_active_all; ...
                E(finite_mask) ...
            ]; %#ok<AGROW>

        end


        %% ----------------------------------------------------------------
        %  Store all y values for common y-axis scaling
        % -----------------------------------------------------------------
        y_all = [ ...
            y_all; ...
            P(isfinite(P)) ...
        ]; %#ok<AGROW>

    end


    %% ====================================================================
    %  X-axis label
    % =====================================================================
    xlabel( ...
        ax, ...
        'Photon energy $\hbar\Omega$ (meV)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


    %% ====================================================================
    %  Y-axis label
    %
    %  The reference OAP spectra use arbitrary units.
    % =====================================================================
    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        ylabel( ...
            ax, ...
            'Normalized MOAP', ...
            'Interpreter', 'latex', ...
            'FontSize', 13);


    else

        ylabel( ...
            ax, ...
            'MOAP (arb. units)', ...
            'Interpreter', 'latex', ...
            'FontSize', 13);

    end


    %% ====================================================================
    %  X-axis limits
    %
    %  Use the UNION of the physically active regions of the three
    %  displayed sweep spectra.
    %
    %  This is required because changing B, T, Lz, U0, etc. can shift
    %  resonance peaks substantially to the left or right.
    % =====================================================================
    x_active_all = ...
        x_active_all(isfinite(x_active_all));


    if numel(x_active_all) >= 2

        xmin_active = ...
            min(x_active_all);

        xmax_active = ...
            max(x_active_all);


        active_width = ...
            xmax_active - xmin_active;


        if active_width > 0

            xpad = ...
                0.035 * active_width;

        else

            xpad = ...
                0;

        end


        all_E_min = ...
            inf;

        all_E_max = ...
            -inf;


        for i = 1:n_cases

            Ei = ...
                Ecell{i};


            Ei = ...
                Ei(isfinite(Ei));


            if isempty(Ei)

                continue;

            end


            all_E_min = ...
                min( ...
                    all_E_min, ...
                    min(Ei));


            all_E_max = ...
                max( ...
                    all_E_max, ...
                    max(Ei));

        end


        xmin = ...
            max( ...
                all_E_min, ...
                xmin_active - xpad);


        xmax = ...
            min( ...
                all_E_max, ...
                xmax_active + xpad);


        if isfinite(xmin) && ...
                isfinite(xmax) && ...
                xmax > xmin

            xlim( ...
                ax, ...
                [xmin xmax]);

        end

    end


    %% ====================================================================
    %  Y-axis limits
    %
    %  Start, mid and end share exactly the same y-axis.
    %
    %  MOAP is non-negative; therefore zero is retained as the physical
    %  baseline whenever negative values are only numerical round-off.
    % =====================================================================
    y_all = ...
        y_all(isfinite(y_all));


    if ~isempty(y_all)

        ymax = ...
            max(y_all);

        ymin = ...
            min(y_all);


        if isfinite(ymax) && ...
                ymax > 0

            negative_tolerance = ...
                1.0e-10 * ymax;


            if ymin >= -negative_tolerance

                ylim( ...
                    ax, ...
                    [ ...
                        0, ...
                        1.055 * ymax ...
                    ]);


            else

                yrange = ...
                    ymax - ymin;

                ypad = ...
                    0.04 * yrange;


                ylim( ...
                    ax, ...
                    [ ...
                        ymin - ypad, ...
                        ymax + ypad ...
                    ]);

            end

        end

    end


    %% ====================================================================
    %  Legend
    %
    %  Parameter values, rather than generic "case 1", "case 2", etc.,
    %  are used in the same way as the original sweep spectra.
    % =====================================================================
    if ~isempty(legend_handles)

        lgd = legend( ...
            ax, ...
            legend_handles, ...
            legend_labels, ...
            'Location', 'best', ...
            'Interpreter', 'latex', ...
            'FontSize', 9.5, ...
            'Box', 'on');


        if numel(legend_labels) > 5

            lgd.NumColumns = ...
                2;

        end


        lgd.ItemTokenSize = ...
            [24 10];

    end


    %% ====================================================================
    %  Axes formatting
    %
    %  - closed rectangular frame
    %  - inward major/minor ticks
    %  - light major grid
    %  - no minor grid
    % =====================================================================
    set( ...
        ax, ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'LineWidth', 0.9, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'TickDir', 'in', ...
        'TickLength', [0.018 0.018], ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLabelInterpreter', 'latex');


    grid(ax, 'on');


    ax.GridLineStyle = ...
        '--';

    ax.GridAlpha = ...
        0.20;

    ax.XMinorGrid = ...
        'off';

    ax.YMinorGrid = ...
        'off';


    %% ====================================================================
    %  No title
    %
    %  The sweep parameter and fixed physical parameters should be stated
    %  in the figure caption, as in the reference papers.
    % =====================================================================


    %% ====================================================================
    %  Compact publication layout
    % =====================================================================
    ax.LooseInset = ...
        max( ...
            ax.TightInset, ...
            0.02);


end



% =========================================================================
%  Local helper: LaTeX legend label for sweep parameter
% =========================================================================
function label = sweep_legend_label(parameter, value)

    switch lower(parameter)

        case 'b_t'

            label = sprintf( ...
                '$B = %.3g\\,\\mathrm{T}$', ...
                value);


        case 't_k'

            label = sprintf( ...
                '$T = %.3g\\,\\mathrm{K}$', ...
                value);


        case 'lz_nm'

            label = sprintf( ...
                '$L_z = %.3g\\,\\mathrm{nm}$', ...
                value);


        case 'e_kvcm'

            label = sprintf( ...
                '$E = %.3g\\,\\mathrm{kV/cm}$', ...
                value);


        case 'nd_sheet_cm2'

            label = sprintf( ...
                '$N_D = %.3g\\,\\mathrm{cm}^{-2}$', ...
                value);


        case 'manning_prefactor'

            label = sprintf( ...
                '$\\nu = %.3g$', ...
                value);


        case 'u0_mev'

            label = sprintf( ...
                '$U_0 = %.3g\\,\\mathrm{meV}$', ...
                value);


        otherwise

            label = sprintf( ...
                '$%s = %.3g$', ...
                strrep(parameter, '_', '\_'), ...
                value);

    end

end