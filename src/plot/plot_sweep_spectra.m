function fig = plot_sweep_spectra(cases, labels, cfg, mechanism)
%PLOT_SWEEP_SPECTRA
% Publication-style overlay of MOAP spectra for a parameter sweep.
%
% Visual convention:
%
%   Color + marker -> sweep-parameter value
%   Line style     -> solid for every case
%
% This follows the spectrum-sweep figures in the reference papers,
% especially Fig. 3, Fig. 5 and Fig. 7 of:
%
%   Micro and Nanostructures 198 (2025) 208062
%
% Main plotting rules:
%
%   1) every sweep case is plotted on the SAME vertical scale,
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


    n_cases = numel(cases);


    %% ====================================================================
    %  Verify spectra and collect raw data
    %
    %  IMPORTANT:
    %
    %  For a sweep figure, relative peak heights between different cases
    %  are physically meaningful.
    %
    %  normalize_spectrum_for_plot() currently normalizes EACH case
    %  independently when plot_normalization = 'global_max'.
    %
    %  Therefore this function intentionally reconstructs the plotted
    %  curves from total_raw and, if requested, applies ONE common
    %  normalization factor to all sweep cases.
    % =====================================================================
    Ecell = cell(1, n_cases);

    Praw = cell(1, n_cases);


    common_max = 0;


    for i = 1:n_cases

        spectrum = cases{i}.spectrum;


        if ~isfield(spectrum.mechanism, mk)

            error( ...
                'Mechanism "%s" not found in sweep case %d.', ...
                mk, ...
                i);

        end


        E_i = spectrum.energy_meV(:);

        P_i = spectrum.mechanism.(mk).total_raw(:);


        if numel(E_i) ~= numel(P_i)

            error( ...
                'Photon-energy and MOAP array sizes differ in case %d.', ...
                i);

        end


        Ecell{i} = E_i;

        Praw{i} = P_i;


        finite_P = P_i(isfinite(P_i));


        if ~isempty(finite_P)

            Pi_max = max(abs(finite_P));

            if isfinite(Pi_max)

                common_max = max( ...
                    common_max, ...
                    Pi_max);

            end

        end

    end


    %% ====================================================================
    %  Common normalization
    %
    %  One scale for ALL sweep cases.
    % =====================================================================
    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        common_scale = common_max;


        if ~isfinite(common_scale) || common_scale <= 0

            common_scale = 1.0;

        end


    elseif strcmpi(cfg.oap.plot_normalization, 'none')

        common_scale = 1.0;


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


    ax = axes(fig);

    hold(ax, 'on');


    %% ====================================================================
    %  Sweep colors
    %
    %  The first three colors are chosen close to the visual convention
    %  used in the reference OAP sweep figures:
    %
    %       magenta -> blue -> dark red
    %
    %  Additional colors allow sweeps containing more than three cases.
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
    legend_handles = gobjects(0);

    legend_labels = {};

    x_active_all = [];

    y_all = [];


    %% ====================================================================
    %  Plot every sweep case
    % =====================================================================
    for i = 1:n_cases

        E = Ecell{i};

        P = Praw{i} ./ common_scale;


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
        %  All lines remain SOLID, as in the reference parameter-sweep
        %  spectra. Line style is therefore not overloaded with another
        %  physical meaning.
        % -----------------------------------------------------------------
        h = plot( ...
            ax, ...
            E, ...
            P, ...
            '-', ...
            'Color', curve_color, ...
            'LineWidth', 1.65);


        legend_handles(end+1) = h; %#ok<AGROW>


        %% ----------------------------------------------------------------
        %  Legend label
        % -----------------------------------------------------------------
        sweep_value = NaN;


        if isfield(cfg, 'sweep') && ...
                isfield(cfg.sweep, 'values') && ...
                numel(cfg.sweep.values) >= i

            sweep_value = cfg.sweep.values(i);

        end


        if isfinite(sweep_value)

            legend_labels{end+1} = ...
                sweep_legend_label( ...
                    cfg.sweep.parameter, ...
                    sweep_value); %#ok<AGROW>

        elseif nargin >= 2 && ...
                numel(labels) >= i

            % Compatibility fallback.
            legend_labels{end+1} = ...
                labels{i}; %#ok<AGROW>

        else

            legend_labels{end+1} = ...
                sprintf('Case %d', i); %#ok<AGROW>

        end


        %% ----------------------------------------------------------------
        %  Determine active spectral region for THIS case
        %
        %  Threshold is relative to the maximum of each case separately.
        %
        %  This ensures that a weak spectrum is not discarded merely
        %  because another sweep case has a much larger absorption peak.
        % -----------------------------------------------------------------
        P_abs = abs(P);

        finite_mask = ...
            isfinite(E) & ...
            isfinite(P);


        Pfinite = P_abs(finite_mask);


        if ~isempty(Pfinite)

            Pmax_case = max(Pfinite);

        else

            Pmax_case = 0;

        end


        if isfinite(Pmax_case) && ...
                Pmax_case > 0

            active_threshold = ...
                1.0e-4 * Pmax_case;


            idx_active = find( ...
                finite_mask & ...
                P_abs >= active_threshold);


        else

            idx_active = [];

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


            marker_pos = unique(round(linspace( ...
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
    %  Use the UNION of all physically active sweep spectra.
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

            xpad = 0;

        end


        all_E_min = inf;

        all_E_max = -inf;


        for i = 1:n_cases

            Ei = Ecell{i};

            Ei = Ei(isfinite(Ei));


            if isempty(Ei)
                continue;
            end


            all_E_min = ...
                min(all_E_min, min(Ei));

            all_E_max = ...
                max(all_E_max, max(Ei));

        end


        xmin = max( ...
            all_E_min, ...
            xmin_active - xpad);

        xmax = min( ...
            all_E_max, ...
            xmax_active + xpad);


        if isfinite(xmin) && ...
                isfinite(xmax) && ...
                xmax > xmin

            xlim(ax, [xmin xmax]);

        end

    end


    %% ====================================================================
    %  Y-axis limits
    %
    %  All sweep cases share exactly the same y-axis.
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

                ylim(ax, [ ...
                    0, ...
                    1.055 * ymax ...
                ]);


            else

                yrange = ...
                    ymax - ymin;

                ypad = ...
                    0.04 * yrange;


                ylim(ax, [ ...
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
    %  are used in the same way as the reference sweep spectra.
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

            lgd.NumColumns = 2;

        end


        lgd.ItemTokenSize = [24 10];

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

    ax.GridLineStyle = '--';

    ax.GridAlpha = 0.20;

    ax.XMinorGrid = 'off';

    ax.YMinorGrid = 'off';


    %% ====================================================================
    %  No title
    %
    %  The sweep parameter and the fixed physical parameters should be
    %  stated in the figure caption, as in the reference papers.
    % =====================================================================


    %% ====================================================================
    %  Compact publication layout
    % =====================================================================
    ax.LooseInset = ...
        max(ax.TightInset, 0.02);


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