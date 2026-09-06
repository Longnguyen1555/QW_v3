function fig = plot_moap_contributions(spectrum, cfg, tag, mechanism)
%PLOT_MOAP_CONTRIBUTIONS
% Publication-style photon-order resolved MOAP contributions.
%
% The plot distinguishes:
%
%   color        -> photon absorption order: 1PA, 2PA, 3PA
%   line style   -> phonon process:
%                   solid   = phonon emission
%                   dashed  = phonon absorption
%   marker       -> photon absorption order
%
% Plotting principles follow the MOAP/OAP figures in the reference papers:
%   - photon energy written as hbar*Omega,
%   - closed rectangular axes,
%   - major grid + minor ticks,
%   - numerical points shown on calculated curves,
%   - no title inside the axes,
%   - compact axis limits around the physically relevant spectrum.


    %% ====================================================================
    %  Input validation
    % =====================================================================
    if nargin < 4
        error('plot_moap_contributions requires a mechanism.');
    end

    mk = lower(mechanism);

    if ~isfield(spectrum.mechanism, mk)
        error('Unknown MOAP mechanism: %s', mechanism);
    end


    switch mk

        case 'optical'
            % Electron-LO phonon

        case 'piezoelectric'
            % Electron-piezoelectric acoustic phonon

        otherwise
            error('Unsupported mechanism: %s', mechanism);

    end


    %% ====================================================================
    %  Data
    % =====================================================================
    E = spectrum.energy_meV(:);

    photon_orders = cfg.oap.photon_orders(:).';

    n_orders = numel(photon_orders);


    %% ====================================================================
    %  Figure
    % =====================================================================
    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 720 520]);


    ax = axes(fig);

    hold(ax, 'on');


    %% ====================================================================
    %  Plot style
    %
    %  One color for each photon order.
    %
    %  Emission:
    %       solid line + filled marker
    %
    %  Absorption:
    %       dashed line + open marker
    % =====================================================================

    order_colors = [ ...
        0.000 0.300 0.800; ...    % 1PA - blue
        0.850 0.100 0.100; ...    % 2PA - red
        0.100 0.550 0.100; ...    % 3PA - green
        0.650 0.100 0.650; ...    % higher order if required
        0.100 0.600 0.650  ...    % higher order if required
    ];


    order_markers = { ...
        'o', ...     % 1PA
        's', ...     % 2PA
        '^', ...     % 3PA
        'd', ...
        'v' ...
    };


    %% ====================================================================
    %  Storage for axis limits
    %
    %  The active x-range is calculated from EACH individual contribution,
    %  not from the strongest 1PA curve only.
    %
    %  This is important because 2PA and especially 3PA may be much
    %  weaker than 1PA but still physically relevant.
    % =====================================================================
    x_active_all = [];

    y_all = [];

    legend_handles = gobjects(0);

    labels = {};


    %% ====================================================================
    %  Photon-order contributions
    % =====================================================================
    for io = 1:n_orders

        order = photon_orders(io);

        ok = sprintf('order_%d', order);


        if ~isfield(spectrum.mechanism.(mk), ok)
            continue;
        end


        P_em = ...
            spectrum.mechanism.(mk).(ok).emission_plot(:);

        P_ab = ...
            spectrum.mechanism.(mk).(ok).absorption_plot(:);


        if numel(P_em) ~= numel(E) || ...
           numel(P_ab) ~= numel(E)

            error( ...
                'MOAP contribution size mismatch for photon order %d.', ...
                order);

        end


        %% ----------------------------------------------------------------
        %  Style for this order
        % -----------------------------------------------------------------
        color_idx = mod(io - 1, size(order_colors, 1)) + 1;

        marker_idx_style = ...
            mod(io - 1, numel(order_markers)) + 1;


        curve_color = ...
            order_colors(color_idx, :);

        curve_marker = ...
            order_markers{marker_idx_style};


        %% ----------------------------------------------------------------
        %  Phonon emission
        % -----------------------------------------------------------------
        h_em = plot( ...
            ax, ...
            E, ...
            P_em, ...
            '-', ...
            'Color', curve_color, ...
            'LineWidth', 1.65);


        legend_handles(end+1) = h_em; %#ok<AGROW>

        labels{end+1} = sprintf( ...
            '$%d\\mathrm{PA}$ emission', ...
            order); %#ok<AGROW>


        %% ----------------------------------------------------------------
        %  Markers for phonon emission
        %
        %  Markers are distributed over the physically active part of
        %  the individual curve rather than over the whole energy domain.
        % -----------------------------------------------------------------
        Pem_abs = abs(P_em);

        Pem_max = max(Pem_abs(isfinite(Pem_abs)));


        if ~isempty(Pem_max) && ...
                isfinite(Pem_max) && ...
                Pem_max > 0

            threshold_em = 1.0e-4 * Pem_max;

            idx_em_active = find( ...
                isfinite(P_em) & ...
                Pem_abs >= threshold_em);


            if numel(idx_em_active) >= 2

                n_marker_em = ...
                    min(18, numel(idx_em_active));

                marker_pos = unique(round(linspace( ...
                    1, ...
                    numel(idx_em_active), ...
                    n_marker_em)));


                marker_idx_em = ...
                    idx_em_active(marker_pos);


                plot( ...
                    ax, ...
                    E(marker_idx_em), ...
                    P_em(marker_idx_em), ...
                    curve_marker, ...
                    'LineStyle', 'none', ...
                    'Color', curve_color, ...
                    'MarkerEdgeColor', curve_color, ...
                    'MarkerFaceColor', curve_color, ...
                    'MarkerSize', 3.4, ...
                    'HandleVisibility', 'off');

            end


            x_active_all = [ ...
                x_active_all; ...
                E(idx_em_active) ...
            ]; %#ok<AGROW>

        end


        %% ----------------------------------------------------------------
        %  Phonon absorption
        % -----------------------------------------------------------------
        h_ab = plot( ...
            ax, ...
            E, ...
            P_ab, ...
            '--', ...
            'Color', curve_color, ...
            'LineWidth', 1.65);


        legend_handles(end+1) = h_ab; %#ok<AGROW>

        labels{end+1} = sprintf( ...
            '$%d\\mathrm{PA}$ absorption', ...
            order); %#ok<AGROW>


        %% ----------------------------------------------------------------
        %  Markers for phonon absorption
        %
        %  Open markers distinguish absorption from emission while
        %  retaining the same photon-order marker shape.
        % -----------------------------------------------------------------
        Pab_abs = abs(P_ab);

        Pab_max = max(Pab_abs(isfinite(Pab_abs)));


        if ~isempty(Pab_max) && ...
                isfinite(Pab_max) && ...
                Pab_max > 0

            threshold_ab = 1.0e-4 * Pab_max;

            idx_ab_active = find( ...
                isfinite(P_ab) & ...
                Pab_abs >= threshold_ab);


            if numel(idx_ab_active) >= 2

                n_marker_ab = ...
                    min(18, numel(idx_ab_active));

                marker_pos = unique(round(linspace( ...
                    1, ...
                    numel(idx_ab_active), ...
                    n_marker_ab)));


                marker_idx_ab = ...
                    idx_ab_active(marker_pos);


                plot( ...
                    ax, ...
                    E(marker_idx_ab), ...
                    P_ab(marker_idx_ab), ...
                    curve_marker, ...
                    'LineStyle', 'none', ...
                    'Color', curve_color, ...
                    'MarkerEdgeColor', curve_color, ...
                    'MarkerFaceColor', 'w', ...
                    'MarkerSize', 3.4, ...
                    'HandleVisibility', 'off');

            end


            x_active_all = [ ...
                x_active_all; ...
                E(idx_ab_active) ...
            ]; %#ok<AGROW>

        end


        %% ----------------------------------------------------------------
        %  Y data for automatic axis limit
        % -----------------------------------------------------------------
        y_all = [ ...
            y_all; ...
            P_em(isfinite(P_em)); ...
            P_ab(isfinite(P_ab)) ...
        ]; %#ok<AGROW>

    end


    %% ====================================================================
    %  Labels
    % =====================================================================
    xlabel( ...
        ax, ...
        'Photon energy $\hbar\Omega$ (meV)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


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
    %  Use the union of the active regions of ALL six contributions.
    %
    %  This avoids a common plotting error:
    %  cropping according to the strongest 1PA contribution can remove
    %  the weaker 2PA / 3PA resonance peaks at lower photon energies.
    % =====================================================================
    x_active_all = ...
        x_active_all(isfinite(x_active_all));


    if numel(x_active_all) >= 2

        xmin_active = min(x_active_all);

        xmax_active = max(x_active_all);

        active_width = ...
            xmax_active - xmin_active;


        if active_width > 0

            xpad = 0.035 * active_width;

        else

            xpad = 0.02 * ...
                (max(E) - min(E));

        end


        xmin = max( ...
            min(E), ...
            xmin_active - xpad);

        xmax = min( ...
            max(E), ...
            xmax_active + xpad);


        if xmax > xmin

            xlim(ax, [xmin xmax]);

        else

            xlim(ax, [min(E) max(E)]);

        end

    else

        xlim(ax, [min(E) max(E)]);

    end


    %% ====================================================================
    %  Y-axis limits
    %
    %  MOAP contributions should be non-negative.
    % =====================================================================
    y_all = y_all(isfinite(y_all));


    if ~isempty(y_all)

        ymax = max(y_all);

        ymin = min(y_all);


        if isfinite(ymax) && ymax > 0

            negative_tolerance = ...
                1.0e-10 * ymax;


            if ymin >= -negative_tolerance

                ylim(ax, [ ...
                    0, ...
                    1.06 * ymax ...
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
    %  The paired order is intentional:
    %
    %       1PA emission
    %       1PA absorption
    %       2PA emission
    %       2PA absorption
    %       3PA emission
    %       3PA absorption
    %
    %  With several contributions, two columns produce a much more
    %  compact journal-style legend.
    % =====================================================================
    if ~isempty(legend_handles)

        lgd = legend( ...
            ax, ...
            legend_handles, ...
            labels, ...
            'Location', 'best', ...
            'Interpreter', 'latex', ...
            'FontSize', 9.5, ...
            'Box', 'on');


        if numel(labels) >= 4

            lgd.NumColumns = 2;

        end


        lgd.ItemTokenSize = [22 10];

    end


    %% ====================================================================
    %  Axes formatting
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


    %% ====================================================================
    %  Grid
    %
    %  Major grid only.
    %
    %  Minor ticks are useful for reading resonance positions, while a
    %  minor grid would make six overlapping spectra unnecessarily busy.
    % =====================================================================
    grid(ax, 'on');

    ax.GridLineStyle = '--';

    ax.GridAlpha = 0.22;

    ax.XMinorGrid = 'off';

    ax.YMinorGrid = 'off';


    %% ====================================================================
    %  No title
    %
    %  Journal figures normally put the mechanism / parameter information
    %  in the caption rather than inside the plotting area.
    %
    %  Keep "tag" in the interface because run_single_case currently
    %  passes it to this function.
    % =====================================================================
    %#ok<NASGU>
    tag = tag;


    %% ====================================================================
    %  Compact layout
    % =====================================================================
    ax.LooseInset = max(ax.TightInset, 0.02);


end