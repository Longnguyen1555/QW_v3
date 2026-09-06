function fig = plot_linewidth_sweep(linewidth, cfg, mechanism)
%PLOT_LINEWIDTH_SWEEP
% Publication-style FWHM dependence on a sweep parameter.
%
% Visual convention:
%
%   Color        -> photon absorption order
%                   1PA : blue
%                   2PA : red
%                   3PA : magenta
%
%   Marker       -> photon absorption order
%                   1PA : circle
%                   2PA : square
%                   3PA : triangle
%
%   Line style   -> phonon process
%                   emission   : solid
%                   absorption : dashed
%
% This follows the plotting convention of the reference FWHM/HWHM
% calculation code, while adding markers and publication-style axes.


    %% ====================================================================
    %  Input validation
    % =====================================================================
    if nargin < 3
        error('plot_linewidth_sweep requires a mechanism.');
    end


    mk = lower(mechanism);


    switch mk

        case 'optical'
            % Electron-LO phonon linewidth.

        case 'piezoelectric'
            % Electron-piezoelectric phonon linewidth.

        otherwise

            error( ...
                'Unsupported mechanism: %s', ...
                mechanism);

    end


    if ~isfield(linewidth, mk)

        error( ...
            'Linewidth data do not contain mechanism "%s".', ...
            mk);

    end


    %% ====================================================================
    %  Sweep parameter
    % =====================================================================
    x = linewidth.values(:);


    if isempty(x)

        error('Linewidth sweep contains no parameter values.');

    end


    photon_orders = cfg.oap.photon_orders(:).';

    n_orders = numel(photon_orders);


    %% ====================================================================
    %  Figure
    % =====================================================================
    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 700 510]);


    ax = axes(fig);

    hold(ax, 'on');


    %% ====================================================================
    %  Style
    %
    %  The original reference code used approximately:
    %
    %      l = 1 -> red
    %      l = 2 -> blue
    %      l = 3 -> magenta
    %
    %  Here the colors are slightly softened for publication output,
    %  while preserving a clear color identity for each photon order.
    % =====================================================================
    order_colors = [ ...
        0.000 0.300 0.800; ...   % 1PA
        0.850 0.100 0.100; ...   % 2PA
        0.700 0.000 0.650; ...   % 3PA
        0.100 0.550 0.100; ...   % higher order
        0.100 0.600 0.650  ...   % higher order
    ];


    order_markers = { ...
        'o', ...     % 1PA
        's', ...     % 2PA
        '^', ...     % 3PA
        'd', ...
        'v' ...
    };


    %% ====================================================================
    %  Plot storage
    % =====================================================================
    legend_handles = gobjects(0);

    legend_labels = {};

    y_all = [];


    %% ====================================================================
    %  Photon-order resolved FWHM
    % =====================================================================
    for io = 1:n_orders

        order = photon_orders(io);

        ok = sprintf('order_%d', order);


        if ~isfield(linewidth.(mk), ok)

            continue;

        end


        %% ----------------------------------------------------------------
        %  Data
        % -----------------------------------------------------------------
        y_em = ...
            linewidth.(mk).(ok).emission_fwhm(:);

        y_ab = ...
            linewidth.(mk).(ok).absorption_fwhm(:);


        if numel(y_em) ~= numel(x) || ...
           numel(y_ab) ~= numel(x)

            error( ...
                'FWHM array size mismatch for photon order %d.', ...
                order);

        end


        %% ----------------------------------------------------------------
        %  Style for this photon order
        % -----------------------------------------------------------------
        color_idx = ...
            mod(io - 1, size(order_colors, 1)) + 1;

        marker_idx = ...
            mod(io - 1, numel(order_markers)) + 1;


        curve_color = ...
            order_colors(color_idx, :);

        curve_marker = ...
            order_markers{marker_idx};


        %% ================================================================
        %  Phonon emission
        %
        %  Reference code:
        %
        %      Widthpx -> solid line
        %
        %  Filled marker = emission.
        % =================================================================
        valid_em = ...
            isfinite(x) & ...
            isfinite(y_em);


        if any(valid_em)

            h_em = plot( ...
                ax, ...
                x(valid_em), ...
                y_em(valid_em), ...
                '-', ...
                'Color', curve_color, ...
                'LineWidth', 1.60, ...
                'Marker', curve_marker, ...
                'MarkerSize', 5.0, ...
                'MarkerEdgeColor', curve_color, ...
                'MarkerFaceColor', curve_color);


            legend_handles(end+1) = h_em; %#ok<AGROW>

            legend_labels{end+1} = ...
                sprintf( ...
                    '$\\ell=%d$ emission', ...
                    order); %#ok<AGROW>


            y_all = [ ...
                y_all; ...
                y_em(valid_em) ...
            ]; %#ok<AGROW>

        end


        %% ================================================================
        %  Phonon absorption
        %
        %  Reference code:
        %
        %      Widthht -> dashed line
        %
        %  Open marker = absorption.
        % =================================================================
        valid_ab = ...
            isfinite(x) & ...
            isfinite(y_ab);


        if any(valid_ab)

            h_ab = plot( ...
                ax, ...
                x(valid_ab), ...
                y_ab(valid_ab), ...
                '--', ...
                'Color', curve_color, ...
                'LineWidth', 1.60, ...
                'Marker', curve_marker, ...
                'MarkerSize', 5.0, ...
                'MarkerEdgeColor', curve_color, ...
                'MarkerFaceColor', 'w');


            legend_handles(end+1) = h_ab; %#ok<AGROW>

            legend_labels{end+1} = ...
                sprintf( ...
                    '$\\ell=%d$ absorption', ...
                    order); %#ok<AGROW>


            y_all = [ ...
                y_all; ...
                y_ab(valid_ab) ...
            ]; %#ok<AGROW>

        end

    end


    %% ====================================================================
    %  X-axis label
    %
    %  Use physical notation rather than exposing internal config names
    %  such as "B_T" or "Lz_nm".
    % =====================================================================
    xlabel_text = ...
        linewidth_sweep_xlabel( ...
            cfg.sweep.parameter);


    xlabel( ...
        ax, ...
        xlabel_text, ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


    %% ====================================================================
    %  Y-axis label
    % =====================================================================
    ylabel( ...
        ax, ...
        'FWHM (meV)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


    %% ====================================================================
    %  X limits
    %
    %  Sweep plots should not contain unused horizontal space.
    % =====================================================================
    x_finite = x(isfinite(x));


    if ~isempty(x_finite)

        xmin = min(x_finite);

        xmax = max(x_finite);


        if xmax > xmin

            xlim(ax, [xmin xmax]);

        else

            % Safety for a one-point sweep.
            dx = max( ...
                0.05 * abs(xmin), ...
                1.0);

            xlim(ax, [ ...
                xmin - dx, ...
                xmax + dx ...
            ]);

        end

    end


    %% ====================================================================
    %  Y limits
    %
    %  FWHM is a positive physical quantity.
    %
    %  We retain zero as the physical lower baseline and leave only a
    %  small margin above the largest calculated linewidth.
    % =====================================================================
    y_all = y_all(isfinite(y_all));


    if ~isempty(y_all)

        ymax = max(y_all);


        if ymax > 0

            ylim(ax, [ ...
                0, ...
                1.06 * ymax ...
            ]);

        end

    end


    %% ====================================================================
    %  Legend
    %
    %  We use l rather than "1PA/2PA/3PA" because the theoretical
    %  equations and the reference linewidth figures denote the photon
    %  order by l.
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


        lgd.ItemTokenSize = [24 11];

    end


    %% ====================================================================
    %  Axes publication formatting
    %
    %  Requirements discussed for the whole plotting module:
    %
    %     - closed rectangular frame,
    %     - inward ticks,
    %     - minor divisions,
    %     - light major grid,
    %     - no minor grid.
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

    ax.GridAlpha = 0.22;

    ax.XMinorGrid = 'off';

    ax.YMinorGrid = 'off';


    %% ====================================================================
    %  No title
    %
    %  The mechanism and fixed physical parameters should be given in
    %  the figure caption when the plot is used in a paper.
    % =====================================================================


    %% ====================================================================
    %  Compact layout
    % =====================================================================
    ax.LooseInset = ...
        max(ax.TightInset, 0.02);


end



% =========================================================================
%  Local helper: physical x-axis label
% =========================================================================
function label = linewidth_sweep_xlabel(parameter)

    switch lower(parameter)

        case 'b_t'

            label = ...
                'Magnetic field $B$ (T)';


        case 't_k'

            label = ...
                'Temperature $T$ (K)';


        case 'lz_nm'

            label = ...
                'Well width $L_z$ (nm)';


        case 'e_kvcm'

            label = ...
                'Electric field $E$ (kV/cm)';


        case 'nd_sheet_cm2'

            label = ...
                'Sheet doping density $N_D$ (cm$^{-2}$)';


        case 'alpha'

            label = ...
                'Structural parameter $\alpha$';


        case 'u0_mev'

            label = ...
                'Potential depth $U_0$ (meV)';


        otherwise

            % Safe fallback for a newly added sweep parameter.
            label = ...
                strrep(parameter, '_', '\_');

    end

end