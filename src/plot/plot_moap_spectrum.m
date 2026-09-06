function fig = plot_moap_spectrum(spectrum, cfg, tag, mechanism)
%PLOT_MOAP_SPECTRUM
% Publication-style total magneto-optical absorption power spectrum.
%
% Style is based mainly on:
%   - Fig. 3, Micro and Nanostructures 198 (2025) 208062
%   - Fig. 4, Physics Letters A 565 (2026) 131162
%   - Fig. 5, Musa & Dakhlaoui, Results in Engineering 28 (2025)
%
% Main plotting rules:
%   1) photon energy is written as hbar*Omega,
%   2) numerical markers are shown along the calculated curve,
%   3) major grid only; minor ticks are enabled,
%   4) closed axes box,
%   5) no title inside the plotting area,
%   6) axes are cropped to the physically relevant spectral region.


    %% ====================================================================
    %  Input validation
    % =====================================================================
    if nargin < 4
        error('plot_moap_spectrum requires a mechanism.');
    end

    mk = lower(mechanism);

    if ~isfield(spectrum.mechanism, mk)
        error('Unknown MOAP mechanism: %s', mechanism);
    end


    switch mk

        case 'optical'

            % Used only for choosing the visual identity of the curve.
            curve_color = [0.35 0.00 0.55];

        case 'piezoelectric'

            curve_color = [0.75 0.10 0.10];

        otherwise

            error('Unsupported mechanism: %s', mechanism);

    end


    %% ====================================================================
    %  Data
    % =====================================================================
    E = spectrum.energy_meV(:);

    P = spectrum.mechanism.(mk).total_plot(:);


    if numel(E) ~= numel(P)
        error('Photon-energy and MOAP arrays must have the same length.');
    end


    %% ====================================================================
    %  Figure
    % =====================================================================
    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 680 500]);


    ax = axes(fig);

    hold(ax, 'on');


    %% ====================================================================
    %  Main spectrum
    %
    %  The reference PLA plots use small numerical points directly on
    %  the calculated spectral curves. We keep the solid curve and
    %  overlay sparse markers so that a dense energy grid does not turn
    %  into a thick band of markers.
    % =====================================================================
    hSpectrum = plot( ...
        ax, ...
        E, ...
        P, ...
        '-', ...
        'Color', curve_color, ...
        'LineWidth', 1.75);


    %% ====================================================================
    %  Numerical markers
    %
    %  Do not mark every point because the current photon-energy grid
    %  can contain thousands of points.
    %
    %  Around 30-40 markers gives a visual appearance close to the
    %  spectra in the Physics Letters A reference.
    % =====================================================================
    n_marker = min(36, numel(E));

    marker_idx = unique( ...
        round(linspace(1, numel(E), n_marker)));


    plot( ...
        ax, ...
        E(marker_idx), ...
        P(marker_idx), ...
        'o', ...
        'LineStyle', 'none', ...
        'Color', curve_color, ...
        'MarkerFaceColor', curve_color, ...
        'MarkerEdgeColor', curve_color, ...
        'MarkerSize', 3.0, ...
        'HandleVisibility', 'off');


    %% ====================================================================
    %  Axis labels
    %
    %  The notation follows the reference OAP/MOAP figures:
    %
    %       photon energy = hbar Omega
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

        % The reference papers report the spectral magnitude in
        % arbitrary units rather than claiming an absolute calibrated
        % experimental power.
        ylabel( ...
            ax, ...
            'MOAP (arb. units)', ...
            'Interpreter', 'latex', ...
            'FontSize', 13);

    end


    %% ====================================================================
    %  X limits
    %
    %  Remove unused photon-energy regions.
    %
    %  Lorentzian broadening generally produces nonzero tails over the
    %  entire numerical domain, so xlim([min(E) max(E)]) would often
    %  leave a large amount of visually useless space.
    %
    %  We identify the active spectrum relative to the largest peak.
    % =====================================================================
    valid = isfinite(E) & isfinite(P);

    E_valid = E(valid);
    P_valid = P(valid);


    if isempty(E_valid)

        warning('No finite MOAP data are available for plotting.');

        return;

    end


    P_abs = abs(P_valid);

    Pmax = max(P_abs);


    if isfinite(Pmax) && Pmax > 0

        % Small enough to retain weak satellite / multiphoton peaks,
        % but large enough to remove negligible Lorentzian tails.
        active_threshold = 1.0e-5 * Pmax;

        active = find(P_abs >= active_threshold);


        if numel(active) >= 2

            xmin_active = E_valid(active(1));
            xmax_active = E_valid(active(end));

            active_width = xmax_active - xmin_active;


            if active_width > 0

                % Small amount of whitespace, similar to journal plots.
                xpad = 0.035 * active_width;

            else

                xpad = 0.02 * ...
                    (max(E_valid) - min(E_valid));

            end


            xmin = max( ...
                min(E_valid), ...
                xmin_active - xpad);

            xmax = min( ...
                max(E_valid), ...
                xmax_active + xpad);


            if xmax > xmin

                xlim(ax, [xmin xmax]);

            else

                xlim(ax, [min(E_valid) max(E_valid)]);

            end

        else

            xlim(ax, [min(E_valid) max(E_valid)]);

        end

    else

        xlim(ax, [min(E_valid) max(E_valid)]);

    end


    %% ====================================================================
    %  Y limits
    %
    %  Optical absorption power should be non-negative.
    %
    %  For normal spectra, therefore start exactly at zero and leave
    %  only a small amount of space above the highest resonance peak.
    % =====================================================================
    Pfinite = P_valid(isfinite(P_valid));


    if ~isempty(Pfinite)

        ymax = max(Pfinite);
        ymin = min(Pfinite);


        if isfinite(ymax) && ymax > 0

            % Numerical calculations can occasionally create extremely
            % small negative round-off values. These should not force
            % the plotting window below zero.
            negative_tolerance = 1.0e-10 * ymax;


            if ymin >= -negative_tolerance

                ylim(ax, [0, 1.055*ymax]);

            else

                yrange = ymax - ymin;

                ypad = 0.04 * yrange;

                ylim(ax, [ ...
                    ymin - ypad, ...
                    ymax + ypad ...
                ]);

            end

        end

    end


    %% ====================================================================
    %  Axes formatting
    %
    %  Combination of:
    %    - PLA Fig. 4: light dashed major grid + markers
    %    - Musa Fig. 5: closed rectangular frame
    %    - user's requested minor divisions
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
    %  Minor grid is intentionally disabled because resonance spectra
    %  become visually crowded when both major and minor grids are shown.
    % =====================================================================
    grid(ax, 'on');

    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.22;

    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';


    %% ====================================================================
    %  No title
    %
    %  The reference journal figures identify the physical case in the
    %  caption / panel legend rather than using an axes title.
    %
    %  Keep the input "tag" in the function interface because it is used
    %  elsewhere in the current project workflow.
    % =====================================================================
    %#ok<NASGU>
    tag = tag;


    %% ====================================================================
    %  Tight publication layout
    % =====================================================================
    ax.LooseInset = max(ax.TightInset, 0.02);


end