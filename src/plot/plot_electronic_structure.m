function fig = plot_electronic_structure(sp, cfg)
%PLOT_ELECTRONIC_STRUCTURE
% Publication-style electronic-structure plot.
%
% The figure follows the style used in:
%   - Fig. 2 of the reference paper [6]
%   - Fig. 2 of Musa et al.
%
% It shows:
%   1) the self-consistent effective confinement potential U_eff(z),
%   2) the first four eigenenergies,
%   3) the corresponding probability densities |psi_n(z)|^2,
%      vertically shifted to their eigenenergies,
%   4) the Fermi level E_F, when finite.
%
% INPUT
%   sp  : output structure from solve_schrodinger_poisson
%   cfg : project configuration
%
% OUTPUT
%   fig : MATLAB figure handle


    %% ====================================================================
    %  Figure and axes
    % =====================================================================
    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 720 520]);

    ax = axes(fig);

    hold(ax, 'on');


    %% ====================================================================
    %  Data
    % =====================================================================
    z = sp.z_nm(:);

    % IMPORTANT:
    % Eigenenergies were obtained from the Hamiltonian containing
    %
    % U_eff = U_conf + U_H + U_B + U_E
    %
    % therefore the effective potential must be plotted here rather
    % than only the bare confinement potential.
    Ueff = sp.Veff_meV(:);

    E = sp.E_meV(:);

    Psi = sp.Psi;


    %% ====================================================================
    %  Number of displayed states
    %
    %  Fig. 2 in the reference papers displays the lowest four states.
    % =====================================================================
    n_states_plot = min([ ...
        4, ...
        numel(E), ...
        size(Psi, 2)]);

    Eplot = E(1:n_states_plot);


    %% ====================================================================
    %  Publication-style colors
    %
    %  Similar visual hierarchy to the reference figures:
    %       potential       : blue
    %       ground state    : red
    %       first excited   : green
    %       second excited  : black
    %       third excited   : magenta
    % =====================================================================
    potential_color = [0.000 0.300 0.800];

    state_colors = [ ...
        0.850 0.100 0.100; ...   % n = 0
        0.100 0.550 0.100; ...   % n = 1
        0.050 0.050 0.050; ...   % n = 2
        0.750 0.000 0.650  ...   % n = 3
    ];

    state_styles = { ...
        '-', ...
        '--', ...
        '-.', ...
        ':' ...
    };

    state_markers = { ...
        'o', ...
        's', ...
        '^', ...
        'd' ...
    };


    %% ====================================================================
    %  Scale probability densities
    %
    %  We plot
    %
    %       E_n + A |psi_n(z)|^2
    %
    %  as done conventionally in electronic-structure figures.
    %
    %  A is chosen from the energy-level spacing so that the probability
    %  distributions remain visible but do not strongly overlap adjacent
    %  states.
    % =====================================================================

    if n_states_plot >= 2

        dE = diff(sort(Eplot));

        dE = dE(isfinite(dE) & dE > 1.0e-9);

    else

        dE = [];

    end


    if ~isempty(dE)

        % Use a fraction of the characteristic subband spacing.
        density_scale = 0.35 * median(dE);

    else

        potential_span = max(Ueff) - min(Ueff);

        density_scale = 0.08 * potential_span;

    end


    % Numerical safety
    if ~isfinite(density_scale) || density_scale <= 0

        density_scale = 10.0;

    end


    %% ====================================================================
    %  Sparse markers
    %
    %  Do NOT place a marker at every finite-difference grid point.
    %  That would make a dense spectrum look like a thick band.
    % =====================================================================
    n_marker = min(24, numel(z));

    marker_idx = unique(round(linspace(1, numel(z), n_marker)));


    %% ====================================================================
    %  Effective confinement potential
    % =====================================================================
    hPotential = plot( ...
        ax, ...
        z, ...
        Ueff, ...
        '-', ...
        'Color', potential_color, ...
        'LineWidth', 2.0, ...
        'DisplayName', ...
        '$U_{\mathrm{eff}}(z)$');


    % Sparse numerical points on the potential
    plot( ...
        ax, ...
        z(marker_idx), ...
        Ueff(marker_idx), ...
        'o', ...
        'LineStyle', 'none', ...
        'Color', potential_color, ...
        'MarkerSize', 3.4, ...
        'MarkerFaceColor', 'w', ...
        'HandleVisibility', 'off');


    %% ====================================================================
    %  Energy levels + probability densities
    % =====================================================================
    hState = gobjects(n_states_plot, 1);

    y_state_max = zeros(n_states_plot, 1);
    y_state_min = zeros(n_states_plot, 1);


    for i = 1:n_states_plot

        n = i - 1;


        % ---------------------------------------------------------------
        % Probability density
        % ---------------------------------------------------------------
        rho = abs(Psi(:, i)).^2;

        rho_max = max(rho);


        if isfinite(rho_max) && rho_max > 0

            rho = rho ./ rho_max;

        else

            rho(:) = 0;

        end


        % ---------------------------------------------------------------
        % Shift probability density to the eigenenergy
        %
        % y_n(z) = E_n + A |psi_n(z)|^2
        % ---------------------------------------------------------------
        y_state = E(i) + density_scale .* rho;


        y_state_min(i) = min(y_state);
        y_state_max(i) = max(y_state);


        % ---------------------------------------------------------------
        % Horizontal eigenenergy reference
        % ---------------------------------------------------------------
        yline( ...
            ax, ...
            E(i), ...
            ':', ...
            'Color', state_colors(i, :), ...
            'LineWidth', 0.75, ...
            'HandleVisibility', 'off');


        % ---------------------------------------------------------------
        % Probability-density curve
        % ---------------------------------------------------------------
        hState(i) = plot( ...
            ax, ...
            z, ...
            y_state, ...
            state_styles{i}, ...
            'Color', state_colors(i, :), ...
            'LineWidth', 1.55, ...
            'DisplayName', sprintf('$n=%d$', n));


        % ---------------------------------------------------------------
        % Sparse numerical markers
        % ---------------------------------------------------------------
        plot( ...
            ax, ...
            z(marker_idx), ...
            y_state(marker_idx), ...
            state_markers{i}, ...
            'LineStyle', 'none', ...
            'Color', state_colors(i, :), ...
            'MarkerSize', 3.3, ...
            'MarkerFaceColor', 'w', ...
            'HandleVisibility', 'off');

    end


    %% ====================================================================
    %  Fermi level
    % =====================================================================
    hEF = gobjects(0);

    if isfield(sp, 'EF_meV') && ...
            isscalar(sp.EF_meV) && ...
            isfinite(sp.EF_meV)

        hEF = yline( ...
            ax, ...
            sp.EF_meV, ...
            '--', ...
            'Color', [0.85 0.55 0.00], ...
            'LineWidth', 1.40, ...
            'DisplayName', '$E_{\mathrm{F}}$');

    end


    %% ====================================================================
    %  Axis labels
    % =====================================================================
    xlabel( ...
        ax, ...
        '$z$ (nm)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


    ylabel( ...
        ax, ...
        'Energy (meV)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);


    %% ====================================================================
    %  Axis limits
    %
    %  x-axis:
    %  exact numerical domain -- no unnecessary horizontal whitespace.
    %
    %  y-axis:
    %  automatically fits potential, states and Fermi level with only
    %  a small margin.
    % =====================================================================
    xlim(ax, [min(z), max(z)]);


    all_y = [ ...
        Ueff(:); ...
        y_state_min(:); ...
        y_state_max(:); ...
        Eplot(:) ...
    ];


    if ~isempty(hEF)

        all_y = [all_y; sp.EF_meV];

    end


    all_y = all_y(isfinite(all_y));


    if ~isempty(all_y)

        ymin = min(all_y);
        ymax = max(all_y);

        yrange = ymax - ymin;


        if yrange <= 0 || ~isfinite(yrange)

            yrange = max(1.0, abs(ymax));

        end


        % Small article-style margin only
        ypad = 0.045 * yrange;

        ylim(ax, [ ...
            ymin - ypad, ...
            ymax + ypad ...
        ]);

    end


    %% ====================================================================
    %  Legend
    % =====================================================================
    legend_handles = [hPotential; hState(:)];


    if ~isempty(hEF)

        legend_handles = [legend_handles; hEF];

    end


    lgd = legend( ...
        ax, ...
        legend_handles, ...
        'Location', 'best', ...
        'Interpreter', 'latex', ...
        'FontSize', 10, ...
        'Box', 'on');

    lgd.ItemTokenSize = [24 12];


    %% ====================================================================
    %  Axes publication formatting
    %
    %  - inward ticks
    %  - major + minor ticks
    %  - light major grid
    %  - no minor grid
    %  - closed rectangular frame
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

    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.22;

    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';

    ax.LooseInset = max(ax.TightInset, 0.02);

end