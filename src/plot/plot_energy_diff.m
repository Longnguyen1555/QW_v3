function fig = plot_energy_diff(sweep_result, cfg)
%PLOT_ENERGY_DIFF Plot the unnormalized lowest-subband gap E2-E1.

    x = sweep_result.sweep_value(:);
    y = sweep_result.deltaE21_meV(:);

    if isempty(x)
        error('Energy-gap sweep contains no parameter values.');
    end

    if numel(x) ~= numel(y)
        error('Energy-gap sweep values and energy differences must match.');
    end

    valid = isfinite(x) & isfinite(y);

    if ~any(valid)
        error('Energy-gap sweep contains no finite points to plot.');
    end

    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 700 510]);

    ax = axes(fig);

    plot( ...
        ax, ...
        x(valid), ...
        y(valid), ...
        '-o', ...
        'Color', [0.000 0.300 0.800], ...
        'LineWidth', 1.6, ...
        'MarkerSize', 5.0, ...
        'MarkerEdgeColor', [0.000 0.300 0.800], ...
        'MarkerFaceColor', 'w');

    xlabel( ...
        ax, ...
        energy_gap_sweep_xlabel(cfg.sweep.parameter), ...
        'Interpreter', 'latex', ...
        'FontSize', 13);

    ylabel( ...
        ax, ...
        '$\Delta E_{21}=E_2-E_1$ (meV)', ...
        'Interpreter', 'latex', ...
        'FontSize', 13);

    x_finite = x(valid);
    xmin = min(x_finite);
    xmax = max(x_finite);

    if xmax > xmin
        xlim(ax, [xmin xmax]);
    else
        dx = max(0.05 * abs(xmin), 1.0);
        xlim(ax, [xmin - dx, xmax + dx]);
    end

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
    ax.LooseInset = max(ax.TightInset, 0.02);
end


% =========================================================================
% Physical x-axis label for every supported sweep parameter.
% =========================================================================
function label = energy_gap_sweep_xlabel(parameter)

    switch lower(parameter)
        case 'b_t'
            label = 'Magnetic field $B$ (T)';
        case 't_k'
            label = 'Temperature $T$ (K)';
        case 'lz_nm'
            label = 'Well width $L_z$ (nm)';
        case 'e_kvcm'
            label = 'Electric field $E$ (kV/cm)';
        case 'nd_sheet_cm2'
            label = 'Sheet doping density $N_D$ (cm$^{-2}$)';
        case 'manning_prefactor'
            label = 'Manning prefactor $\nu$';
        case 'u0_mev'
            label = 'Potential depth $U_0$ (meV)';
        otherwise
            label = strrep(parameter, '_', '\_');
    end
end
