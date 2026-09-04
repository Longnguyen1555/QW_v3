function fig = plot_sweep_spectra(cases, labels, cfg, mechanism)

    if nargin < 4
        error('plot_sweep_spectra requires a mechanism.');
    end

    mk = lower(mechanism);

    switch mk
        case 'optical'
            mechanism_title = 'Electron-LO phonon';

        case 'piezoelectric'
            mechanism_title = 'Electron-piezoelectric phonon';

        otherwise
            error('Unsupported mechanism: %s', mechanism);
    end

    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w');

    hold on;

    % ================================================================
    % Plot every sweep case
    % ================================================================
    for i = 1:numel(cases)

        spectrum = cases{i}.spectrum;

        if ~isfield(spectrum.mechanism, mk)
            error( ...
                'Mechanism "%s" not found in sweep case %d.', ...
                mk, i);
        end

        plot( ...
            spectrum.energy_meV, ...
            spectrum.mechanism.(mk).total_plot, ...
            'LineWidth', 1.8);

    end

    % ================================================================
    % Axis
    % ================================================================
    xlabel('Photon energy (meV)');

    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        ylabel('Optical absorption power (normalized)');

    else

        ylabel('Optical absorption power (raw SI-derived scale)');

    end

    % ================================================================
    % Title
    % ================================================================
    title_text = sprintf( ...
        '%s dependence on %s', ...
        mechanism_title, ...
        strrep(cfg.sweep.parameter, '_', '\_'));

    title(title_text);

    % ================================================================
    % Legend
    % ================================================================
    legend(labels, 'Location', 'best');

    grid on;
    box on;

end