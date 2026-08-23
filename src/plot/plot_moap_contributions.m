function fig = plot_moap_contributions(spectrum, cfg, tag, mechanism)

    if nargin < 4
        error('plot_moap_contributions requires a mechanism.');
    end

    mk = lower(mechanism);

    if ~isfield(spectrum.mechanism, mk)
        error('Unknown MOAP mechanism: %s', mechanism);
    end


    switch mk
        case 'optical'
            mechanism_title = 'Electron-LO phonon';

        case 'piezoelectric'
            mechanism_title = 'Electron-piezoelectric phonon';

        otherwise
            error('Unsupported mechanism: %s', mechanism);
    end


    fig = figure('Visible', cfg.output.visible, 'Color', 'w');

    hold on;

    E = spectrum.energy_meV;

    labels = {};


    for order = cfg.oap.photon_orders

        ok = sprintf('order_%d', order);

        plot(E, spectrum.mechanism.(mk).(ok).emission_plot, '-', 'LineWidth', 1.6);

        labels{end+1} = sprintf('%dPA emission', order); %#ok<AGROW>

        plot(E, spectrum.mechanism.(mk).(ok).absorption_plot, '--', 'LineWidth', 1.6);

        labels{end+1} = sprintf('%dPA absorption', order); %#ok<AGROW>

    end

    xlabel('Photon energy (meV)');


    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        ylabel('Optical absorption power (normalized)');

    else

        ylabel('Optical absorption power (raw SI-derived scale)');

    end


    title_text = sprintf('%s contributions - %s', mechanism_title, tag);

    title(strrep(title_text, '_', '\_'));

    legend(labels, 'Location', 'best');


    grid on;
    box on;

end