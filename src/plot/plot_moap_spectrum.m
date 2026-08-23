function fig = plot_moap_spectrum(spectrum, cfg, tag, mechanism)


    if nargin < 4
        error('plot_moap_spectrum requires a mechanism.');
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

    plot( E, spectrum.mechanism.(mk).total_plot, '-', 'LineWidth', 2.0);

    xlabel('Photon energy \hbar\Omega (meV)');

    if strcmpi(cfg.oap.plot_normalization, 'global_max')

        ylabel('Optical absorption power (normalized)');

    else

        ylabel('Optical absorption power (raw SI-derived scale)');

    end


    title_text = sprintf( '%s spectrum - %s', mechanism_title, tag);

    title(strrep(title_text, '_', '\_'));

    grid on;
    box on;

end