function fig = plot_linewidth_sweep(linewidth, cfg, mechanism)

    if nargin < 3
        error('plot_linewidth_sweep requires a mechanism.');
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


    if ~isfield(linewidth, mk)
        error( ...
            'Linewidth data do not contain mechanism "%s".', ...
            mk);
    end


    fig = figure( ...
        'Visible', cfg.output.visible, ...
        'Color', 'w');

    hold on;


    % ================================================================
    % Sweep parameter
    % ================================================================
    x = linewidth.values;

    labels = {};


    % ================================================================
    % Photon-order resolved linewidth
    % ================================================================
    for order = cfg.oap.photon_orders

        ok = sprintf('order_%d', order);


        % ------------------------------------------------------------
        % Phonon emission
        % ------------------------------------------------------------
        plot( ...
            x, ...
            linewidth.(mk).(ok).emission_fwhm, ...
            '-o', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4);

        labels{end+1} = sprintf( ...
            '%dPA emission', order); %#ok<AGROW>


        % ------------------------------------------------------------
        % Phonon absorption
        % ------------------------------------------------------------
        plot( ...
            x, ...
            linewidth.(mk).(ok).absorption_fwhm, ...
            '--s', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4);

        labels{end+1} = sprintf( ...
            '%dPA absorption', order); %#ok<AGROW>

    end


    % ================================================================
    % Axis
    % ================================================================
    xlabel(strrep(cfg.sweep.parameter, '_', '\_'));

    ylabel('FWHM (meV)');


    % ================================================================
    % Title
    % ================================================================
    title(sprintf( ...
        '%s linewidth', ...
        mechanism_title));


    % ================================================================
    % Legend
    % ================================================================
    legend(labels, 'Location', 'best');


    grid on;
    box on;

end