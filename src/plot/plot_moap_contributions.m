function fig = plot_moap_contributions(spectrum, cfg, tag)
%PLOT_MOAP_CONTRIBUTIONS Photon-order and phonon-process contributions.

    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    E = spectrum.energy_meV;

    nrows = 1;
    ncols = 2;
    mech_names = {'optical','piezoelectric'};
    mech_titles = {'Electron-LO phonon','Electron-piezoelectric phonon'};

    for im = 1:2
        subplot(nrows,ncols,im);
        hold on;
        mk = mech_names{im};
        labels = {};

        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            plot(E, spectrum.mechanism.(mk).(ok).emission_plot, ...
                 '-', 'LineWidth', 1.4);
            labels{end+1} = sprintf('%dPA emission',order); %#ok<AGROW>
            plot(E, spectrum.mechanism.(mk).(ok).absorption_plot, ...
                 '--', 'LineWidth', 1.4);
            labels{end+1} = sprintf('%dPA absorption',order); %#ok<AGROW>
        end

        xlabel('Photon energy \hbar\Omega (meV)');
        ylabel('Normalized OAP');
        title(mech_titles{im});
        legend(labels,'Location','best');
        grid on; box on;
    end

    if exist('sgtitle','file')
        sgtitle(strrep(sprintf('Resolved contributions - %s',tag),'_','\_'));
    end
end
