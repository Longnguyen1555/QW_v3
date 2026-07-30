function fig = plot_linewidth_sweep(linewidth, cfg)
%PLOT_LINEWIDTH_SWEEP FWHM curves analogous to the article-[6] scripts.

    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    x = linewidth.values;
    mech_names = {'optical','piezoelectric'};
    titles = {'Optical phonon','Piezoelectric phonon'};

    for im = 1:2
        subplot(1,2,im);
        hold on;
        labels = {};
        mk = mech_names{im};

        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            plot(x, linewidth.(mk).(ok).emission_fwhm, ...
                '-o','LineWidth',1.5,'MarkerSize',4);
            labels{end+1} = sprintf('%dPA emission',order); %#ok<AGROW>
            plot(x, linewidth.(mk).(ok).absorption_fwhm, ...
                '--s','LineWidth',1.5,'MarkerSize',4);
            labels{end+1} = sprintf('%dPA absorption',order); %#ok<AGROW>
        end

        xlabel(strrep(cfg.sweep.parameter,'_','\_'));
        ylabel('FWHM (meV)');
        title(titles{im});
        legend(labels,'Location','best');
        grid on; box on;
    end
end
