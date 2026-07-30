function fig = plot_sweep_spectra(cases, labels, cfg)
%PLOT_SWEEP_SPECTRA Overlay total spectra as in the reference article.

    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    hold on;

    for i = 1:numel(cases)
        plot(cases{i}.spectrum.energy_meV, ...
             cases{i}.spectrum.total_plot, 'LineWidth', 1.8);
    end

    xlabel('Photon energy \hbar\Omega (meV)');
    ylabel('Optical absorption power (normalized per case)');
    title(sprintf('MOAP dependence on %s', strrep(cfg.sweep.parameter,'_','\_')));
    legend(labels,'Location','best');
    grid on; box on;
end
