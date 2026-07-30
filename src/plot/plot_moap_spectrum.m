function fig = plot_moap_spectrum(spectrum, cfg, tag)
%PLOT_MOAP_SPECTRUM Total and mechanism-resolved spectra.

    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    E = spectrum.energy_meV;

    plot(E, spectrum.total_plot, 'k-', 'LineWidth', 2.2); hold on;
    plot(E, spectrum.mechanism.optical.total_plot, '-', 'LineWidth', 1.8);
    plot(E, spectrum.mechanism.piezoelectric.total_plot, '--', 'LineWidth', 1.8);

    xlabel('Photon energy \hbar\Omega (meV)');
    if strcmpi(cfg.oap.plot_normalization,'global_max')
        ylabel('Optical absorption power (normalized)');
    else
        ylabel('Optical absorption power (raw SI-derived scale)');
    end
    title(strrep(sprintf('MOAP spectrum - %s',tag),'_','\_'));
    legend({'Total','Optical phonon','Piezoelectric phonon'}, ...
           'Location','best');
    grid on; box on;
end
