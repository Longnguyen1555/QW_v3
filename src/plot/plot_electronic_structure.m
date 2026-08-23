function fig = plot_electronic_structure(sp, cfg)
%PLOT_ELECTRONIC_STRUCTURE Potential, energy levels and probability density.

    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    hold on;

    plot(sp.z_nm, sp.Vconf_meV, '--', 'LineWidth', 1.4);
    plot(sp.z_nm, sp.EF_meV, 'k-', 'LineWidth', 1.8);

    density_scale = 0.10 * (max(sp.Veff_meV)-min(sp.Veff_meV));
    if density_scale <= 0
        density_scale = 10;
    end

    max_states = min(size(sp.Psi,2), 4);
    for i = 1:max_states
        rho = abs(sp.Psi(:,i)).^2;
        rho = rho / max(rho);
        plot(sp.z_nm, sp.E_meV(i)+density_scale*rho, 'LineWidth', 1.5);
        yline(sp.E_meV(i), ':', sprintf('E_%d',i-1), ...
            'LabelHorizontalAlignment','left');
    end

    xlabel('z (nm)');
    ylabel('Energy (meV)');
    title('Self-consistent anharmonic quantum well');
    legend({'U_A(z)','EF(z)','|\psi_i|^2 + E_i'}, ...
           'Location','best');
    grid on; box on;
end
