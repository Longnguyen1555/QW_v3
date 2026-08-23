function fig = plot_electronic_structure(sp, cfg)


    fig = figure('Visible', cfg.output.visible, 'Color', 'w');
    hold on;

    % Confining potential
    hV = plot(sp.z_nm, sp.Vconf_meV, '--', 'LineWidth', 1.4, 'DisplayName', 'V_{conf}(z)');

    % Fermi level: EF is a scalar -> horizontal line
    hEF = yline(sp.EF_meV, 'k-.', 'LineWidth', 1.8, 'DisplayName', 'E_F');

    density_scale = 0.10 * (max(sp.Veff_meV) - min(sp.Veff_meV));

    if ~isfinite(density_scale) || density_scale <= 0
        density_scale = 10;
    end

    max_states = min(size(sp.Psi,2), 4);

    hPsi = gobjects(1,1);

    for i = 1:max_states

        rho = abs(sp.Psi(:,i)).^2;

        rho_max = max(rho);
        if rho_max > 0
            rho = rho / rho_max;
        end

        if i == 1
            hPsi = plot(sp.z_nm, sp.E_meV(i) + density_scale*rho, 'LineWidth', 1.5, 'DisplayName', 'E_i + scaled |\psi_i|^2');
        else
            plot(sp.z_nm, sp.E_meV(i) + density_scale*rho, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end

        yline(sp.E_meV(i), ':', sprintf('E_%d', i-1), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    end

    xlabel('z (nm)');
    ylabel('Energy (meV)');

    title('Self-consistent electronic structure');

    legend([hV, hEF, hPsi], {'V_{conf}(z)', 'E_F', 'E_i'}, 'Location', 'best');

    grid on;
    box on;

end