function sweep_result = run_energy_gap_sweep(cfg)
%RUN_ENERGY_GAP_SWEEP Sweep E2-E1 using the existing self-consistent solver.

    values = cfg.sweep.values(:);

    if isempty(values)
        error('Energy-gap sweep requires at least one sweep value.');
    end

    if ~isnumeric(values) || any(~isfinite(values))
        error('Energy-gap sweep values must be finite numeric values.');
    end

    n_cases = numel(values);

    E1_meV = nan(n_cases, 1);
    E2_meV = nan(n_cases, 1);
    deltaE21_meV = nan(n_cases, 1);
    converged = false(n_cases, 1);
    iterations = nan(n_cases, 1);
    VH_error_meV = nan(n_cases, 1);
    EF_error_meV = nan(n_cases, 1);

    for i = 1:n_cases
        % Start from the original configuration at every point and change
        % exactly the selected sweep parameter through the shared setter.
        cfg_i = set_sweep_parameter( ...
            cfg, ...
            cfg.sweep.parameter, ...
            values(i));

        sp = solve_schrodinger_poisson(cfg_i);
        energies_meV = sp.E_meV(:);

        if numel(energies_meV) < 2 || ...
                any(~isfinite(energies_meV(1:2)))
            error( ...
                ['The Schrodinger-Poisson solver did not return two ' ...
                 'finite eigenenergies at sweep point %d.'], ...
                i);
        end

        E1_meV(i) = energies_meV(1);
        E2_meV(i) = energies_meV(2);
        deltaE21_meV(i) = E2_meV(i) - E1_meV(i);
        converged(i) = sp.converged;
        iterations(i) = sp.iterations;
        VH_error_meV(i) = sp.VH_error_meV;
        EF_error_meV(i) = sp.EF_error_meV;

        if cfg.run.verbose
            fprintf( ...
                ['[ENERGY GAP] %s = %.6g: ' ...
                 'E1 = %.8f meV, E2 = %.8f meV, ' ...
                 'DeltaE21 = %.8f meV, converged = %d\n'], ...
                cfg.sweep.parameter, ...
                values(i), ...
                E1_meV(i), ...
                E2_meV(i), ...
                deltaE21_meV(i), ...
                converged(i));
        end
    end

    sweep_result = struct();
    sweep_result.parameter = cfg.sweep.parameter;
    sweep_result.values = values;
    sweep_result.sweep_value = values;
    sweep_result.E1_meV = E1_meV;
    sweep_result.E2_meV = E2_meV;
    sweep_result.deltaE21_meV = deltaE21_meV;
    sweep_result.converged = converged;
    sweep_result.iterations = iterations;
    sweep_result.VH_error_meV = VH_error_meV;
    sweep_result.EF_error_meV = EF_error_meV;

    fig = plot_energy_diff(sweep_result, cfg);

    save_sweep_outputs( ...
        sweep_result, ...
        cfg, ...
        'energy_gap_sweep', ...
        {fig}, ...
        {'energy_gap_sweep'});
end
