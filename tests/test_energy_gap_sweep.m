function test_energy_gap_sweep()
%TEST_ENERGY_GAP_SWEEP Verify direct extraction of E2-E1 from the SP solver.

    cfg = apply_numerical_profile(default_config(), 'quick');
    cfg.run.task = 'energy_gap_sweep';
    cfg.run.verbose = false;
    cfg.sweep.parameter = 'T_K';
    cfg.sweep.values = [80 100];
    cfg.output.save_figures = false;
    cfg.output.save_csv = false;
    cfg.output.save_mat = false;
    cfg.output.visible = 'off';

    cleanup = onCleanup(@() close('all'));

    results = run_qw_moap(cfg);
    gap = results.energy_gap;

    assert(all(isfinite(gap.E1_meV)), 'Non-finite E1 in energy-gap sweep.');
    assert(all(isfinite(gap.E2_meV)), 'Non-finite E2 in energy-gap sweep.');
    assert(all(isfinite(gap.deltaE21_meV)), ...
        'Non-finite DeltaE21 in energy-gap sweep.');
    assert(all(gap.E2_meV >= gap.E1_meV), ...
        'Energy-gap sweep returned E2 < E1.');
    assert(all(gap.deltaE21_meV >= 0), ...
        'Energy-gap sweep returned a negative DeltaE21.');
    assert(all(gap.converged), ...
        'SP solver did not converge in the quick energy-gap sweep test.');
    assert(all(abs(gap.deltaE21_meV - ...
        (gap.E2_meV - gap.E1_meV)) <= eps(max(gap.E2_meV))), ...
        'DeltaE21 was not computed directly as E2-E1.');

    cfg_single = set_sweep_parameter( ...
        cfg, ...
        cfg.sweep.parameter, ...
        cfg.sweep.values(1));
    sp_single = solve_schrodinger_poisson(cfg_single);

    tolerance_meV = 1.0e-8 * max(1.0, max(abs(sp_single.E_meV(1:2))));

    assert(abs(gap.E1_meV(1) - sp_single.E_meV(1)) <= tolerance_meV, ...
        'Swept E1 does not match the same single-case SP calculation.');
    assert(abs(gap.E2_meV(1) - sp_single.E_meV(2)) <= tolerance_meV, ...
        'Swept E2 does not match the same single-case SP calculation.');

    fprintf('  PASS: energy-gap sweep\n');
end
