function test_wavefunction_normalization()
    cfg = apply_numerical_profile(default_config(),'quick');
    cfg.sp.max_iter = 120;
    sp = solve_schrodinger_poisson(cfg);
    norms = trapz(sp.z_m, abs(sp.Psi).^2, 1);
    assert(max(abs(norms-1)) < 2e-3, ...
        'Wavefunction normalization test failed.');
    fprintf('  PASS: wavefunction normalization\n');
end
