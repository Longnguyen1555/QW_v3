function test_physics_configuration()
%TEST_PHYSICS_CONFIGURATION Focused regression tests for source controls.

    cfg = default_config();
    c = cfg.constants;

    z = (-2:0.01:2).' * c.nm;
    U1 = anharmonic_potential(z, cfg);
    cfg.structure.manning_prefactor = 4.0;
    U2 = anharmonic_potential(z, cfg);
    assert(max(abs(U1-U2)) > 0, ...
        'Manning prefactor does not change the confinement potential.');

    alpha_failed = false;
    try
        set_sweep_parameter(cfg, 'alpha', 0.3);
    catch ME
        alpha_failed = strcmp(ME.identifier, 'QW:UnsupportedAlphaSweep');
    end
    assert(alpha_failed, 'Unsupported alpha sweep did not fail clearly.');

    cfg_states = apply_numerical_profile(default_config(), 'quick');
    cfg_states.structure.n_states = 3;
    cfg_states.fields.E_kVcm = 0;
    sp_states = solve_schrodinger_poisson(cfg_states);
    assert(numel(sp_states.E_J) == 3, ...
        'SP solver did not honor cfg.structure.n_states.');

    omega1 = 20*c.meV/c.hbar;
    omega2 = 40*c.meV/c.hbar;
    mstar = cfg.material.mstar_rel*c.m0;
    a01 = c.e*(cfg.laser.E0_kVcm*1e5)/(mstar*omega1^2);
    a02 = c.e*(cfg.laser.E0_kVcm*1e5)/(mstar*omega2^2);
    assert(abs(a01/a02 - 4) < 1e-12, 'Dynamic a0 does not scale as Omega^-2.');

    m0 = profile_fwhm(0:0.1:1, zeros(1,11));
    assert(isnan(m0.peak_x) && isnan(m0.fwhm) && isnan(m0.hwhm), ...
        'An all-zero spectrum must have NaN metrics.');
    md = profile_fwhm(0:0.1:0.2, [0 1 0]);
    assert(isnan(md.fwhm), 'A one-bin delta-like peak must be unresolved.');

    fprintf('  PASS: physics configuration\n');
end
