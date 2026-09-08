function test_full_pipeline()
    cfg = apply_numerical_profile(default_config(),'quick');
    cfg.output.save_figures = false;
    cfg.output.save_csv = false;
    cfg.output.save_mat = false;
    cfg.run.task = 'single';
    % The default static-field placeholder is intentionally overridden by
    % a finite physical fixture.  With an exact Dirac root, a transition
    % outside the scanned photon-energy window correctly gives zero rather
    % than the artificial Lorentzian tail used by the old implementation.
    cfg.fields.E_kVcm = 0;
    cfg.structure.Lz_nm = 5;
    cfg.structure.U0_meV = 220;

    sp = solve_schrodinger_poisson(cfg);
    td = compute_transition_data(sp,cfg);
    s = compute_moap_spectrum(sp,td,cfg);

    assert(sp.converged, 'SP solver did not converge in full-pipeline test.');
    assert(all(isfinite(s.total_raw)), 'Non-finite spectrum.');
    assert(max(s.total_raw)>0, 'Zero spectrum.');
    assert(all(isfinite(s.mechanism.optical.total_raw)), ...
        'Non-finite optical spectrum.');
    assert(all(isfinite(s.mechanism.piezoelectric.total_raw)), ...
        'Non-finite piezoelectric spectrum.');
    assert(max(s.mechanism.optical.total_raw) > 0, ...
        'Zero optical spectrum.');
    assert(max(s.mechanism.piezoelectric.total_raw) > 0, ...
        'Zero piezoelectric spectrum.');
    assert(s.meta.gamma_ignored, ...
        'direct_q_integral must not use phenomenological gamma.');
    fprintf('  PASS: full pipeline\n');
end
