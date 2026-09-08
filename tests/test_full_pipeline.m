function test_full_pipeline()
    cfg = apply_numerical_profile(default_config(),'quick');
    cfg.output.save_figures = false;
    cfg.output.save_csv = false;
    cfg.output.save_mat = false;
    cfg.run.task = 'single';
    cfg.fields.E_kVcm = 0;

    sp = solve_schrodinger_poisson(cfg);
    td = compute_transition_data(sp,cfg);
    s = compute_moap_spectrum(sp,td,cfg);

    assert(sp.converged, 'SP solver did not converge in full-pipeline test.');
    assert(all(isfinite(s.total_raw)), 'Non-finite spectrum.');
    assert(max(s.total_raw)>0, 'Zero spectrum.');
    fprintf('  PASS: full pipeline\n');
end
