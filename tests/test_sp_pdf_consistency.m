function test_sp_pdf_consistency()
%TEST_SP_PDF_CONSISTENCY Regression for HARTREE PDF Eqs. (1), (5), (6).

    cfg = default_config();
    c = cfg.constants;

    z = [-2; 0.5; 3]*c.nm;
    u = z/(cfg.structure.Lz_nm*c.nm);
    expected_U = cfg.structure.manning_prefactor * ...
        cfg.structure.U0_meV*c.meV .* ...
        (sech(u).^4-sech(u).^2);
    potential_error = max(abs(anharmonic_potential(z,cfg)-expected_U));
    assert(potential_error <= 20*eps(max(abs(expected_U))), ...
        'Manning potential does not reproduce PDF Eq. (1).');

    % Isolate the magnetic term in U_eff=(hbar*kx+eBz)^2/(2m*).
    cfg.fields.B_T = 7;
    cfg.fields.kx_inv_m = 1.7e8;
    dz = 0.25*c.nm;
    H = build_sp_hamiltonian(z,dz,zeros(size(z)),zeros(size(z)),cfg);
    cfg0 = cfg;
    cfg0.fields.B_T = 0;
    H0 = build_sp_hamiltonian(z,dz,zeros(size(z)),zeros(size(z)),cfg0);
    mstar = cfg.material.mstar_rel*c.m0;
    expected_delta = ((c.hbar*cfg.fields.kx_inv_m + ...
        c.e*cfg.fields.B_T*z).^2 - ...
        (c.hbar*cfg.fields.kx_inv_m).^2)/(2*mstar);
    actual_delta = diag(H-H0);
    magnetic_scale = max(abs(expected_delta));
    assert(max(abs(actual_delta-expected_delta)) < ...
        1.0e-12*magnetic_scale, ...
        'Parallel-field term does not reproduce PDF Eq. (5)-(6).');

    % The loop must honor n_states and converge cleanly in the neutral,
    % zero-Hartree limit instead of leaving dEF=Inf forever.
    empty_cfg = apply_numerical_profile(default_config(),'quick');
    empty_cfg.structure.Nz = 180;
    empty_cfg.structure.n_states = 3;
    empty_cfg.doping.Nd_sheet_cm2 = 0;
    empty_cfg.sp.max_iter = 8;
    sp = solve_schrodinger_poisson(empty_cfg);
    assert(sp.converged && sp.iterations == 1, ...
        'Empty-well self-consistent loop did not converge in one iteration.');
    assert(numel(sp.E_J) == 3 && size(sp.Psi,2) == 3, ...
        'SP solver ignored cfg.structure.n_states.');
    assert(all(sp.subband_sheet_m2 == 0) && ...
        max(abs(sp.VH_J)) == 0 && isnan(sp.EF_J), ...
        'Zero-doping Hartree limit is inconsistent.');

    fprintf('  PASS: HARTREE PDF SP equations and loop\n');
end
