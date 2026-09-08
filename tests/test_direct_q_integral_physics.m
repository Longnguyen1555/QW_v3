function test_direct_q_integral_physics()
%TEST_DIRECT_Q_INTEGRAL_PHYSICS Physics tests for Eq. (2)/(18) direct path.

    cfg = default_config();
    cfg.oap.plot_normalization = 'none';
    cfg.oap.photon_energy_meV = 10:2:250;
    cfg.oap.qperp_max_inv_nm = 1.2;
    cfg.oap.qz_max_inv_nm = 2.0;
    cfg.oap.Nqperp = 48;
    cfg.oap.Nqz = 36;

    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    q = 0.20 / c.nm;
    Einitial = 0;
    deltaE = 10*c.meV;
    hw0 = cfg.material.LO_phonon_meV*c.meV;

    % TEST 1: root domain and exact Jacobian.
    [w_forbidden, ~, B_forbidden] = direct_k_delta_weight( ...
        Einitial, 20*c.meV, deltaE, 20*c.meV, 1, +1, hw0, q, cfg);
    assert(B_forbidden > 0 && w_forbidden == 0, ...
        'B>0 must have exactly zero direct k contribution.');

    [w_allowed, kstar, B_allowed, fstar] = direct_k_delta_weight( ...
        Einitial, 20*c.meV, deltaE, 80*c.meV, 1, +1, hw0, q, cfg);
    expected = mstar*kstar/(c.hbar^2*q)*fstar;
    assert(B_allowed < 0 && kstar > 0, ...
        'B<0 must produce a positive kstar.');
    assert(abs(w_allowed-expected) <= 20*eps(max(expected,1)), ...
        'Exact delta-root Jacobian does not match its analytic value.');

    [w_qzero, k_qzero] = direct_k_delta_weight( ...
        Einitial, 20*c.meV, deltaE, 80*c.meV, 1, +1, hw0, 0, cfg);
    assert(w_qzero == 0 && k_qzero == 0, ...
        'q_perp=0 endpoint handling must be finite and zero.');

    % TEST 2: actual Fermi-Dirac dependence on EF.
    w_low_EF = direct_k_delta_weight(Einitial, -20*c.meV, deltaE, ...
        80*c.meV, 1, +1, hw0, q, cfg);
    w_high_EF = direct_k_delta_weight(Einitial, 30*c.meV, deltaE, ...
        80*c.meV, 1, +1, hw0, q, cfg);
    assert(w_high_EF > w_low_EF, ...
        'Increasing EF must increase the occupied initial-state weight.');

    % TEST 3: recoil is explicitly present in B(q_perp).
    q_pair = [0.10; 0.20] / c.nm;
    [~, ~, B_pair] = direct_k_delta_weight(Einitial, 20*c.meV, ...
        deltaE, 80*c.meV, 1, +1, hw0, q_pair, cfg);
    expected_recoil_change = c.hbar^2*(q_pair(2)^2-q_pair(1)^2)/(2*mstar);
    assert(abs((B_pair(2)-B_pair(1))-expected_recoil_change) <= ...
           50*eps(max(abs(expected_recoil_change), realmin)), ...
        'B(q_perp) is missing or corrupting the recoil term.');

    [qd, ~] = debye_wavevector(struct(), cfg);

    % TEST 4: optical Bose process convention.
    optical = direct_phonon_terms('optical', q_pair, qd, cfg);
    N0 = 1/expm1(hw0/(c.kB*cfg.temperature_K));
    assert(max(abs(optical.hw_J-hw0)) == 0, ...
        'Optical phonon energy must be constant.');
    assert(max(abs(optical.absorption_C2V-optical.C2V*N0)) <= ...
           20*eps(max(max(optical.absorption_C2V), realmin)), ...
        'Optical absorption must use N0.');
    assert(max(abs(optical.emission_C2V-optical.C2V*(N0+1))) <= ...
           20*eps(max(max(optical.emission_C2V), realmin)), ...
        'Optical emission must use N0+1.');

    % TEST 5: physical acoustic piezo dispersion and q-dependent Bose factor.
    q_ac = [0; 0.10; 0.30] / c.nm;
    piezo = direct_phonon_terms('piezoelectric', q_ac, qd, cfg);
    assert(max(abs(piezo.hw_J-c.hbar*cfg.material.sound_speed_mps*q_ac)) == 0, ...
        'Piezoelectric phonons must use hw_q=hbar*c_s*q_perp.');
    assert(isinf(piezo.Nph(1)) && piezo.Nph(2) > piezo.Nph(3), ...
        'Acoustic Bose population must diverge at zero and vary with q.');
    assert(piezo.absorption_C2V(1) == 0 && piezo.emission_C2V(1) == 0, ...
        'The complete screened acoustic q->0 endpoint limit must be zero.');

    % TEST 6: direct path supports only ell=1,2 and positive dynamic Eph.
    validate_config(cfg);
    cfg_bad_order = cfg;
    cfg_bad_order.oap.photon_orders = [1 2 3];
    assert_throws(@() validate_config(cfg_bad_order), 'QW:DirectPhotonOrder');
    cfg_bad_energy = cfg;
    cfg_bad_energy.oap.photon_energy_meV(1) = 0;
    assert_throws(@() validate_config(cfg_bad_energy), 'QW:DirectPhotonEnergy');

    [sp, td] = synthetic_case(cfg);
    s = compute_moap_direct(sp, td, cfg);

    % TEST 7: direct line shape is independent of gamma and population_mode.
    cfg_gamma = cfg;
    cfg_gamma.oap.gamma_optical_meV = 1234;
    cfg_gamma.oap.gamma_piezo_meV = 987;
    cfg_gamma.oap.population_mode = 'normalized';
    s_gamma = compute_moap_direct(sp, td, cfg_gamma);
    assert(isequal(s.total_raw, s_gamma.total_raw), ...
        'gamma/population_mode must not affect direct_q_integral.');
    cfg_block = cfg;
    cfg_block.oap.direct.energy_block_size = 7;
    s_block = compute_moap_direct(sp, td, cfg_block);
    assert(isequal(s.total_raw, s_block.total_raw), ...
        'Energy blocking must not change the direct spectrum.');

    % TEST 8: q_z affects only the confinement form-factor integral.
    td_scaled = td;
    td_scaled.I2 = 3*td.I2;
    s_scaled = compute_moap_direct(sp, td_scaled, cfg);
    scale_error = norm(s_scaled.total_raw-3*s.total_raw) / ...
                  max(norm(3*s.total_raw), realmin);
    assert(scale_error < 1e-12, ...
        'q_z form factor must enter linearly and separately from q_perp coupling.');
    assert(strcmp(s.meta.q_coupling_approximation, 'q_equals_qperp'), ...
        'Direct metadata must record q=q_perp coupling approximation.');

    % TEST 9: finite nonzero outputs and dynamic a0.
    assert(all(isfinite(s.total_raw)) && max(s.total_raw) > 0, ...
        'Direct optical/piezo spectrum must be finite and nonzero.');
    expected_a0 = c.e*(cfg.laser.E0_kVcm*1e5) ./ ...
        (mstar*(s.energy_J/c.hbar).^2);
    assert(max(abs(s.meta.a0_m-expected_a0)) <= ...
           20*eps(max(expected_a0)), ...
        'direct_q_integral must use dynamic a0=eE0/(m*Omega^2).');
    assert(s.meta.a0_m(1) > s.meta.a0_m(end), ...
        'Dynamic a0 must vary with photon energy.');

    % TEST 10: representative two-grid convergence.
    cfg_fine = cfg;
    cfg_fine.oap.Nqperp = 96;
    cfg_fine.oap.Nqz = 72;
    [sp_fine, td_fine] = synthetic_case(cfg_fine);
    s_fine = compute_moap_direct(sp_fine, td_fine, cfg_fine);
    coarse_metric = profile_fwhm(s.energy_meV, s.total_raw);
    fine_metric = profile_fwhm(s_fine.energy_meV, s_fine.total_raw);
    coarse_weight = trapz(s.energy_meV, s.total_raw);
    fine_weight = trapz(s_fine.energy_meV, s_fine.total_raw);
    rel_weight_change = abs(fine_weight-coarse_weight) / ...
                        max(abs(fine_weight), realmin);
    assert(abs(fine_metric.peak_x-coarse_metric.peak_x) <= 4.0, ...
        'Representative direct peak is not stable under q-grid refinement.');
    assert(rel_weight_change < 0.35, ...
        'Representative integrated direct weight is not q-grid converged.');

    fprintf(['  PASS: direct q physics (peak %.3g -> %.3g meV, ', ...
             'integrated-weight change %.3g%%)\n'], ...
        coarse_metric.peak_x, fine_metric.peak_x, 100*rel_weight_change);
end

function [sp, td] = synthetic_case(cfg)
    c = cfg.constants;
    sp = struct();
    sp.E_J = [0; 50*c.meV];
    sp.EF_J = 20*c.meV;

    qz = linspace(0, cfg.oap.qz_max_inv_nm, cfg.oap.Nqz) / c.nm;
    I2 = exp(-(2*c.nm*qz).^2);
    td = struct();
    td.initial = 1;
    td.final = 2;
    td.deltaE_J = sp.E_J(2)-sp.E_J(1);
    td.deltaE_meV = td.deltaE_J/c.meV;
    td.qz_inv_m = qz;
    td.I2 = I2;
end

function assert_throws(fun, expected_id)
    did_throw = false;
    try
        fun();
    catch ME
        did_throw = strcmp(ME.identifier, expected_id);
    end
    assert(did_throw, 'Expected error %s was not thrown.', expected_id);
end
