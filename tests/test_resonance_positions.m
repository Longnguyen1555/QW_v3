function test_resonance_positions()
%TEST_RESONANCE_POSITIONS Verify the exact k_parallel delta-function root.

    cfg = default_config();
    c = cfg.constants;
    mstar = cfg.material.mstar_rel*c.m0;
    qperp = 0.42/c.nm;
    deltaE = 31.0*c.meV;
    hw = cfg.material.LO_phonon_meV*c.meV;
    Eph = 78.0*c.meV;
    ell = 1;
    A = c.hbar^2/(2*mstar);

    balance_em = ell*Eph - deltaE - hw;
    k_em = inplane_delta_root(qperp, balance_em, cfg);
    residual_em = deltaE + A*qperp^2 + 2*A*k_em*qperp + hw - ell*Eph;
    assert(abs(residual_em) < 1e-12*c.meV, ...
        'Emission delta root does not satisfy energy conservation.');

    balance_ab = ell*Eph - deltaE + hw;
    k_ab = inplane_delta_root(qperp, balance_ab, cfg);
    residual_ab = deltaE + A*qperp^2 + 2*A*k_ab*qperp - hw - ell*Eph;
    assert(abs(residual_ab) < 1e-12*c.meV, ...
        'Absorption delta root does not satisfy energy conservation.');
    assert(abs(k_em-k_ab) > 0, ...
        'Optical phonons must not collapse to a q-independent resonance.');
    fprintf('  PASS: delta-root energy conservation\n');
end
