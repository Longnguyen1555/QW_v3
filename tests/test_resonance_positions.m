function test_resonance_positions()
% For a constant LO phonon, expected peaks are
% (DeltaE +/- hbar*omega_LO)/ell within Lorentzian/grid tolerance.

    cfg = apply_numerical_profile(default_config(),'quick');
    cfg.oap.mechanisms = {'optical'};
    cfg.oap.photon_orders = [1 2];
    sp = solve_schrodinger_poisson(cfg);
    td = compute_transition_data(sp,cfg);
    s = compute_moap_spectrum(sp,td,cfg);

    hw = cfg.material.LO_phonon_meV;
    dE = td(1).deltaE_meV;
    tolerance = 2.0;

    for ell = cfg.oap.photon_orders
        ok = sprintf('order_%d',ell);
        em = profile_fwhm(s.energy_meV, ...
             s.mechanism.optical.(ok).emission_raw);
        expected_em = (dE+hw)/ell;
        assert(abs(em.peak_x-expected_em) < tolerance, ...
            'Optical emission resonance-position test failed.');

        expected_ab = (dE-hw)/ell;
        if expected_ab > min(s.energy_meV)
            ab = profile_fwhm(s.energy_meV, ...
                 s.mechanism.optical.(ok).absorption_raw);
            assert(abs(ab.peak_x-expected_ab) < tolerance, ...
                'Optical absorption resonance-position test failed.');
        end
    end
    fprintf('  PASS: resonance positions\n');
end
