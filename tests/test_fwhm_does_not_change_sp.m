function test_fwhm_does_not_change_sp()
%TEST_FWHM_DOES_NOT_CHANGE_SP Broadening mode cannot affect SP outputs.

    for B = [0 5 10 15]
        cfg_none = apply_numerical_profile(default_config(), 'quick');
        cfg_none.fields.B_T = B;
        cfg_none.oap.broadening.mode = 'none';
        cfg_legacy = cfg_none;
        cfg_legacy.oap.broadening.mode = 'legacy_lorentzian';

        sp_none = solve_schrodinger_poisson(cfg_none);
        sp_legacy = solve_schrodinger_poisson(cfg_legacy);

        assert(max(abs(sp_none.E_J - sp_legacy.E_J)) < 1e-18, ...
            'Broadening mode changed subband energies.');
        assert(abs(sp_none.EF_J - sp_legacy.EF_J) < 1e-18, ...
            'Broadening mode changed the Fermi level.');
        assert(max(abs(sp_none.VH_J - sp_legacy.VH_J)) < 1e-18, ...
            'Broadening mode changed the Hartree potential.');
        overlap = abs(trapz(sp_none.z_m, conj(sp_none.Psi).*sp_legacy.Psi, 1));
        assert(min(overlap) > 1 - 1e-10, ...
            'Broadening mode changed SP wavefunctions.');
    end
    fprintf('  PASS: FWHM settings leave SP unchanged\n');
end
