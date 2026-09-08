function test_analytical_series()
%TEST_ANALYTICAL_SERIES Formula and regression tests for Eq. (5)/(7).

    cfg = analytical_fixture();
    c = cfg.constants;
    mstar = cfg.material.mstar_rel*c.m0;
    kBT = c.kB*cfg.temperature_K;

    % A/B: integer-order optical and half-integer piezo Bessel factors.
    D = [5 20]*c.meV;
    cases = { ...
        'optical',       1, 1, 0; ...
        'optical',       3, 2, 2; ...
        'piezoelectric', 1, 1, 0; ...
        'piezoelectric', 2, 2, 0};
    for ic = 1:size(cases,1)
        kind = cases{ic,1};
        s = cases{ic,2};
        ell = cases{ic,3};
        eta = cases{ic,4};
        actual = analytical_exp_bessel_factor( ...
            D, s, eta, ell, kind, mstar, c.hbar, kBT);
        expected = direct_expJ(D, s, eta, ell, kind, mstar, c.hbar, kBT);
        assert_close(actual, expected, 2e-13, ...
            sprintf('%s Bessel orders are incorrect.', kind));
    end

    % The derived closed form is used only on D>0.  D=0 and D<0 are
    % deliberately masked, while non-finite D is reported rather than zeroed.
    [masked, valid] = analytical_exp_bessel_factor( ...
        [5 0 -5]*c.meV, 1, 0, 1, 'optical', mstar, c.hbar, kBT);
    assert(isequal(valid, [true false false]) && ...
           masked(2) == 0 && masked(3) == 0 && isfinite(masked(1)), ...
        'Analytical D>0 domain handling is incorrect.');
    assert_throws(@() analytical_exp_bessel_factor( ...
        NaN, 1, 0, 1, 'optical', mstar, c.hbar, kBT), ...
        'QW:AnalyticalNonfiniteD');

    [sp, td] = synthetic_case(cfg);
    base = compute_moap_analytical_series(sp, td, cfg);
    assert(all(isfinite(base.total_raw)), ...
        'Analytical-series fixture produced non-finite output.');

    % C: population_mode must be irrelevant; Eq. (5)/(7) already contains
    % occupation through exp[-s(E_n-E_F)/(kBT)].
    cfg_population = cfg;
    cfg_population.oap.population_mode = 'intentionally_unsupported';
    no_population_call = compute_moap_analytical_series(sp, td, cfg_population);
    assert(isequal(base.total_raw, no_population_call.total_raw), ...
        'analytical_series still depends on a separate population model.');

    % D: S is the explicit normalization area.  V0 and both couplings use
    % the same volume, so the complete Eq. (5)/(7) result is linear in S.
    cfg_area = cfg;
    cfg_area.oap.normalization_area_m2 = 3*cfg.oap.normalization_area_m2;
    area_scaled = compute_moap_analytical_series(sp, td, cfg_area);
    assert_close(area_scaled.total_raw, 3*base.total_raw, 5e-13, ...
        'Analytical spectrum does not use normalization_area_m2 correctly.');
    cfg_width = cfg;
    cfg_width.structure.Lz_nm = 4*cfg.structure.Lz_nm;
    width_changed = compute_moap_analytical_series(sp, td, cfg_width);
    assert(isequal(base.total_raw, width_changed.total_raw), ...
        'Quantum-well width is still being used as normalization area.');

    % E: analytical Eq. (5)/(7) supports exactly ell=1,2.
    validate_config(cfg);
    cfg_bad_order = cfg;
    cfg_bad_order.oap.photon_orders = [1 2 3];
    assert_throws(@() validate_config(cfg_bad_order), ...
        'QW:AnalyticalPhotonOrder');
    assert_throws(@() compute_moap_analytical_series(sp, td, cfg_bad_order), ...
        'QW:AnalyticalPhotonOrder');

    % F: Eq. (7) uses the same hw0 in D as Eq. (5).  Sound speed therefore
    % enters analytical piezo power only through its linear Cpiezo factor.
    cfg_piezo = cfg;
    cfg_piezo.oap.mechanisms = {'piezoelectric'};
    piezo = compute_moap_analytical_series(sp, td, cfg_piezo);
    cfg_fast_sound = cfg_piezo;
    cfg_fast_sound.material.sound_speed_mps = ...
        2*cfg_piezo.material.sound_speed_mps;
    piezo_fast = compute_moap_analytical_series(sp, td, cfg_fast_sound);
    assert_close(piezo_fast.total_raw, 2*piezo.total_raw, 5e-13, ...
        'Piezo D incorrectly depends on hbar*c_s*qd.');
    assert(strcmp(piezo.meta.piezo_dispersion, 'literal_Eq_7_omega0'), ...
        'Piezo Eq. (7) convention is missing from analytical metadata.');

    % I: constant a0 gives a0^2 (ell=1) and a0^4/16 (ell=2), while
    % dynamic a0=eE0/(m*Omega^2) is applied point-by-point.
    cfg_constant_2 = cfg;
    cfg_constant_2.laser.a0_nm = 2*cfg.laser.a0_nm;
    constant_2 = compute_moap_analytical_series(sp, td, cfg_constant_2);
    assert_close(constant_2.mechanism.optical.order_1.emission_raw, ...
        4*base.mechanism.optical.order_1.emission_raw, 5e-13, ...
        'One-photon constant-a0 scaling is not quadratic.');
    assert_close(constant_2.mechanism.optical.order_2.emission_raw, ...
        16*base.mechanism.optical.order_2.emission_raw, 5e-13, ...
        'Two-photon constant-a0 scaling is not quartic.');

    cfg_dynamic = cfg;
    cfg_dynamic.laser.a0_mode = 'dynamic';
    dynamic = compute_moap_analytical_series(sp, td, cfg_dynamic);
    expected_a0 = c.e*(cfg.laser.E0_kVcm*1e5) ./ ...
        (mstar*(dynamic.energy_J/c.hbar).^2);
    assert_close(dynamic.meta.a0_m, expected_a0, 20*eps, ...
        'Dynamic a0 does not equal eE0/(m*Omega^2).');
    ratio_1 = (expected_a0./base.meta.a0_m).^2;
    ratio_2 = (expected_a0./base.meta.a0_m).^4;
    assert_close(dynamic.mechanism.optical.order_1.emission_raw, ...
        ratio_1.*base.mechanism.optical.order_1.emission_raw, 1e-12, ...
        'Dynamic one-photon a0 was not applied point-by-point.');
    assert_close(dynamic.mechanism.optical.order_2.emission_raw, ...
        ratio_2.*base.mechanism.optical.order_2.emission_raw, 1e-12, ...
        'Dynamic two-photon a0 was not applied point-by-point.');

    % G/H: an eta-truncated fixture exposes a genuine negative raw value;
    % it must survive unchanged instead of being clipped or abs-valued.
    cfg_negative = cfg;
    cfg_negative.oap.series.s_max = 1;
    cfg_negative.oap.series.eta_max = 2;
    negative = compute_moap_analytical_series(sp, td, cfg_negative);
    assert(any(negative.total_raw < 0) && ...
           negative.meta.negative_raw_point_count == nnz(negative.total_raw < 0), ...
        'Negative truncated-series values were not preserved.');

    fprintf('  PASS: analytical Eq. (5)/(7) series\n');
end

function cfg = analytical_fixture()
    cfg = default_config();
    cfg.oap.model = 'analytical_series';
    cfg.oap.plot_normalization = 'none';
    cfg.oap.photon_energy_meV = [2 4 6 8 10];
    cfg.oap.photon_orders = [1 2];
    cfg.oap.mechanisms = {'optical', 'piezoelectric'};
    cfg.oap.series.s_max = 2;
    cfg.oap.series.eta_max = 2;
    cfg.laser.a0_mode = 'constant';
    cfg.oap.normalization_area_m2 = 2.5;
end

function [sp, td] = synthetic_case(cfg)
    c = cfg.constants;
    sp = struct();
    sp.E_J = [20; 70]*c.meV;
    sp.EF_J = 0;

    td = struct();
    td.initial = 1;
    td.final = 2;
    td.deltaE_J = sp.E_J(2)-sp.E_J(1);
    td.deltaE_meV = td.deltaE_J/c.meV;
    td.Q_half_inv_m = 2.0e8;
end

function value = direct_expJ(D, s, eta, ell, kind, mstar, hbar, kBT)
    zeta = s*abs(D)/(2*kBT);
    base = 2*mstar*abs(D)/hbar^2;
    if strcmp(kind, 'optical')
        power = -eta+ell;
        nu1 = eta-ell+1;
        nu2 = eta-ell;
    else
        power = -eta+ell+0.5;
        nu1 = eta-ell+0.5;
        nu2 = eta-ell-0.5;
    end
    J = hbar^2/(2*mstar) .* base.^power .* ...
        (besselk(nu1, zeta)-besselk(nu2, zeta));
    value = exp(s*D/(2*kBT)).*J;
end

function assert_close(actual, expected, relative_tolerance, message)
    scale = max(abs(expected(:)));
    error_size = max(abs(actual(:)-expected(:)));
    assert(error_size <= relative_tolerance*max(scale, realmin), message);
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
