function test_moap_nhan()
%TEST_MOAP_NHAN Formula/regression checks for HƯỚNG DẪN Eqs. (15)-(27).

    cfg = fixture();
    c = cfg.constants;
    mstar = cfg.material.mstar_rel*c.m0;
    kBT = c.kB*cfg.temperature_K;

    % Integer optical and half-integer piezoelectric Bessel orders.
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
        actual = expJ_from_spectrum(D,s,eta,ell,kind,cfg);
        expected = direct_expJ( ...
            D,s,eta,ell,kind,mstar,c.hbar,kBT);
        assert_close(actual,expected,2e-13, ...
            'Nhan Bessel orders/powers do not match the PDF.');
    end

    masked = expJ_from_spectrum( ...
        [5 0 -5]*c.meV,1,0,1,'optical',cfg);
    assert(masked(1) ~= 0 && masked(2) == 0 && masked(3) == 0, ...
        'Nhan formula was extrapolated outside the printed D>0 domain.');

    [sp,td] = synthetic_case(cfg);
    base = compute_moap_nhan(sp,td,cfg);
    assert(all(isfinite(base.total_raw)), ...
        'Nhan spectrum contains non-finite values.');

    % PDF occupation is already in exp[-s(Ei-EF)/kBT].
    cfg_population = cfg;
    cfg_population.oap.population_mode = 'intentionally_unsupported';
    no_extra_population = compute_moap_nhan(sp,td,cfg_population);
    assert(isequal(base.total_raw,no_extra_population.total_raw), ...
        'Nhan spectrum still applies a second population weight.');

    % S is the in-plane area. Lz cannot be substituted for it.
    cfg_area = cfg;
    cfg_area.oap.normalization_area_m2 = 3*cfg.oap.normalization_area_m2;
    area_scaled = compute_moap_nhan(sp,td,cfg_area);
    assert_close(area_scaled.total_raw,3*base.total_raw,5e-13, ...
        'Nhan spectrum does not scale linearly with normalization area.');
    cfg_width = cfg;
    cfg_width.structure.Lz_nm = 4*cfg.structure.Lz_nm;
    width_changed = compute_moap_nhan(sp,td,cfg_width);
    assert(isequal(base.total_raw,width_changed.total_raw), ...
        'Nhan spectrum still uses Lz^2 as its normalization area.');

    % PDF Eq. (26) keeps the same hw0 in D; cs occurs only in Cpiezo.
    cfg_piezo = cfg;
    cfg_piezo.oap.mechanisms = {'piezoelectric'};
    piezo = compute_moap_nhan(sp,td,cfg_piezo);
    cfg_fast = cfg_piezo;
    cfg_fast.material.sound_speed_mps = 2*cfg_piezo.material.sound_speed_mps;
    piezo_fast = compute_moap_nhan(sp,td,cfg_fast);
    assert_close(piezo_fast.total_raw,2*piezo.total_raw,5e-13, ...
        'Nhan D incorrectly depends on hbar*c_s*qd.');

    % Auto mode must fall back to the parent integral outside Ei>EF.
    invalid_sp = sp;
    invalid_sp.EF_J = sp.E_J(1)+c.meV;
    fallback = compute_moap_nhan(invalid_sp,td,cfg);
    assert(all(isfinite(fallback.total_raw)) && ...
        fallback.meta.series_diagnostic.quadrature_points > 0, ...
        'Nhan auto mode did not replace the divergent Fermi series.');
    strict_cfg = cfg;
    strict_cfg.oap.series.evaluation = 'series';
    assert_throws(@() compute_moap_nhan(invalid_sp,td,strict_cfg), ...
        'QW:Nhan:FermiDomain');

    % Signed negative photon energies are retained exactly and remain finite.
    negative_cfg = cfg;
    negative_cfg.oap.photon_energy_meV = [-40 -20 -5];
    negative_cfg.laser.a0_mode = 'dynamic';
    negative = compute_moap_nhan(invalid_sp,td,negative_cfg);
    assert(isequal(negative.energy_meV,negative_cfg.oap.photon_energy_meV) && ...
        all(isfinite(negative.total_raw)), ...
        'Nhan lost or invalidated the negative photon-energy branch.');
    zero_cfg = negative_cfg;
    zero_cfg.oap.photon_energy_meV = [-1 0 1];
    assert_throws(@() compute_moap_nhan(invalid_sp,td,zero_cfg), ...
        'QW:NhanPhotonEnergyZero');

    fprintf('  PASS: HƯỚNG DẪN Nhan Taylor-Bessel formulas\n');
end

function cfg = fixture()
    cfg = default_config();
    cfg.oap.model = 'nhan';
    cfg.oap.plot_normalization = 'none';
    cfg.oap.photon_energy_meV = [2 4 6 8 10];
    cfg.oap.photon_orders = [1 2];
    cfg.oap.mechanisms = {'optical','piezoelectric'};
    cfg.oap.series.s_max = 2;
    cfg.oap.series.eta_max = 2;
    cfg.laser.a0_mode = 'constant';
    cfg.oap.normalization_area_m2 = 2.5;
end

function [sp,td] = synthetic_case(cfg)
    c = cfg.constants;
    sp = struct('E_J',[20;70]*c.meV,'EF_J',0);
    td = struct('initial',1,'final',2, ...
        'deltaE_J',50*c.meV,'deltaE_meV',50, ...
        'Q_half_inv_m',2.0e8);
end

function expJ = expJ_from_spectrum(D,s,eta,ell,kind,cfg)
% Isolate one local (s,eta) term through the public Nhan spectrum API.
    c = cfg.constants;
    cfg.oap.series.evaluation = 'series';
    deltaE = 20*c.meV;
    hw0 = cfg.material.LO_phonon_meV*c.meV;
    cfg.oap.photon_energy_meV = (deltaE+hw0-D)/(ell*c.meV);
    cfg.oap.photon_orders = ell;
    cfg.oap.mechanisms = {kind};

    % qd near the thermal momentum keeps all eta terms well scaled.
    test_qd = 1.0e8;
    ne = test_qd^2*cfg.material.eps_static*c.eps0 * ...
        c.kB*cfg.temperature_K/c.e^2;
    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = ne/1e6;

    sp = struct('E_J',[20;40]*c.meV,'EF_J',0);
    td = struct('initial',1,'final',2, ...
        'deltaE_J',deltaE,'deltaE_meV',deltaE/c.meV, ...
        'Q_half_inv_m',2.0e8);

    rectangle = series_rectangle(sp,td,cfg,s,eta,ell,kind);
    rectangle_s = series_rectangle(sp,td,cfg,s-1,eta,ell,kind);
    rectangle_eta = series_rectangle(sp,td,cfg,s,eta-1,ell,kind);
    rectangle_both = series_rectangle(sp,td,cfg,s-1,eta-1,ell,kind);
    term = rectangle-rectangle_s-rectangle_eta+rectangle_both;

    [qd,~] = debye_wavevector(sp,cfg);
    fermi = exp(-s*(sp.E_J(1)-sp.EF_J)/(c.kB*cfg.temperature_K));
    screening = (-1)^(s+eta)*(eta+1)*qd^(2*eta);
    expJ = term/(screening*fermi);
end

function value = series_rectangle(sp,td,cfg,smax,etamax,ell,kind)
    if smax < 1 || etamax < 0
        value = zeros(size(cfg.oap.photon_energy_meV));
        return;
    end

    cfg.oap.series.s_max = smax;
    cfg.oap.series.eta_max = etamax;
    spectrum = compute_moap_nhan(sp,td,cfg);

    c = cfg.constants;
    E0 = cfg.laser.E0_kVcm*1e5;
    mstar = cfg.material.mstar_rel*c.m0;
    kBT = c.kB*cfg.temperature_K;
    hw0 = cfg.material.LO_phonon_meV*c.meV;
    omega0 = hw0/c.hbar;
    N0 = 1/expm1(hw0/kBT);
    area = cfg.oap.normalization_area_m2;
    volume = area*cfg.oap.active_thickness_nm*c.nm;
    if strcmp(kind,'optical')
        Cmech = E0^2*sqrt(cfg.material.eps_static)/8 * ...
            c.e^2*omega0/(volume*c.eps0) * ...
            (1/cfg.material.eps_high-1/cfg.material.eps_static);
    else
        Cmech = cfg.material.piezo_kappa2*c.e^2*E0^2 * ...
            cfg.material.sound_speed_mps*sqrt(cfg.material.eps_static) / ...
            (8*volume*cfg.material.eps_static*c.eps0);
    end
    common = area*volume*mstar^2/(32*pi^3*c.hbar^4);
    a0 = cfg.laser.a0_nm*c.nm;
    if ell == 1
        field_factor = a0^2;
    else
        field_factor = a0^4/16;
    end
    key = sprintf('order_%d',ell);
    raw = spectrum.mechanism.(kind).(key).emission_raw;
    value = raw/(Cmech*common*td.Q_half_inv_m*(N0+1)*field_factor);
end

function value = direct_expJ(D,s,eta,ell,kind,mstar,hbar,kBT)
    zeta = s*abs(D)/(2*kBT);
    base = 2*mstar*abs(D)/hbar^2;
    if strcmp(kind,'optical')
        power = -eta+ell;
        nu1 = eta-ell+1;
        nu2 = eta-ell;
    else
        power = -eta+ell+0.5;
        nu1 = eta-ell+0.5;
        nu2 = eta-ell-0.5;
    end
    value = exp(s*D/(2*kBT)).*hbar^2/(2*mstar).* ...
        base.^power.*(besselk(nu1,zeta)-besselk(nu2,zeta));
end

function assert_close(actual,expected,tolerance,message)
    scale = max(abs(expected(:)));
    error_size = max(abs(actual(:)-expected(:)));
    assert(error_size <= tolerance*max(scale,realmin),message);
end

function assert_throws(fun,identifier)
    try
        fun();
    catch err
        assert(strcmp(err.identifier,identifier), ...
            'Unexpected error: %s',err.message);
        return;
    end
    error('Expected error was not raised: %s',identifier);
end
