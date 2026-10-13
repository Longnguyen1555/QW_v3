function test_direct_approx_qspace()
%TEST_DIRECT_APPROX_QSPACE Physical-q quadrature and td.I2 reuse regression.

    cfg = default_config();
    cfg.oap.model = 'direct_approx';
    cfg.oap.plot_normalization = 'none';
    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.mechanisms = {'piezoelectric'};
    cfg.oap.photon_orders = 1;
    cfg.oap.photon_energy_meV = 10;
    cfg.oap.qperp_max_inv_nm = 1.2;
    cfg.oap.Nqperp = 7;
    cfg.oap.qz_max_inv_nm = 1.2;
    cfg.oap.Nqz = 13;

    c = cfg.constants;
    sp = struct('E_J',[20 40]*c.meV,'EF_J',40*c.meV);
    qz_positive = linspace(0,1.2,13)/c.nm;
    I2_positive = exp(-(qz_positive*c.nm/0.7).^2);
    td = struct('initial',1,'final',2,'deltaE_J',20*c.meV, ...
        'deltaE_meV',20,'dipole_m',0.5*c.nm, ...
        'Q2_ab_m2',(0.5*c.nm)^2,'qz_inv_m',qz_positive, ...
        'I2',I2_positive);

    spectrum = compute_moap_direct_approx(sp,td,cfg);
    [qd,~] = debye_wavevector(sp,cfg);
    order_factor = (cfg.laser.a0_nm*c.nm)^2/4;
    Omega = cfg.oap.photon_energy_meV*c.meV/c.hbar;
    E0 = cfg.laser.E0_kVcm*1e5;
    me = cfg.material.mstar_rel*c.m0;
    kappa0 = cfg.material.eps_static;
    P0 = E0^4*c.e^4*me^2*cfg.material.piezo_kappa2 * ...
        (c.kB*cfg.temperature_K)*cfg.oap.normalization_area_m2 * ...
        sqrt(kappa0)/(256*pi^3*c.hbar^7*kappa0*c.eps0*Omega);
    actual = spectrum.mechanism.piezoelectric.order_1.total_raw;
    integral = reference_qspace_integral( ...
        sp.E_J(1),sp.EF_J,td.deltaE_J- ...
        cfg.oap.photon_energy_meV*c.meV,1,qd, ...
        qz_positive,I2_positive,cfg);
    expected = P0*td.Q2_ab_m2*order_factor*integral;
    relative_error = abs(actual-expected)/expected;
    assert(relative_error < 1.0e-13, ...
        'direct_approx disagrees with the physical qperp-qz reference.');

    scaled_td = td;
    scaled_td.I2 = 3*td.I2;
    scaled = compute_moap_direct_approx(sp,scaled_td,cfg);
    scale_error = abs(scaled.total_raw-3*spectrum.total_raw) / ...
        (3*spectrum.total_raw);
    assert(scale_error < 1.0e-14, ...
        'direct_approx does not reuse transition_data.I2 linearly.');
    cfg.oap.mechanisms = {'optical','piezoelectric'};
    cfg.oap.photon_orders = [1 2];
    up_to_two = compute_moap_direct_approx(sp,td,cfg);
    cfg.oap.photon_orders = [1 2 3 4];
    with_higher_orders = compute_moap_direct_approx(sp,td,cfg);
    assert(isequaln(with_higher_orders.total_raw,up_to_two.total_raw) && ...
        isequaln(with_higher_orders.mechanism.optical.total_raw, ...
        up_to_two.mechanism.optical.total_raw) && ...
        isequaln(with_higher_orders.mechanism.piezoelectric.total_raw, ...
        up_to_two.mechanism.piezoelectric.total_raw), ...
        'Photon orders above 2 changed the direct_approx spectrum.');
    for order = [3 4]
        key = sprintf('order_%d',order);
        assert(with_higher_orders.mechanism.optical.(key).total_raw == 0 && ...
            with_higher_orders.mechanism.piezoelectric.(key).total_raw == 0, ...
            'direct_approx computed a photon order above 2.');
    end

    fprintf('  PASS: direct_approx physical q-space quadrature\n');
end

function value = reference_qspace_integral( ...
        eps_alpha,EF,D,ell,qd,qz_positive,I2_positive,cfg)

    c = cfg.constants;
    me = cfg.material.mstar_rel*c.m0;
    qperp = linspace(0,cfg.oap.qperp_max_inv_nm, ...
        cfg.oap.Nqperp)/c.nm;
    qz = qz_positive(:);
    I2 = I2_positive(:);

    q2 = qperp.^2+qz.^2;
    E_qperp = c.hbar^2*qperp.^2/(2*me);
    epsilon_star = eps_alpha+E_qperp/4+D/2;
    nonzero = qperp ~= 0;
    epsilon_star(nonzero) = epsilon_star(nonzero)+ ...
        D^2./(4*E_qperp(nonzero));
    epsilon_star(~nonzero) = Inf;
    occupation = stable_fermi_reference( ...
        (epsilon_star-EF)/(c.kB*cfg.temperature_K));
    integrand = q2.^(ell+0.5)./(q2+qd^2).^2 .* ...
        abs(E_qperp+D).*occupation.*I2;
    value = trapz(qperp,trapz(qz,integrand,1),2);
end

function F = stable_fermi_reference(x)
    F = zeros(size(x));
    positive = x >= 0;
    t = exp(-x(positive));
    F(positive) = t./(1+t);
    t = exp(x(~positive));
    F(~positive) = 1./(1+t);
end
