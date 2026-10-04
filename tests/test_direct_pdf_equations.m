function test_direct_pdf_equations()
%TEST_DIRECT_PDF_EQUATIONS Refined PDF Eqs. (19)/(25) direct regression.

    cfg = default_config();
    cfg.oap.model = 'direct_q_integral';
    cfg.oap.plot_normalization = 'none';
    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.direct.optical.qperp_bounds_inv_nm = [0 50];
    cfg.oap.direct.optical.qperp_step_inv_nm = 0.1;
    cfg.oap.direct.optical.qz_bounds_inv_nm = [-50 50];
    cfg.oap.direct.optical.qz_step_inv_nm = 0.1;
    cfg.oap.direct.piezo.qperp_bounds_inv_nm = [0 50];
    cfg.oap.direct.piezo.qperp_step_inv_nm = 0.1;
    cfg.oap.direct.piezo.qz_bounds_inv_nm = [-50 50];
    cfg.oap.direct.piezo.qz_step_inv_nm = 0.1;
    c = cfg.constants;
    kBT = c.kB*cfg.temperature_K;
    qd = 0.2/c.nm;
    ne = qd^2*cfg.material.eps_static*c.eps0*kBT/c.e^2;
    cfg.oap.fixed_ne_cm3 = ne/1.0e6;

    initial_energy = 50*c.meV;
    transition_energy = 100*c.meV;
    well_width = 10*c.nm;
    z_m = linspace(-well_width/2,well_width/2,121).';
    psi1 = cos(pi*z_m/well_width);
    psi2 = sin(2*pi*z_m/well_width);
    psi1 = psi1/sqrt(trapz(z_m,abs(psi1).^2));
    psi2 = psi2/sqrt(trapz(z_m,abs(psi2).^2));
    sp = struct('E_J',[initial_energy initial_energy+transition_energy], ...
        'EF_J',0,'Nd_sheet_m2',0,'z_m',z_m,'Psi',[psi1 psi2]);
    dipole = trapz(z_m,conj(psi1).*psi2.*z_m);
    td = struct('initial',1,'final',2, ...
        'deltaE_J',transition_energy, ...
        'deltaE_meV',transition_energy/c.meV, ...
        'dipole_m',dipole,'Q_half_inv_m',1.7e8);
    transition_factor = abs(td.dipole_m)^2;

    % Modified Eq. (25): independently evaluate the physical qperp-qz
    % trapz integral for D > 0, D = 0, and D < 0.
    cfg.oap.mechanisms = {'piezoelectric'};
    Dpiezo = [40 0 -40]*c.meV;
    for ell = [1 2 3]
        cfg.oap.photon_orders = ell;
        cfg.oap.photon_energy_meV = ...
            (transition_energy-Dpiezo)/(ell*c.meV);
        spectrum = compute_moap_spectrum(sp,td,cfg);
        key = sprintf('order_%d',ell);
        order_factor = (cfg.laser.a0_nm*c.nm)^(2*ell) / ...
            (2^(2*ell)*factorial(ell)^2);
        actual = spectrum.mechanism.piezoelectric.(key).total_raw ./ ...
            (spectrum.meta.P0_piezoelectric*transition_factor*order_factor);
        expected = arrayfun(@(D) reference_full_q_integral_trapz( ...
            initial_energy,0,D,ell,qd,sp,cfg.oap.direct.piezo,cfg), ...
            Dpiezo);
        relative_error = abs(actual-expected)./expected;
        fprintf('    Eq. (25) ell=%d relative errors:',ell);
        fprintf(' %.3e',relative_error);
        fprintf('\n');
        assert(all(relative_error < 1.0e-12), ...
            'Direct piezoelectric Eq. (25) disagrees with quadrature.');
        assert(relative_spread(spectrum.meta.P0_piezoelectric .* ...
            spectrum.meta.Omega_rad_s) < 1.0e-14, ...
            'Eq. (25)/(27) piezoelectric prefactor must scale as 1/Omega.');
    end

    changed_K = td;
    changed_K.Q_half_inv_m = 11*td.Q_half_inv_m;
    changed_K.K_ab_inv_m = 17*td.Q_half_inv_m;
    changed_cfg = cfg;
    changed_cfg.oap.direct.piezo.u_block_size = 37;
    changed = compute_moap_spectrum(sp,changed_K,changed_cfg);
    assert(isequal(changed.mechanism.piezoelectric.total_raw, ...
        spectrum.mechanism.piezoelectric.total_raw), ...
        ['Modified Eq. (25) depends on an external K_ab or the ', ...
         'configured u block size.']);
    assert(strcmp(spectrum.meta.piezo_quadrature,'trapz_uv') && ...
        spectrum.meta.piezo_base_Nu == 501 && ...
        spectrum.meta.piezo_Nv == 1001 && ...
        spectrum.meta.piezo_form_factor_cached && ...
        strcmp(spectrum.meta.piezo_qz_domain, ...
            'signed_full_no_symmetry_factor') && ...
        strcmp(spectrum.meta.piezo_matrix_layout,'Nv_by_u_block'), ...
        'Modified Eq. (25) quadrature metadata is incorrect.');

    % Modified Eq. (19): verify the full q=sqrt(qperp^2+qz^2) integral,
    % signed qz range, in-integral form factor, and absence of external K_ab.
    cfg.oap.mechanisms = {'optical'};
    ell = 1;
    cfg.oap.photon_orders = ell;
    hw0 = cfg.material.LO_phonon_meV*c.meV;
    cfg.oap.photon_energy_meV = (transition_energy+hw0)/c.meV;
    validate_uv_transformation(qd,ell,initial_energy,0,-40*c.meV,cfg);
    spectrum = compute_moap_spectrum(sp,td,cfg);
    key = sprintf('order_%d',ell);
    order_factor = (cfg.laser.a0_nm*c.nm)^(2*ell) / ...
        (2^(2*ell)*factorial(ell)^2);
    N0 = 1/expm1(hw0/kBT);
    Eph = cfg.oap.photon_energy_meV*c.meV;
    Dplus = transition_energy+hw0-Eph;
    Dminus = transition_energy-hw0-Eph;
    Q2 = abs(td.dipole_m)^2;
    actual_plus = spectrum.mechanism.optical.(key).emission_raw / ...
        (spectrum.meta.P0_optical*Q2*order_factor*(N0+1));
    actual_minus = spectrum.mechanism.optical.(key).absorption_raw / ...
        (spectrum.meta.P0_optical*Q2*order_factor*N0);
    expected_plus = reference_optical_integral_2d( ...
        initial_energy,0,Dplus,ell,qd,sp,cfg);
    expected_minus = reference_optical_integral_2d( ...
        initial_energy,0,Dminus,ell,qd,sp,cfg);
    optical_plus_error = abs(actual_plus-expected_plus)/expected_plus;
    optical_minus_error = abs(actual_minus-expected_minus)/expected_minus;
    fprintf('    Eq. (19) optical relative errors: %.3e %.3e\n', ...
        optical_plus_error,optical_minus_error);
    assert(optical_plus_error < 5.0e-4, ...
        'Modified optical emission term disagrees with Eq. (19).');
    assert(optical_minus_error < 5.0e-4, ...
        'Modified optical absorption term disagrees with Eq. (19).');

    changed_K = td;
    changed_K.Q_half_inv_m = 11*td.Q_half_inv_m;
    changed_cfg = cfg;
    changed_cfg.oap.direct.optical.u_block_size = 37;
    changed = compute_moap_spectrum(sp,changed_K,changed_cfg);
    assert(isequal(changed.mechanism.optical.total_raw, ...
        spectrum.mechanism.optical.total_raw), ...
        ['Modified Eq. (19) depends on an external K_ab or the ', ...
         'configured u block size.']);
    isolated_cfg = cfg;
    isolated_cfg.oap.direct.piezo.qperp_bounds_inv_nm = [0 3];
    isolated_cfg.oap.direct.piezo.qperp_step_inv_nm = 0.3;
    isolated_cfg.oap.direct.piezo.qz_bounds_inv_nm = [-4 4];
    isolated_cfg.oap.direct.piezo.qz_step_inv_nm = 0.4;
    isolated_cfg.oap.direct.piezo.u_block_size = 11;
    isolated = compute_moap_spectrum(sp,td,isolated_cfg);
    assert(isequal(isolated.mechanism.optical.total_raw, ...
        spectrum.mechanism.optical.total_raw), ...
        'The piezoelectric grid contaminated the optical calculation.');
    assert(isequal(spectrum.meta.optical_qperp_bounds_inv_nm,[0 50]) && ...
        isequal(spectrum.meta.optical_qz_bounds_inv_nm,[-50 50]), ...
        'Modified Eq. (19) integration bounds are incorrect.');
    assert(strcmp(spectrum.meta.optical_quadrature,'trapz_uv') && ...
        spectrum.meta.optical_base_Nu == 501 && ...
        spectrum.meta.optical_Nv == 1001 && ...
        strcmp(spectrum.meta.piezo_quadrature,'trapz_uv') && ...
        spectrum.meta.optical_form_factor_cached && ...
        spectrum.meta.optical_u_block_size == ...
            cfg.oap.direct.optical.u_block_size && ...
        strcmp(spectrum.meta.optical_matrix_layout,'Nv_by_u_block'), ...
        'Modified Eq. (19) quadrature metadata is incorrect.');
    expected_cutoff = 50/(qd*c.nm);
    cutoff_tolerance = 10*eps(expected_cutoff);
    assert(max(abs(spectrum.meta.optical_u_bounds- ...
        [0 expected_cutoff])) < cutoff_tolerance && ...
        max(abs(spectrum.meta.optical_v_bounds- ...
        [-expected_cutoff expected_cutoff])) < cutoff_tolerance, ...
        'Physical q cutoffs were not transformed to u-v limits.');

    combined_cfg = cfg;
    combined_cfg.oap.mechanisms = {'optical','piezoelectric'};
    combined = compute_moap_spectrum(sp,td,combined_cfg);
    assert(combined.meta.piezo_reuses_optical_form_factor, ...
        'Matching optical/piezo qz grids did not reuse the form-factor cache.');

    assert(strcmp(spectrum.meta.model,'direct_pdf_eq19_eq25') && ...
        isequal(spectrum.meta.source_equations,[19 25]) && ...
        ~spectrum.meta.uses_bessel, ...
        'Direct-kernel source-equation metadata is incorrect.');

    zero_cfg = cfg;
    zero_cfg.oap.photon_energy_meV = 0;
    expect_error(@() compute_moap_direct(sp,td,zero_cfg), ...
        'QW:DirectPhotonEnergyZero');

    fprintf('  PASS: direct Refined PDF Eqs. (19)/(25)\n');
end

function value = reference_full_q_integral_trapz( ...
        eps_alpha,EF,D,ell,qd,sp,grid_cfg,cfg)
% Independent physical-q implementation of the full Eq. (19)/(25) kernel.

    c = cfg.constants;
    kBT = c.kB*cfg.temperature_K;
    me = cfg.material.mstar_rel*c.m0;
    qp_nm = build_reference_grid(grid_cfg.qperp_bounds_inv_nm, ...
        grid_cfg.qperp_step_inv_nm);
    if D < 0
        qp_zero_nm = sqrt(-2*me*D/c.hbar^2)*c.nm;
        if qp_zero_nm > qp_nm(1) && qp_zero_nm < qp_nm(end)
            half_width = 5*grid_cfg.qperp_step_inv_nm;
            local_bounds = [max(qp_nm(1),qp_zero_nm-half_width) ...
                min(qp_nm(end),qp_zero_nm+half_width)];
            local_qp = build_reference_grid( ...
                local_bounds,grid_cfg.qperp_step_inv_nm/20);
            qp_nm = sort(unique([qp_nm qp_zero_nm local_qp]));
        end
    end
    qz_nm = build_reference_grid( ...
        grid_cfg.qz_bounds_inv_nm,grid_cfg.qz_step_inv_nm).';
    qperp = qp_nm/c.nm;
    qz = qz_nm/c.nm;
    q2 = qperp.^2+qz.^2;

    product = conj(sp.Psi(:,1)).*sp.Psi(:,2);
    I2 = form_factor_at_qz(qz,sp.z_m,product);
    integrand = zeros(size(q2));
    nonzero_qperp = qperp > 0;
    qp = qperp(nonzero_qperp);
    epsilon_star = eps_alpha+c.hbar^2*qp.^2/(8*me)+D/2+ ...
        me*D^2./(2*c.hbar^2*qp.^2);
    occupation = stable_fermi((epsilon_star-EF)/kBT);
    integrand(:,nonzero_qperp) = ...
        q2(:,nonzero_qperp).^((2*ell+1)/2) ./ ...
        (q2(:,nonzero_qperp)+qd^2).^2 .* ...
        abs(c.hbar^2*qp.^2/(2*me)+D) .* occupation .* I2;
    value = trapz(qperp,trapz(qz,integrand,1),2);
end

function validate_uv_transformation(qd,ell,eps_alpha,EF,D,cfg)
% Acceptance checks A--F for the exact qperp=qd*u, qz=qd*v transform.
    c = cfg.constants;
    me = cfg.material.mstar_rel*c.m0;
    E_d = c.hbar^2*qd^2/(2*me);
    u = 0.73;
    v = -1.21;
    qperp = qd*u;
    qz = qd*v;
    rho2 = u^2+v^2;
    q2 = qperp^2+qz^2;
    assert(abs(qperp/qd-u) < 10*eps && abs(qz/qd-v) < 10*eps, ...
        'u-v wave-vector definitions are inconsistent.');
    assert(abs(q2-qd^2*rho2)/q2 < 10*eps, ...
        'Full q^2 transformation is inconsistent.');

    physical_kernel_with_jacobian = ...
        q2^((2*ell+1)/2)/(q2+qd^2)^2*qd^2;
    transformed_kernel = qd^(2*ell-1)* ...
        rho2^(ell+0.5)/(1+rho2)^2;
    assert(abs(physical_kernel_with_jacobian-transformed_kernel) / ...
        abs(transformed_kernel) < 1.0e-14, ...
        'Transformed q kernel or Jacobian is incorrect.');

    B_physical = c.hbar^2*qperp^2/(2*me)+D;
    B_transformed = E_d*u^2+D;
    assert(abs(B_physical-B_transformed) < 10*eps(abs(B_transformed)), ...
        'Transformed B(u) is incorrect.');

    epsilon_physical = eps_alpha+c.hbar^2*qperp^2/(8*me)+D/2+ ...
        me*D^2/(2*c.hbar^2*qperp^2);
    epsilon_transformed = eps_alpha+E_d*u^2/4+D/2+ ...
        D^2/(4*E_d*u^2);
    assert(abs(epsilon_physical-epsilon_transformed) < ...
        20*eps(abs(epsilon_transformed)), ...
        'Transformed epsilon_star(u) is incorrect.');
    assert(isfinite((epsilon_transformed-EF)/(c.kB*cfg.temperature_K)), ...
        'Transformed Fermi exponent is invalid.');
end

function value = reference_optical_integral_2d( ...
        eps_alpha,EF,D,ell,qd,sp,cfg)
    c = cfg.constants;
    me = cfg.material.mstar_rel*c.m0;
    kBT = c.kB*cfg.temperature_K;
    q_unit = 1/c.nm;
    product = conj(sp.Psi(:,1)).*sp.Psi(:,2);
    f = @(qp_nm,qz_nm) reference_optical_integrand( ...
        qp_nm,qz_nm,eps_alpha,EF,D,ell,qd,me,kBT,c.hbar, ...
        q_unit,sp.z_m,product);
    qp_bounds = [0 50];
    qz_bounds = [-50 50];
    value = q_unit^2*integral2(f,qp_bounds(1),qp_bounds(2), ...
        qz_bounds(1),qz_bounds(2),'Method','iterated', ...
        'RelTol',2.0e-8,'AbsTol',0.0);
end

function y = reference_optical_integrand(qp_nm,qz_nm,eps_alpha,EF,D, ...
        ell,qd,me,kBT,hbar,q_unit,z_m,product)
    qp = qp_nm*q_unit;
    qz = qz_nm*q_unit;
    q2 = qp.^2+qz.^2;
    form_factor = form_factor_at_qz(qz,z_m,product);
    y = zeros(size(qp));
    mask = qp > 0;
    qpm = qp(mask);
    eps_star = eps_alpha+hbar^2*qpm.^2/(8*me)+D/2+ ...
        me*D^2./(2*hbar^2*qpm.^2);
    occupation = stable_fermi((eps_star-EF)/kBT);
    y(mask) = q2(mask).^((2*ell+1)/2)./(q2(mask)+qd^2).^2 .* ...
        abs(hbar^2*qpm.^2/(2*me)+D).*occupation.*form_factor(mask);
end

function I2 = form_factor_at_qz(qz,z_m,product)
    output_size = size(qz);
    phase = exp(1i*(z_m(:)*qz(:).'));
    I = trapz(z_m(:),product(:).*phase,1);
    I2 = reshape(abs(I).^2,output_size);
end

function grid = build_reference_grid(bounds,step)
    grid = bounds(1):step:bounds(2);
    if grid(end) < bounds(2)
        grid = [grid bounds(2)];
    end
end

function value = stable_fermi(x)
    value = zeros(size(x));
    pos = x >= 0;
    t = exp(-x(pos));
    value(pos) = t./(1+t);
    t = exp(x(~pos));
    value(~pos) = 1./(1+t);
end

function spread = relative_spread(values)
    spread = (max(values)-min(values))/max(abs(values));
end

function expect_error(f,identifier)
    try
        f();
    catch err
        assert(strcmp(err.identifier,identifier), ...
            'Unexpected error: %s',err.message);
        return;
    end
    error('Expected error was not raised: %s',identifier);
end
