function spectrum = compute_moap_direct_approx(sp, td, cfg)
%COMPUTE_MOAP_DIRECT Direct quadrature of Refined PDF Eqs. (19) and (25).

    c = cfg.constants;

    hbar = c.hbar;
    e    = c.e;
    kBT  = c.kB * cfg.temperature_K;
    eps0 = c.eps0;
    me = cfg.material.mstar_rel * c.m0;

    kappa0   = cfg.material.eps_static;
    kappainf = cfg.material.eps_high;

   
    kappa_LO = 1.0 / (1.0/kappainf - 1.0/kappa0);

    hw0 = cfg.material.LO_phonon_meV * c.meV;
    omega0 = hw0 / hbar;

    % Bose population of LO phonon
    N0 = 1.0 / expm1(hw0 / kBT);

    % Strong EM field
    E0 = cfg.laser.E0_kVcm * 1.0e5;

    % Normalization area S from Eqs. (19) and (25).
    S = cfg.oap.normalization_area_m2;

    orders = cfg.oap.photon_orders(:).';
    use_optical = any(strcmpi(cfg.oap.mechanisms,'optical'));
    use_piezoelectric = any(strcmpi( ...
        cfg.oap.mechanisms,'piezoelectric'));

    
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);

    if any(Ephot == 0)
        error('QW:DirectPhotonEnergyZero', ...
            'Refined PDF Eqs. (19) and (25) are singular at hbar*Omega = 0.');
    end

    [qd, ne3d] = debye_wavevector(sp, cfg);
    optical_cache = prepare_optical_cache(sp,td,qd,cfg,use_optical);
    piezo_cache = prepare_piezo_cache( ...
        sp,td,qd,cfg,use_piezoelectric,optical_cache,use_optical);

    a0 = (cfg.laser.a0_nm*c.nm) .* ones(size(Ephot));

    spectrum = initialize_spectrum_struct(Ephot, cfg);

    spectrum.meta.model = 'direct_pdf_eq19_eq25';
    spectrum.meta.source_equations = [19 25];
    spectrum.meta.uses_bessel = false;
    spectrum.meta.piezo_processes_combined = true;
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd * c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.normalization_area_m2 = S;

    % Keep these arrays for later audit.
    %spectrum.meta.a0_m = a0;
    spectrum.meta.Omega_rad_s = zeros(1, nE);
    spectrum.meta.P0_optical = zeros(1, nE);
    spectrum.meta.P0_piezoelectric = zeros(1, nE);
    optical_qperp_bounds = cfg.oap.direct.optical.qperp_bounds_inv_nm;
    optical_qz_bounds = cfg.oap.direct.optical.qz_bounds_inv_nm;
    optical_qperp_step = cfg.oap.direct.optical.qperp_step_inv_nm;
    optical_qz_step = cfg.oap.direct.optical.qz_step_inv_nm;
    spectrum.meta.optical_qperp_bounds_inv_nm = optical_qperp_bounds;
    spectrum.meta.optical_qz_bounds_inv_nm = optical_qz_bounds;
    spectrum.meta.optical_u_bounds = optical_qperp_bounds/(qd*c.nm);
    spectrum.meta.optical_v_bounds = optical_qz_bounds/(qd*c.nm);
    spectrum.meta.optical_qperp_step_inv_nm = optical_qperp_step;
    spectrum.meta.optical_qz_step_inv_nm = optical_qz_step;
    spectrum.meta.optical_u_step = optical_qperp_step/(qd*c.nm);
    spectrum.meta.optical_v_step = optical_qz_step/(qd*c.nm);
    spectrum.meta.optical_base_Nu = numel(optical_cache.u_base);
    spectrum.meta.optical_Nv = numel(optical_cache.v);
    spectrum.meta.optical_quadrature = 'trapz_uv';
    spectrum.meta.optical_u_block_size = optical_cache.u_block_size;
    spectrum.meta.optical_matrix_layout = 'Nv_by_u_block';
    spectrum.meta.optical_form_factor_cached = use_optical;
    piezo_qperp_bounds = cfg.oap.direct.piezo.qperp_bounds_inv_nm;
    piezo_qz_bounds = cfg.oap.direct.piezo.qz_bounds_inv_nm;
    piezo_qperp_step = cfg.oap.direct.piezo.qperp_step_inv_nm;
    piezo_qz_step = cfg.oap.direct.piezo.qz_step_inv_nm;
    spectrum.meta.piezo_qperp_bounds_inv_nm = piezo_qperp_bounds;
    spectrum.meta.piezo_qz_bounds_inv_nm = piezo_qz_bounds;
    spectrum.meta.piezo_u_bounds = piezo_qperp_bounds/(qd*c.nm);
    spectrum.meta.piezo_v_bounds = piezo_qz_bounds/(qd*c.nm);
    spectrum.meta.piezo_qperp_step_inv_nm = piezo_qperp_step;
    spectrum.meta.piezo_qz_step_inv_nm = piezo_qz_step;
    spectrum.meta.piezo_u_step = piezo_qperp_step/(qd*c.nm);
    spectrum.meta.piezo_v_step = piezo_qz_step/(qd*c.nm);
    spectrum.meta.piezo_base_Nu = numel(piezo_cache.u_base);
    spectrum.meta.piezo_Nv = numel(piezo_cache.v);
    spectrum.meta.piezo_quadrature = 'trapz_uv';
    spectrum.meta.piezo_u_block_size = piezo_cache.u_block_size;
    spectrum.meta.piezo_matrix_layout = 'Nv_by_u_block';
    spectrum.meta.piezo_form_factor_cached = use_piezoelectric;
    spectrum.meta.piezo_reuses_optical_form_factor = ...
        piezo_cache.reuses_optical_form_factor;
    spectrum.meta.piezo_qz_domain = 'signed_full_no_symmetry_factor';

    for iE = 1:nE

        Eph = Ephot(iE);

        % hbar Omega = photon energy
        Omega = Eph / hbar;

        spectrum.meta.Omega_rad_s(iE) = Omega;

        %a0 = e * E0 / (me * Omega^2);
    
        P0_LO = E0^4 * e^4 * me^2 * omega0 * S * sqrt(kappa0)/ (256.0 * pi^3 * hbar^6 * eps0 * kappa_LO * Omega^2);

        
        % Eq. (25), collected as Eq. (27): the acoustic prefactor is 1/Omega.
        P0_PE = E0^4 * e^4 * me^2 * cfg.material.piezo_kappa2 * kBT * S * sqrt(kappa0) / (256.0 * pi^3 * hbar^7 * kappa0 * eps0 * Omega);

        spectrum.meta.P0_optical(iE) = P0_LO;
        spectrum.meta.P0_piezoelectric(iE) = P0_PE;

        for it = 1:numel(td)

            alpha = td(it).initial;

            eps_alpha = sp.E_J(alpha);
            EF = sp.EF_J;

            DeltaE = td(it).deltaE_J;

            % Modified Eqs. (19) and (25) keep |Q_ab|^2 outside their 2D
            % integrals and evaluate |I_ab(qz)|^2 inside each integral.
            if isfield(td(it), 'Q2_ab_m2')
                Q2 = td(it).Q2_ab_m2;
            else
                Q2 = abs(td(it).dipole_m)^2;
            end

            optical_transition_factor = Q2;
            piezo_transition_factor = Q2;

            transition_total_LO = 0.0;
            transition_total_PE = 0.0;

            for ell = orders

                order_key = sprintf('order_%d', ell);

                % Exactly:
                %
                % ell = 1 -> a0^2 / 4
                % ell = 2 -> a0^4 / 64
                %
                order_factor = a0(iE)^(2*ell) / ...
                    (2^(2*ell) * factorial(ell)^2);
                %order_factor = a0^(2*ell) / ...
                %   (2^(2*ell) * factorial(ell)^2);

                % ========================================================
                % 1. ELECTRON - LO PHONON
                % ========================================================
               
                if use_optical
                    Dplus = calc_optical_D(DeltaE,hw0,ell,Eph,+1);
                    Dminus = calc_optical_D(DeltaE,hw0,ell,Eph,-1);

                    Tplus = compute_T_optical_2d( ...
                        eps_alpha,EF,Dplus,ell,qd,optical_cache.u_base, ...
                        optical_cache.v,optical_cache.I2{it},cfg);
                    Tminus = compute_T_optical_2d( ...
                        eps_alpha,EF,Dminus,ell,qd,optical_cache.u_base, ...
                        optical_cache.v,optical_cache.I2{it},cfg);

                    Pplus = P0_LO*optical_transition_factor* ...
                        order_factor*(N0+1.0)*Tplus;
                    Pminus = P0_LO*optical_transition_factor* ...
                        order_factor*N0*Tminus;

                    spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                    spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                    spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                    transition_total_LO = transition_total_LO+Pplus+Pminus;
                end

                % ========================================================
                % 2. ELECTRON - PIEZOELECTRIC PHONON
                % ========================================================
                
                if use_piezoelectric
                    Dpiezo = DeltaE - ell * Eph;
                    Tpiezo = compute_T_piezo_2d( ...
                        eps_alpha,EF,Dpiezo,ell,qd,piezo_cache.u_base, ...
                        piezo_cache.v,piezo_cache.I2{it},cfg);
                    Ppiezo = P0_PE*piezo_transition_factor* ...
                        order_factor*Tpiezo;
                    spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + Ppiezo;
                    transition_total_PE = transition_total_PE+Ppiezo;
                end
            end

            tr_key = sprintf('transition_%d_%d', ...
                             td(it).initial, td(it).final);

            if ~isfield(spectrum.mechanism.optical.transitions, tr_key)
                spectrum.mechanism.optical.transitions.(tr_key) = zeros(1, nE);

                spectrum.mechanism.piezoelectric.transitions.(tr_key) = zeros(1, nE);
            end

            spectrum.mechanism.optical.transitions.(tr_key)(iE) = transition_total_LO;

            spectrum.mechanism.piezoelectric.transitions.(tr_key)(iE) = transition_total_PE;
        end
    end

    % ====================================================================
    % Sum mechanisms
    % ====================================================================

    for ell = orders
        order_key = sprintf('order_%d', ell);

        spectrum.mechanism.optical.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.optical.(order_key).total_raw;

        spectrum.mechanism.piezoelectric.total_raw = spectrum.mechanism.piezoelectric.total_raw + spectrum.mechanism.piezoelectric.(order_key).total_raw;
    end

    spectrum.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.piezoelectric.total_raw;

    % plot_normalization must be 'none'.
    
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end

function cache = prepare_piezo_cache(sp,td,qd,cfg,use_piezoelectric, ...
        optical_cache,use_optical)
% Cache the finite Eq. (25) u-v grid and transition form factors.

    c = cfg.constants;
    piezo_cfg = cfg.oap.direct.piezo;
    qperp_grid_inv_nm = build_trapz_grid( ...
        piezo_cfg.qperp_bounds_inv_nm,piezo_cfg.qperp_step_inv_nm);
    qz_grid_inv_nm = build_trapz_grid( ...
        piezo_cfg.qz_bounds_inv_nm,piezo_cfg.qz_step_inv_nm);

    cache.u_base = qperp_grid_inv_nm/(qd*c.nm);
    cache.v = (qz_grid_inv_nm/(qd*c.nm)).';
    cache.u_block_size = piezo_cfg.u_block_size;
    cache.I2 = cell(numel(td),1);
    cache.reuses_optical_form_factor = false;

    if ~use_piezoelectric
        return;
    end
    if use_optical && isequal(cache.v,optical_cache.v)
        cache.I2 = optical_cache.I2;
        cache.reuses_optical_form_factor = true;
        return;
    end

    z_column = sp.z_m(:);
    phase = exp(1i*(z_column*(qd*cache.v.')));
    for it = 1:numel(td)
        wavefunction_product = conj(sp.Psi(:,td(it).initial)) .* ...
            sp.Psi(:,td(it).final);
        cache.I2{it} = calc_I2_v( ...
            z_column,wavefunction_product,phase);
    end
end

function cache = prepare_optical_cache(sp,td,qd,cfg,use_optical)
% Reuse the fixed u-v grids and transition form factors across all energies.

    c = cfg.constants;
    optical_cfg = cfg.oap.direct.optical;
    qperp_grid_inv_nm = build_trapz_grid( ...
        optical_cfg.qperp_bounds_inv_nm,optical_cfg.qperp_step_inv_nm);
    qz_grid_inv_nm = build_trapz_grid( ...
        optical_cfg.qz_bounds_inv_nm,optical_cfg.qz_step_inv_nm);

    cache.u_base = qperp_grid_inv_nm/(qd*c.nm);
    cache.v = (qz_grid_inv_nm/(qd*c.nm)).';
    cache.u_block_size = optical_cfg.u_block_size;
    cache.I2 = cell(numel(td),1);

    if ~use_optical
        return;
    end

    z_column = sp.z_m(:);
    phase = exp(1i*(z_column*(qd*cache.v.')));
    for it = 1:numel(td)
        wavefunction_product = conj(sp.Psi(:,td(it).initial)) .* ...
            sp.Psi(:,td(it).final);
        cache.I2{it} = calc_I2_v( ...
            z_column,wavefunction_product,phase);
    end
end

function Tval = compute_T_optical_2d(eps_alpha_J,EF_J,D_J,ell,qd, ...
        u_base,v,I2_v,cfg)
% Dimensionless u-v form of the modified optical Eq. (19).

    c = cfg.constants;
    hbar = c.hbar;
    kBT = c.kB*cfg.temperature_K;
    me = cfg.material.mstar_rel*c.m0;

    validate_integral_inputs(ell,qd);

    u = u_base;
    E_d = hbar^2*qd^2/(2.0*me);
    if D_J < 0.0
        u_B_zero = sqrt(-D_J/E_d);
        if u_B_zero > u(1) && u_B_zero < u(end)
            base_u_step = cfg.oap.direct.optical.qperp_step_inv_nm / ...
                (qd*c.nm);
            local_half_width = 5.0*base_u_step;
            local_bounds = [max(u(1),u_B_zero-local_half_width) ...
                min(u(end),u_B_zero+local_half_width)];
            local_u = build_trapz_grid(local_bounds,base_u_step/20.0);
            u = sort(unique([u u_B_zero local_u]));
        end
    end

    % qd^(2*ell-1) contains the transformed q kernel and the qd^2
    % Jacobian exactly once.  The full signed v range is integrated, so no
    % additional qz-symmetry factor is present.
    dimensionless_integral = integrate_uv(u,v,I2_v,eps_alpha_J,EF_J, ...
        D_J,ell,qd,E_d,kBT,cfg.oap.direct.optical.u_block_size);
    Tval = qd^(2*ell-1)*dimensionless_integral;
end

function value = integrate_uv(u,v,I2_v,eps_alpha,EF,D,ell,qd,E_d, ...
        kBT,u_block_size)
% Composite trapezoidal quadrature of the full signed u-v domain.

    u = u(:).';
    v = v(:);
    I2_v = I2_v(:);
    integral_over_v = zeros(size(u));
    for first = 1:u_block_size:numel(u)
        indices = first:min(first+u_block_size-1,numel(u));
        integrand = calc_integrand_uv( ...
            u(indices),v,I2_v,eps_alpha,EF,D,ell,qd,E_d,kBT);
        integral_over_v(indices) = trapz(v,integrand,1);
    end
    value = trapz(u,integral_over_v,2);
end

function y = calc_integrand_uv(u,v,I2_v,eps_alpha,EF,D,ell,qd,E_d,kBT)

    qperp = qd*u;
    qz = qd*v;
    rho2 = u.^2+v.^2;
    q2 = qd^2*rho2;
    q = sqrt(q2);
    B = calc_B_u(u,D,E_d);
    epsilon_star = calc_epsilon_star_u(u,D,E_d,eps_alpha);
    occupation = stable_fermi((epsilon_star-EF)/kBT);
    kernel_shape = rho2.^(ell+0.5)./(1.0+rho2).^2;

    % qperp, qz, q2, and q above make the full-vector definition explicit.
    % q is also used in the finite-domain mask; it is never replaced by
    % qperp in the transformed kernel.
    mask = isfinite(q) & isfinite(qperp) & isfinite(qz);
    y = kernel_shape.*abs(B).*occupation.*I2_v;
    y(~mask) = 0.0;
end

function D = calc_optical_D(deltaE,hw0,ell,Eph,phonon_sign)
    D = deltaE+phonon_sign*hw0-ell*Eph;
end

function B = calc_B_u(u,D,E_d)
    B = E_d*u.^2+D;
end

function epsilon_star = calc_epsilon_star_u(u,D,E_d,eps_alpha)
    epsilon_star = eps_alpha+E_d*u.^2/4.0+D/2.0;
    nonzero = u ~= 0;
    epsilon_star(nonzero) = epsilon_star(nonzero)+ ...
        D^2./(4.0*E_d*u(nonzero).^2);
    if D ~= 0
        epsilon_star(~nonzero) = Inf;
    else
        epsilon_star(~nonzero) = eps_alpha;
    end
end

function I2 = calc_I2_v(z_column,wavefunction_product,phase)
% Evaluate cached I_ab(qd*v); the phase exponent is dimensionless.
    product_column = wavefunction_product(:);
    I = trapz(z_column,product_column.*phase,1);
    I2 = abs(I(:)).^2;
end

function grid = build_trapz_grid(bounds,step)
% Include both requested finite endpoints in a uniform trapz grid.
    grid = bounds(1):step:bounds(2);
    if grid(end) < bounds(2)
        grid = [grid bounds(2)];
    end
end

function Tval = compute_T_piezo_2d(eps_alpha_J,EF_J,D_J,ell,qd, ...
        u_base,v,I2_v,cfg)
% Dimensionless u-v form of modified piezoelectric Eq. (25).

    c = cfg.constants;
    hbar = c.hbar;
    kBT = c.kB*cfg.temperature_K;
    me = cfg.material.mstar_rel*c.m0;

    validate_integral_inputs(ell,qd);

    u = u_base;
    E_d = hbar^2*qd^2/(2.0*me);
    if D_J < 0.0
        u_B_zero = sqrt(-D_J/E_d);
        if u_B_zero > u(1) && u_B_zero < u(end)
            base_u_step = cfg.oap.direct.piezo.qperp_step_inv_nm / ...
                (qd*c.nm);
            local_half_width = 5.0*base_u_step;
            local_bounds = [max(u(1),u_B_zero-local_half_width) ...
                min(u(end),u_B_zero+local_half_width)];
            local_u = build_trapz_grid(local_bounds,base_u_step/20.0);
            u = sort(unique([u u_B_zero local_u]));
        end
    end

    % qd^(2*ell-1) contains the transformed full-q kernel and Jacobian.
    % The signed v interval is integrated once, without a symmetry factor.
    dimensionless_integral = integrate_uv(u,v,I2_v,eps_alpha_J,EF_J, ...
        D_J,ell,qd,E_d,kBT,cfg.oap.direct.piezo.u_block_size);
    Tval = qd^(2*ell-1)*dimensionless_integral;
end

function validate_integral_inputs(ell,qd)
    if ~(isscalar(ell) && isfinite(ell) && ell >= 1 && ell == round(ell))
        error('Photon order ell must be a positive integer.');
    end
    if ~(isfinite(qd) && qd > 0)
        error('Debye wavevector qd must be finite and positive.');
    end
end

function F = stable_fermi(x)

    F = zeros(size(x));

    pos = (x >= 0);

    t = exp(-x(pos));
    F(pos) = t ./ (1.0 + t);

    t = exp(x(~pos));
    F(~pos) = 1.0 ./ (1.0 + t);
end
