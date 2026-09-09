function spectrum = compute_moap_direct(sp, td, cfg)
%COMPUTE_MOAP_DIRECT Direct q-space evaluation of the second-order model.
%
% This evaluates the continuum transition-probability integral corresponding
% to the target paper's Eqs. (10), (13), (21), (22), and (24)-(26), with:
%   - Dirac delta -> Lorentzian collision broadening;
%   - cylindrical q integration;
%   - Eq. (12) in-plane magnetic form factor;
%   - q_perp in the laser-dressing factor (a0*q_perp)^(2*l).
%
% Process convention:
%   emission   : l*hOmega = DeltaE + hbar*omega_q, Bose factor N_q+1
%   absorption : l*hOmega = DeltaE - hbar*omega_q, Bose factor N_q

    c = cfg.constants;
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);
    dE = Ephot(2)-Ephot(1);

    qperp = linspace(0, cfg.oap.qperp_max_inv_nm, ...
                     cfg.oap.Nqperp) / c.nm;
    wperp = trapezoid_weights(qperp);

    [QP, QZ] = meshgrid(qperp, td(1).qz_inv_m);
    [WP, WZ] = meshgrid(wperp, trapezoid_weights(td(1).qz_inv_m));

    q = sqrt(QP.^2 + QZ.^2);
    J2 = landau_form_factor_squared(QP, cfg);

    % Cylindrical integral:
    % int_{-inf}^{inf}dqz int_0^inf 2*pi*qperp*dqperp /(2*pi)^3
    measure = (2.0/(2*pi)^2) .* QP .* WP .* WZ .* J2;

    populations = electron_populations_for_oap(sp, cfg);
    E0 = cfg.laser.E0_kVcm * 1.0e5;
    pref_oap = E0^2*sqrt(cfg.material.eps_static)/(8*pi) * (2*pi/c.hbar);

    if strcmpi(cfg.laser.a0_mode, 'constant')
        a0 = (cfg.laser.a0_nm*c.nm) .* ones(size(Ephot));
    elseif strcmpi(cfg.laser.a0_mode, 'dynamic')
        omega = Ephot/c.hbar;
        mstar = cfg.material.mstar_rel*c.m0;
        a0 = c.e*E0 ./ (mstar*omega.^2);
    else
        error('Unknown laser.a0_mode.');
    end

    [qd, ne3d] = debye_wavevector(sp, cfg);

    spectrum = initialize_spectrum_struct(Ephot, cfg);
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd*c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.qperp_inv_m = qperp;

    for im = 1:numel(cfg.oap.mechanisms)
        mechanism = cfg.oap.mechanisms{im};
        [C2V, hw, Nph] = phonon_coupling_density(mechanism, q, qd, cfg);

        if strcmpi(mechanism, 'optical')
            gamma = cfg.oap.gamma_optical_meV*c.meV;
            mech_key = 'optical';
        else
            gamma = cfg.oap.gamma_piezo_meV*c.meV;
            mech_key = 'piezoelectric';
        end

        mechanism_total = zeros(1,nE);

        for it = 1:numel(td)
            i = td(it).initial;
            deltaE = td(it).deltaE_J;
            I2grid = repmat(td(it).I2(:), 1, numel(qperp));

            base = measure .* C2V .* I2grid;
            dipole2 = abs(td(it).dipole_m)^2;
            Mrad2 = c.e^2*E0^2*dipole2/4.0;

            transition_total = zeros(1,nE);

            for order = cfg.oap.photon_orders
                dressing_q = QP.^(2*order) ./ ...
                    (2^(2*order) * factorial(order)^2);

                w_em = base .* dressing_q .* (Nph + 1.0);
                w_ab = base .* dressing_q .* Nph;

                centers_em = (deltaE + hw) ./ order;
                centers_ab = (deltaE - hw) ./ order;

                shape_em = broaden_binned_centers(Ephot, centers_em, ...
                    w_em, order, gamma);
                shape_ab = broaden_binned_centers(Ephot, centers_ab, ...
                    w_ab, order, gamma);

                energy_pref = pref_oap .* populations(i) .* Mrad2 .* ...
                    (a0.^(2*order));

                P_em = energy_pref .* shape_em;
                P_ab = energy_pref .* shape_ab;
                P_order = P_em + P_ab;

                order_key = sprintf('order_%d', order);
                spectrum.mechanism.(mech_key).(order_key).emission_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).emission_raw + P_em;
                spectrum.mechanism.(mech_key).(order_key).absorption_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).absorption_raw + P_ab;
                spectrum.mechanism.(mech_key).(order_key).total_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).total_raw + P_order;

                transition_total = transition_total + P_order;
            end

            tr_key = sprintf('transition_%d_%d', td(it).initial, td(it).final);
            spectrum.mechanism.(mech_key).transitions.(tr_key) = transition_total;
            mechanism_total = mechanism_total + transition_total;
        end

        spectrum.mechanism.(mech_key).total_raw = mechanism_total;
        spectrum.total_raw = spectrum.total_raw + mechanism_total;
    end

    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
