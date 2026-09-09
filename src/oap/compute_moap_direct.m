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

    dbg = isfield(cfg.run, 'debug_moap_direct') && ...
      cfg.run.debug_moap_direct;

    if dbg
        fprintf('\n========================================\n');
        fprintf('DEBUG compute_moap_direct\n');
        fprintf('========================================\n');
        fprintf('model          = %s\n', cfg.oap.model);
        fprintf('B              = %.6g T\n', cfg.fields.B_T);
        fprintf('T              = %.6g K\n', cfg.temperature_K);
        fprintf('a0_mode        = %s\n', cfg.laser.a0_mode);
        fprintf('inplane_factor = %s\n', cfg.oap.inplane_factor);
        fprintf('Nqz            = %d\n', cfg.oap.Nqz);
        fprintf('Nqperp         = %d\n', cfg.oap.Nqperp);
    end


    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;

    assert(all(isfinite(Ephot)), ...
        'Ephot contains NaN/Inf.');

    assert(all(Ephot > 0), ...
        'Ephot must be strictly positive.');

    assert(all(diff(Ephot) > 0), ...
        'Ephot must be strictly increasing.');

    if dbg
        fprintf('Ephot = %.6f -> %.6f meV, dE = %.6g meV\n', ...
            Ephot(1)/c.meV, Ephot(end)/c.meV, ...
            (Ephot(2)-Ephot(1))/c.meV);
    end


    nE = numel(Ephot);
 

    qperp = linspace(0, cfg.oap.qperp_max_inv_nm, ...
                     cfg.oap.Nqperp) / c.nm;
    wperp = trapezoid_weights(qperp);

    [QP, QZ] = meshgrid(qperp, td(1).qz_inv_m);
    [WP, WZ] = meshgrid(wperp, trapezoid_weights(td(1).qz_inv_m));

    q = sqrt(QP.^2 + QZ.^2);
    J2 = landau_form_factor_squared(QP, cfg);


    assert(all(isfinite(q(:))), 'q contains NaN/Inf.');
    assert(all(q(:) >= 0), 'q contains negative values.');

    assert(all(isfinite(J2(:))), 'J2 contains NaN/Inf.');
    assert(all(J2(:) >= 0), 'J2 contains negative values.');

    if dbg
        fprintf('\n[q grid]\n');
        fprintf('q_perp max = %.6g nm^-1\n', max(QP(:))*c.nm);
        fprintf('q_z max    = %.6g nm^-1\n', max(QZ(:))*c.nm);
        fprintf('q max      = %.6g nm^-1\n', max(q(:))*c.nm);

        fprintf('\n[in-plane factor]\n');
        fprintf('min(J2) = %.6e\n', min(J2(:)));
        fprintf('max(J2) = %.6e\n', max(J2(:)));
    end

    measure = (2.0/(2*pi)^2) .* QP .* WP .* WZ .* J2;
    %measure = (2.0/(2*pi)^2) .* QP .* WP .* WZ;

    assert(all(isfinite(measure(:))), ...
        'measure contains NaN/Inf.');

    assert(all(measure(:) >= 0), ...
        'measure contains negative values.');

    if dbg
        fprintf('sum(measure) = %.12e\n', sum(measure(:)));
    end

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

        assert(all(isfinite(C2V(:))), ...
            '%s C2V contains NaN/Inf.', mechanism);

        assert(all(C2V(:) >= 0), ...
            '%s C2V contains negative values.', mechanism);

        assert(all(isfinite(hw(:))), ...
            '%s phonon energy contains NaN/Inf.', mechanism);

        assert(all(hw(:) >= 0), ...
            '%s phonon energy contains negative values.', mechanism);

        assert(all(isfinite(Nph(:))), ...
            '%s Bose factor contains NaN/Inf.', mechanism);

        assert(all(Nph(:) >= 0), ...
            '%s Bose factor contains negative values.', mechanism);

        if dbg
            fprintf('\n[%s]\n', mechanism);
            fprintf('C2V min/max = %.6e / %.6e\n', ...
                min(C2V(:)), max(C2V(:)));

            fprintf('hw min/max  = %.6f / %.6f meV\n', ...
                min(hw(:))/c.meV, max(hw(:))/c.meV);

            fprintf('Nph min/max = %.6e / %.6e\n', ...
                min(Nph(:)), max(Nph(:)));
        end

        if strcmpi(mechanism, 'optical')
            gamma = cfg.oap.gamma_optical_meV*c.meV;
            mech_key = 'optical';
        else
            gamma = cfg.oap.gamma_piezo_meV*c.meV;
            mech_key = 'piezoelectric';
        end

        if dbg && strcmpi(mechanism, 'optical')
            spread_hw = max(hw(:))-min(hw(:));

            fprintf('LO hw spread = %.12e meV\n', ...
                spread_hw/c.meV);

            assert(spread_hw/c.meV < 1e-10, ...
                'Optical phonon is unexpectedly dispersive.');
        end

        mechanism_total = zeros(1,nE);

        for it = 1:numel(td)
            i = td(it).initial;
            deltaE = td(it).deltaE_J;
            I2grid = repmat(td(it).I2(:), 1, numel(qperp));

            base = measure .* C2V .* I2grid;

            assert(all(isfinite(I2grid(:))), ...
                'I2grid contains NaN/Inf.');

            assert(all(I2grid(:) >= 0), ...
                'I2grid contains negative values.');

            assert(all(isfinite(base(:))), ...
                'base contains NaN/Inf.');

            assert(all(base(:) >= 0), ...
                'base contains negative values.');

            if dbg
                fprintf('\n[transition %d -> %d]\n', ...
                    td(it).initial, td(it).final);

                fprintf('DeltaE = %.9f meV\n', ...
                    deltaE/c.meV);

                fprintf('dipole = %.9e nm\n', ...
                    td(it).dipole_m/c.nm);

                fprintf('Q_half = %.9e m^-1\n', ...
                    td(it).Q_half_inv_m);

                fprintf('sum(base) = %.12e\n', ...
                    sum(base(:)));
            end

            if dbg
                fprintf('|I(qz=0)|^2 = %.12e\n', td(it).I2(1));
                fprintf('max |I(qz)|^2 = %.12e\n', max(td(it).I2));
            end
            
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

                if dbg
                    cem = centers_em(:)/c.meV;
                    cab = centers_ab(:)/c.meV;

                    fprintf('\n  order = %d\n', order);
                    fprintf('  emission center   = %.9f -> %.9f meV\n', ...
                        min(cem), max(cem));
                    fprintf('  absorption center = %.9f -> %.9f meV\n', ...
                        min(cab), max(cab));
                end

                shape_em = broaden_binned_centers(Ephot, centers_em, ...
                    w_em, order, gamma);
                shape_ab = broaden_binned_centers(Ephot, centers_ab, ...
                    w_ab, order, gamma);

                if dbg
                    EmeV = Ephot/c.meV;

                    ms_em = profile_fwhm(EmeV, shape_em);
                    ms_ab = profile_fwhm(EmeV, shape_ab);

                    expected_fwhm = ...
                        2*(gamma/c.meV)/order;

                    fprintf('  gamma             = %.9f meV\n', ...
                        gamma/c.meV);

                    fprintf('  expected 2G/l     = %.9f meV\n', ...
                        expected_fwhm);

                    fprintf('  shape emission FWHM = %.9f meV\n', ...
                        ms_em.fwhm);

                    fprintf('  shape absorption FWHM = %.9f meV\n', ...
                        ms_ab.fwhm);
                end

                energy_pref = pref_oap .* populations(i) .* Mrad2 .* ...
                    (a0.^(2*order));

                P_em = energy_pref .* shape_em;
                P_ab = energy_pref .* shape_ab;

                if dbg
                    mp_em = profile_fwhm(Ephot/c.meV, P_em);
                    mp_ab = profile_fwhm(Ephot/c.meV, P_ab);

                    fprintf('  final P emission peak = %.9f meV\n', ...
                        mp_em.peak_x);

                    fprintf('  final P emission FWHM = %.9f meV\n', ...
                        mp_em.fwhm);

                    fprintf('  final P absorption peak = %.9f meV\n', ...
                        mp_ab.peak_x);

                    fprintf('  final P absorption FWHM = %.9f meV\n', ...
                        mp_ab.fwhm);
                end


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
