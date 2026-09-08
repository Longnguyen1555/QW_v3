function spectrum = compute_moap_direct(sp, td, cfg)
%COMPUTE_MOAP_DIRECT Source-faithful direct q-space OAP calculation.
%
% The in-plane electron dispersion is epsilon(n,k)=epsilon(n)+hbar^2*k^2/(2m*).
% For q_perp along k_parallel, the conservation argument is
%   DeltaE + hbar^2*q_perp^2/(2*m*) + hbar^2*k_parallel*q_perp/m*
%       +/- hbar*omega_q - ell*hbar*Omega.
% k_parallel is eliminated analytically with the Dirac delta.  The remaining
% k_y Fermi integral is evaluated numerically without a lifetime broadening.
% The optional Lorentzian path is solely a display convolution of that physical
% spectrum; metrics use the unbroadened *_binned_raw fields.

    c = cfg.constants;
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);
    dE = Ephot(2) - Ephot(1);
    E0 = cfg.laser.E0_kVcm * 1.0e5;

    qperp = linspace(0, cfg.oap.qperp_max_inv_nm, cfg.oap.Nqperp) / c.nm;
    qz = td(1).qz_inv_m;
    [QP, QZ] = meshgrid(qperp, qz);
    [WP, WZ] = meshgrid(trapezoid_weights(qperp), trapezoid_weights(qz));

    switch lower(cfg.oap.q_model)
        case 'source_qperp'
            % Source approximation q ~= q_perp; qz occurs only in I_nn'(qz).
            q_coupling = QP;
        case 'generalized_3d'
            % Explicit research/comparison model, not the source approximation.
            q_coupling = sqrt(QP.^2 + QZ.^2);
        otherwise
            error('Unknown cfg.oap.q_model: %s', cfg.oap.q_model);
    end

    % int dqz int 2*pi*qperp*dqperp/(2*pi)^3; qz is integrated on [0,inf)
    % and doubled using the form factor's even modulus squared.
    measure = (2.0/(2*pi)^2) .* QP .* WP .* WZ;
    a0 = laser_displacement(Ephot, E0, cfg);
    [qd, ne3d] = debye_wavevector(sp, cfg);
    pref_oap = E0^2 * sqrt(cfg.material.eps_static) / (8*pi) * (2*pi/c.hbar);

    spectrum = initialize_spectrum_struct(Ephot, cfg);
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd*c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.qperp_inv_m = qperp;
    spectrum.meta.q_model = cfg.oap.q_model;
    spectrum.meta.a0_m = a0;
    spectrum.meta.broadening.mode = cfg.oap.broadening.mode;
    spectrum.meta.broadening.gamma_optical_meV = cfg.oap.gamma_optical_meV;
    spectrum.meta.broadening.gamma_piezo_meV = cfg.oap.gamma_piezo_meV;
    spectrum.meta.linewidth_source = 'binned_raw';
    spectrum.meta.kinematics = ['DeltaE + hbar^2*qperp^2/(2*m*) + ', ...
        'hbar^2*kparallel*qperp/m* +/- hbar*omega_q - ell*hbar*Omega'];
    spectrum.meta.absolute_normalization = 'not_source_validated';

    for im = 1:numel(cfg.oap.mechanisms)
        mechanism = cfg.oap.mechanisms{im};
        [C2V, hw, Nph] = phonon_coupling_density(mechanism, q_coupling, qd, cfg);
        [mech_key, gamma] = mechanism_display_parameters(mechanism, cfg);
        mechanism_total = zeros(1, nE);
        mechanism_binned_total = zeros(1, nE);

        for it = 1:numel(td)
            deltaE = td(it).deltaE_J;
            Ei = sp.E_J(td(it).initial);
            I2grid = repmat(td(it).I2(:), 1, numel(qperp));
            base = measure .* C2V .* I2grid;
            if strcmpi(cfg.oap.q_model, 'source_qperp')
                % The source model separates the qz form-factor integral.
                base = sum(base, 1);
                qkin = qperp;
                hwkin = hw(1,:);
                Nphkin = Nph(1,:);
            else
                qkin = QP;
                hwkin = hw;
                Nphkin = Nph;
            end
            transition_total = zeros(1, nE);

            for order = cfg.oap.photon_orders
                dressing_q = QP.^(2*order) ./ (2^(2*order) * factorial(order)^2);
                if strcmpi(cfg.oap.q_model, 'source_qperp')
                    dressing = dressing_q(1,:);
                else
                    dressing = dressing_q;
                end
                raw_em = delta_integrated_channel(Ephot, qkin, base .* dressing .* ...
                    (Nphkin + 1.0), deltaE, hwkin, +1, order, Ei, sp.EF_J, cfg);
                raw_ab = delta_integrated_channel(Ephot, qkin, base .* dressing .* ...
                    Nphkin, deltaE, hwkin, -1, order, Ei, sp.EF_J, cfg);

                energy_pref = pref_oap .* a0.^(2*order);
                physical_em = energy_pref .* raw_em;
                physical_ab = energy_pref .* raw_ab;
                [P_em, P_binned_em] = display_channel(Ephot, physical_em, order, gamma, ...
                    cfg.oap.broadening.mode, dE);
                [P_ab, P_binned_ab] = display_channel(Ephot, physical_ab, order, gamma, ...
                    cfg.oap.broadening.mode, dE);
                P_order = P_em + P_ab;
                P_binned_order = P_binned_em + P_binned_ab;

                order_key = sprintf('order_%d', order);
                component = spectrum.mechanism.(mech_key).(order_key);
                component.emission_raw = component.emission_raw + P_em;
                component.absorption_raw = component.absorption_raw + P_ab;
                component.total_raw = component.total_raw + P_order;
                component.emission_binned_raw = component.emission_binned_raw + P_binned_em;
                component.absorption_binned_raw = component.absorption_binned_raw + P_binned_ab;
                component.total_binned_raw = component.total_binned_raw + P_binned_order;
                spectrum.mechanism.(mech_key).(order_key) = component;

                transition_total = transition_total + P_order;
                mechanism_binned_total = mechanism_binned_total + P_binned_order;
            end

            tr_key = sprintf('transition_%d_%d', td(it).initial, td(it).final);
            spectrum.mechanism.(mech_key).transitions.(tr_key) = transition_total;
            mechanism_total = mechanism_total + transition_total;
        end

        spectrum.mechanism.(mech_key).total_raw = mechanism_total;
        spectrum.mechanism.(mech_key).total_binned_raw = mechanism_binned_total;
        spectrum.total_raw = spectrum.total_raw + mechanism_total;
        spectrum.total_binned_raw = spectrum.total_binned_raw + mechanism_binned_total;
    end

    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end

function raw = delta_integrated_channel(Ephot, qperp, base, deltaE, hw, phonon_sign, order, Ei, EF, cfg)
    raw = zeros(size(Ephot));
    qvec = qperp(:);
    hwvec = hw(:);
    basevec = base(:);
    % Bound the temporary Fermi-quadrature array while vectorizing photon
    % energies.  This is algebraically identical to the scalar loop.
    block_size = max(1, floor(2.0e6 / (cfg.oap.Nkin_perp*numel(qvec))));
    for first = 1:block_size:numel(Ephot)
        last = min(first + block_size - 1, numel(Ephot));
        Eblock = Ephot(first:last);
        % emission: balance = ell*E - DeltaE - hw;
        % absorption: balance = ell*E - DeltaE + hw.
        balance = order*Eblock - deltaE - phonon_sign*hwvec;
        qblock = repmat(qvec, 1, numel(Eblock));
        kparallel = inplane_delta_root(qblock, balance, cfg);
        fermi_density = inplane_delta_fermi_density(qblock, kparallel, Ei, EF, cfg);
        raw(first:last) = sum(basevec .* fermi_density, 1);
    end
end

function a0 = laser_displacement(Ephot, E0, cfg)
    c = cfg.constants;
    switch lower(cfg.laser.a0_mode)
        case 'constant'
            a0 = (cfg.laser.a0_nm*c.nm) .* ones(size(Ephot));
        case 'dynamic'
            if any(Ephot <= 0)
                error('Dynamic a0 requires a strictly positive photon-energy grid.');
            end
            omega = Ephot / c.hbar;
            mstar = cfg.material.mstar_rel * c.m0;
            a0 = c.e * E0 ./ (mstar * omega.^2);
        otherwise
            error('Unknown laser.a0_mode: %s', cfg.laser.a0_mode);
    end
end

function [key, gamma] = mechanism_display_parameters(mechanism, cfg)
    if strcmpi(mechanism, 'optical')
        key = 'optical';
        gamma = cfg.oap.gamma_optical_meV * cfg.constants.meV;
    elseif any(strcmpi(mechanism, {'piezoelectric','piezo'}))
        key = 'piezoelectric';
        gamma = cfg.oap.gamma_piezo_meV * cfg.constants.meV;
    else
        error('Unsupported phonon mechanism: %s', mechanism);
    end
end

function [displayed, physical] = display_channel(Ephot, physical_input, order, gamma, mode, dE)
    % Reuse the established convolution while preserving the exact physical
    % sampled spectrum in the binned output.
    [displayed, physical] = broaden_binned_centers(Ephot, Ephot, ...
        physical_input .* (order*dE), order, gamma, mode);
end
