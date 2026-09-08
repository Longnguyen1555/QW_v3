function spectrum = compute_moap_direct(sp, td, cfg)
%COMPUTE_MOAP_DIRECT Direct evaluation of the document Eq. (2)/(18).
%
% Energy conservation is evaluated through the exact positive-k root of
%
%   B(q_perp) + hbar^2*q_perp*k_perp/m* = 0,
%
% including recoil, the collinear cross term, and the Fermi-Dirac factor at
% that root.  No resonance-center binning or phenomenological linewidth is
% used in this model.
%
% The source approximation q_perp >> q_z is applied deliberately: q_perp
% enters the screened coupling, phonon dispersion, and laser dressing,
% while q_z appears only in |I_n'n(q_z)|^2.  The q_z grid covers [0,Inf),
% so Q_full below explicitly restores the negative-q_z half-space.

    c = cfg.constants;
    orders = cfg.oap.photon_orders(:).';
    if any(~ismember(orders, [1 2])) || any(orders ~= round(orders))
        error('QW:DirectPhotonOrder', ...
            'direct_q_integral supports photon orders ell = 1 and 2 only.');
    end

    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);
    E0 = cfg.laser.E0_kVcm * 1.0e5;
    mstar = cfg.material.mstar_rel * c.m0;

    switch lower(cfg.laser.a0_mode)
        case 'dynamic'
            if any(~isfinite(Ephot) | Ephot <= 0)
                error('QW:DirectPhotonEnergy', ...
                    ['Dynamic a0 requires a finite, strictly positive ', ...
                     'photon-energy grid.']);
            end
            Omega = Ephot / c.hbar;
            a0 = c.e * E0 ./ (mstar * Omega.^2);
        case 'constant'
            warning('QW:NonphysicalConstantA0', ...
                ['constant a0 is a nonphysical/debug compatibility mode ', ...
                 'for direct_q_integral.']);
            Omega = Ephot / c.hbar;
            a0 = cfg.laser.a0_nm * c.nm .* ones(size(Ephot));
        otherwise
            error('QW:DirectA0Mode', 'Unknown laser.a0_mode: %s', ...
                cfg.laser.a0_mode);
    end
    if any(~isfinite(a0) | a0 <= 0)
        error('QW:DirectA0Overflow', ...
            ['a0=eE0/(m*Omega^2) is non-finite. Increase the strictly ', ...
             'positive lower photon-energy bound.']);
    end

    qperp = linspace(0, cfg.oap.qperp_max_inv_nm, ...
                     cfg.oap.Nqperp).' / c.nm;
    wperp = trapezoid_weights(qperp).';

    [qd, ne3d] = debye_wavevector(sp, cfg);
    if ~(isfinite(qd) && qd >= 0)
        error('QW:DirectDebyeWavevector', ...
            'Debye wavevector must be finite and nonnegative.');
    end

    S = cfg.oap.normalization_area_m2;
    if ~(isscalar(S) && isfinite(S) && S > 0)
        error('QW:DirectNormalizationArea', ...
            'cfg.oap.normalization_area_m2 must be finite and positive.');
    end

    if isfield(cfg.oap, 'direct') && ...
            isfield(cfg.oap.direct, 'energy_block_size')
        energy_block_size = cfg.oap.direct.energy_block_size;
    else
        energy_block_size = 2048;
    end
    if ~(isscalar(energy_block_size) && isfinite(energy_block_size) && ...
            energy_block_size >= 1 && ...
            energy_block_size == round(energy_block_size))
        error('QW:DirectEnergyBlockSize', ...
            'direct.energy_block_size must be a positive integer.');
    end

    spectrum = initialize_spectrum_struct(Ephot, cfg);
    spectrum.meta.model = 'direct_q_integral';
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd * c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.qperp_inv_m = qperp;
    spectrum.meta.Omega_rad_s = Omega;
    spectrum.meta.a0_m = a0;
    spectrum.meta.a0_mode = cfg.laser.a0_mode;
    spectrum.meta.normalization_area_m2 = S;
    spectrum.meta.energy_block_size = energy_block_size;
    spectrum.meta.gamma_ignored = true;
    spectrum.meta.q_coupling_approximation = 'q_equals_qperp';
    spectrum.meta.qz_half_space_factor = 2.0;
    spectrum.meta.collinear_kq = true;
    spectrum.meta.piezo_dispersion = 'omega_q_equals_cs_qperp';
    spectrum.meta.raw_normalization = [ ...
        'Source continuum convention with q-volume cancellation and ', ...
        'an explicit in-plane normalization area; absolute-W units ', ...
        'remain source-convention dependent.'];

    % Preserve the source prefactor convention used by this repository:
    % E0^2*sqrt(kappa0)/(8*pi) times the golden-rule 2*pi/hbar.
    % V|C_q|^2 below cancels V from the 3-D q-state density.  The source
    % mixes electromagnetic and SI conventions, so raw values are retained
    % as source-normalized power rather than asserted to be absolute watts.
    source_prefactor = E0^2 * sqrt(cfg.material.eps_static) / (8*pi) * ...
                       (2*pi/c.hbar);
    k_state_prefactor = S / (2*pi);

    qz_full_integrals = zeros(numel(td), 1);

    for im = 1:numel(cfg.oap.mechanisms)
        mechanism = cfg.oap.mechanisms{im};
        terms = direct_phonon_terms(mechanism, qperp, qd, cfg);

        switch lower(mechanism)
            case 'optical'
                mech_key = 'optical';
            case {'piezoelectric', 'piezo'}
                mech_key = 'piezoelectric';
            otherwise
                error('QW:DirectMechanism', ...
                    'Unsupported phonon mechanism: %s', mechanism);
        end

        mechanism_total = zeros(1, nE);

        for it = 1:numel(td)
            initial = td(it).initial;
            deltaE = td(it).deltaE_J;
            Einitial = sp.E_J(initial);

            % Cylindrical q-state density:
            %   1/(2*pi)^3 int d^3q
            % = Q_full/(4*pi^2) int_0^Inf q_perp dq_perp.
            % Q_full = 2*int_0^Inf |I(q_z)|^2 dq_z because the stored q_z
            % grid covers only the positive half-space.
            wqz = trapezoid_weights(td(it).qz_inv_m);
            Qhalf = sum(td(it).I2(:).' .* wqz);
            Qfull = 2.0 * Qhalf;
            qz_full_integrals(it) = Qfull;
            q_measure = Qfull/(4*pi^2) .* qperp .* wperp;

            transition_total = zeros(1, nE);

            for ell = orders
                order_key = sprintf('order_%d', ell);

                % Exactly a0^2*q_perp^2/4 for ell=1 and
                % a0^4*q_perp^4/64 for ell=2.  a0 varies with Eph in the
                % physical dynamic mode.
                P_em = zeros(1, nE);
                P_ab = zeros(1, nE);
                for first = 1:energy_block_size:nE
                    idx = first:min(first+energy_block_size-1, nE);
                    Eph_block = Ephot(idx);
                    a0_block = a0(idx);
                    dressing = (qperp .* a0_block).^(2*ell) ./ ...
                        (2^(2*ell) * factorial(ell)^2);

                    delta_em = direct_k_delta_weight(Einitial, sp.EF_J, ...
                        deltaE, Eph_block, ell, +1, terms.hw_J, qperp, cfg);
                    delta_ab = direct_k_delta_weight(Einitial, sp.EF_J, ...
                        deltaE, Eph_block, ell, -1, terms.hw_J, qperp, cfg);

                    integrand_em = q_measure .* terms.emission_C2V .* ...
                                   delta_em .* dressing;
                    integrand_ab = q_measure .* terms.absorption_C2V .* ...
                                   delta_ab .* dressing;

                    P_em(idx) = source_prefactor * k_state_prefactor .* ...
                                sum(integrand_em, 1);
                    P_ab(idx) = source_prefactor * k_state_prefactor .* ...
                                sum(integrand_ab, 1);
                end

                if any(~isfinite(P_em)) || any(~isfinite(P_ab))
                    error('QW:DirectNonfiniteSpectrum', ...
                        ['Non-finite direct spectrum. Check the positive ', ...
                         'photon-energy grid and q-grid limits.']);
                end

                P_order = P_em + P_ab;
                spectrum.mechanism.(mech_key).(order_key).emission_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).emission_raw + P_em;
                spectrum.mechanism.(mech_key).(order_key).absorption_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).absorption_raw + P_ab;
                spectrum.mechanism.(mech_key).(order_key).total_raw = ...
                    spectrum.mechanism.(mech_key).(order_key).total_raw + P_order;

                transition_total = transition_total + P_order;
            end

            tr_key = sprintf('transition_%d_%d', ...
                td(it).initial, td(it).final);
            spectrum.mechanism.(mech_key).transitions.(tr_key) = ...
                transition_total;
            mechanism_total = mechanism_total + transition_total;
        end

        spectrum.mechanism.(mech_key).total_raw = mechanism_total;
        spectrum.total_raw = spectrum.total_raw + mechanism_total;
    end

    spectrum.meta.qz_form_factor_full_inv_m = qz_full_integrals;
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
