function spectrum = compute_moap_analytical_series(sp, td, cfg)
%COMPUTE_MOAP_ANALYTICAL_SERIES Literal truncated Eq. (5)/(7) series.
%
% This comparison model encodes the final analytical formulas in the
% calculation note.  It is intentionally separate from direct_q_integral:
% no q integration, collision broadening, or q-dependent acoustic phonon
% dispersion is introduced here.

    c = cfg.constants;
    orders = cfg.oap.photon_orders(:).';
    if any(orders ~= round(orders)) || any(~ismember(orders, [1 2]))
        error('QW:AnalyticalPhotonOrder', ...
            ['analytical_series implements only photon orders ell = 1,2 ', ...
             'according to Eq. (5)/(7).']);
    end
    if ~(isscalar(cfg.oap.series.s_max) && ...
            isfinite(cfg.oap.series.s_max) && cfg.oap.series.s_max >= 1 && ...
            cfg.oap.series.s_max == round(cfg.oap.series.s_max)) || ...
            ~(isscalar(cfg.oap.series.eta_max) && ...
              isfinite(cfg.oap.series.eta_max) && cfg.oap.series.eta_max >= 0 && ...
              cfg.oap.series.eta_max == round(cfg.oap.series.eta_max))
        error('QW:AnalyticalSeriesTruncation', ...
            'series.s_max >= 1 and series.eta_max >= 0 must be integers.');
    end

    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    spectrum = initialize_spectrum_struct(Ephot, cfg);

    [qd, ne3d] = debye_wavevector(sp, cfg);
    if ~(isscalar(qd) && isfinite(qd) && qd >= 0)
        error('QW:AnalyticalDebyeWavevector', ...
            'Debye wavevector must be finite and nonnegative.');
    end

    E0 = cfg.laser.E0_kVcm * 1.0e5;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    hw0 = cfg.material.LO_phonon_meV * c.meV;
    omega0 = hw0 / c.hbar;
    N0 = 1 / expm1(hw0 / kBT);

    % S [m^2] is the electron k-space normalization area, not Lz^2.
    S = cfg.oap.normalization_area_m2;
    if ~(isscalar(S) && isfinite(S) && S > 0)
        error('QW:AnalyticalNormalizationArea', ...
            'cfg.oap.normalization_area_m2 must be finite and positive.');
    end
    active_thickness = cfg.oap.active_thickness_nm * c.nm;
    if ~(isscalar(active_thickness) && isfinite(active_thickness) && ...
            active_thickness > 0)
        error('QW:AnalyticalNormalizationVolume', ...
            'cfg.oap.active_thickness_nm must be finite and positive.');
    end
    V0 = S * active_thickness; % [m^3] phonon normalization volume

    [a0, Omega] = analytical_laser_displacement(Ephot, cfg);

    % The source gives no separate electron/phonon normalization volumes.
    % Use the explicit project normalization volume for V = Vpi = V0.
    Copt = E0^2 * sqrt(cfg.material.eps_static) / 8 * ...
        c.e^2 * omega0 / (V0 * c.eps0) * ...
        (1/cfg.material.eps_high - 1/cfg.material.eps_static);
    Cpiezo = cfg.material.piezo_kappa2 * c.e^2 * E0^2 * ...
        cfg.material.sound_speed_mps * sqrt(cfg.material.eps_static) / ...
        (8 * V0 * cfg.material.eps_static * c.eps0);

    % Common prefactor in both Eq. (5) and Eq. (7).  The photon-order
    % field factors a0^2 and a0^4/16 are applied point-by-point below.
    common = S * V0 * mstar^2 / (32*pi^3*c.hbar^4);

    spectrum.meta.model = 'analytical_series';
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.normalization_area_m2 = S;
    spectrum.meta.normalization_volume_m3 = V0;
    spectrum.meta.Omega_rad_s = Omega;
    spectrum.meta.a0_m = a0;
    spectrum.meta.a0_mode = cfg.laser.a0_mode;
    spectrum.meta.series_s_max = cfg.oap.series.s_max;
    spectrum.meta.series_eta_max = cfg.oap.series.eta_max;
    spectrum.meta.D_domain = 'D > 0 (Eq. (5)/(7) derivation domain)';
    spectrum.meta.piezo_dispersion = 'literal_Eq_7_omega0';
    spectrum.meta.D_positive_point_count = 0;
    spectrum.meta.D_nonpositive_point_count = 0;

    for it = 1:numel(td)
        initial = td(it).initial;
        Theta = td(it).Q_half_inv_m;
        if ~(isscalar(Theta) && isfinite(Theta) && Theta >= 0)
            error('QW:AnalyticalFormFactor', ...
                'Theta = td.Q_half_inv_m must be finite and nonnegative.');
        end

        for ell = orders
            order_key = sprintf('order_%d', ell);
            if ell == 1
                field_factor = a0.^2;
            else
                field_factor = a0.^4 / 16;
            end

            for process = 1:2
                if process == 1
                    xi = +1;
                    bose = N0 + 1;
                    process_key = 'emission_raw';
                else
                    xi = -1;
                    bose = N0;
                    process_key = 'absorption_raw';
                end

                % Eq. (5)/(7): D = DeltaE + xi*hbar*omega0 - ell*Ephot.
                % Eq. (7) is followed literally here even though its
                % introductory coupling model discusses omega_q = c_s*q.
                % The q-dependent acoustic dispersion belongs only to
                % direct_q_integral; do not replace hw0 by hbar*c_s*qd.
                D = td(it).deltaE_J + xi*hw0 - ell*Ephot;
                positive_D = D > 0;
                spectrum.meta.D_positive_point_count = ...
                    spectrum.meta.D_positive_point_count + nnz(positive_D);
                spectrum.meta.D_nonpositive_point_count = ...
                    spectrum.meta.D_nonpositive_point_count + nnz(~positive_D);

                for im = 1:numel(cfg.oap.mechanisms)
                    mechanism = cfg.oap.mechanisms{im};
                    switch lower(mechanism)
                        case 'optical'
                            mechanism_key = 'optical';
                            Cmech = Copt;
                            Jkind = 'optical';
                        case {'piezoelectric', 'piezo'}
                            mechanism_key = 'piezoelectric';
                            Cmech = Cpiezo;
                            Jkind = 'piezoelectric';
                        otherwise
                            error('QW:AnalyticalMechanism', ...
                                'Unsupported phonon mechanism: %s', mechanism);
                    end

                    series_sum = zeros(size(Ephot));
                    compensation = zeros(size(Ephot));
                    for s = 1:cfg.oap.series.s_max
                        % Eq. (5)/(7) electron occupation is present only
                        % through this Taylor-Fermi factor.  Do not add a
                        % separate subband population weight.
                        fermi_factor = exp(-s * ...
                            (sp.E_J(initial)-sp.EF_J) / kBT);
                        if ~isfinite(fermi_factor)
                            error('QW:AnalyticalNonfiniteFermiFactor', ...
                                ['Non-finite Taylor-Fermi factor at ', ...
                                 'transition %d, s=%d.'], it, s);
                        end

                        for eta = 0:cfg.oap.series.eta_max
                            % analytical_exp_bessel_factor evaluates
                            % exp(s*D/(2*kBT))*J with scaled MATLAB besselk.
                            expJ = analytical_exp_bessel_factor( ...
                                D, s, eta, ell, Jkind, mstar, c.hbar, kBT);
                            screening_factor = (-1)^(s+eta) * (eta+1) * ...
                                qd^(2*eta);
                            term = zeros(size(Ephot));
                            term(positive_D) = screening_factor * ...
                                fermi_factor .* expJ(positive_D);
                            if any(~isfinite(term(positive_D)))
                                error('QW:AnalyticalNonfiniteTerm', ...
                                    ['Non-finite Eq. (5)/(7) term for %s, ', ...
                                     'ell=%d, s=%d, eta=%d.'], ...
                                    mechanism_key, ell, s, eta);
                            end

                            % Compensated summation reduces roundoff in the
                            % genuinely alternating s/eta series.
                            corrected = term - compensation;
                            updated = series_sum + corrected;
                            compensation = (updated-series_sum) - corrected;
                            series_sum = updated;
                        end
                    end

                    P = Cmech * common * Theta * bose .* ...
                        field_factor .* series_sum;
                    if any(~isfinite(P))
                        error('QW:AnalyticalNonfiniteSpectrum', ...
                            ['Non-finite analytical spectrum for %s, ', ...
                             'ell=%d.'], mechanism_key, ell);
                    end

                    % Do not clip negative values: they are part of the raw
                    % finite-truncation alternating series.
                    spectrum.mechanism.(mechanism_key).(order_key). ...
                        (process_key) = ...
                        spectrum.mechanism.(mechanism_key).(order_key). ...
                        (process_key) + P;
                end
            end
        end
    end

    mechanism_names = {'optical', 'piezoelectric'};
    for im = 1:numel(mechanism_names)
        mechanism_key = mechanism_names{im};
        for ell = orders
            order_key = sprintf('order_%d', ell);
            P = spectrum.mechanism.(mechanism_key).(order_key).emission_raw + ...
                spectrum.mechanism.(mechanism_key).(order_key).absorption_raw;
            spectrum.mechanism.(mechanism_key).(order_key).total_raw = P;
            spectrum.mechanism.(mechanism_key).total_raw = ...
                spectrum.mechanism.(mechanism_key).total_raw + P;
        end
        spectrum.total_raw = spectrum.total_raw + ...
            spectrum.mechanism.(mechanism_key).total_raw;
    end

    spectrum.meta.negative_raw_point_count = nnz(spectrum.total_raw < 0);
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
