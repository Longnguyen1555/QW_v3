function spectrum = compute_moap_analytical_series(sp, td, cfg)
%COMPUTE_MOAP_ANALYTICAL_SERIES Analytical Taylor-Bessel Eq. (5)/(7).
%
% This implementation follows the supplied calculation note literally:
%   - optical mechanism: Eq. (5);
%   - piezoelectric mechanism: Eq. (7);
%   - photon orders ell = 1,2 only;
%   - the same D_{n,n',ell,zeta} = DeltaE + zeta*hbar*omega0 - ell*hbar*Omega
%     and N0 are used in both Eq. (5) and Eq. (7), exactly as written in
%     the supplied analytical formulas.
%
% The series is truncated at cfg.oap.series.s_max and eta_max. Because the
% Fermi Taylor series used in the derivation requires E_n > E_F, this
% condition is checked when cfg.oap.series.enforce_fermi_condition is true.

    c = cfg.constants;
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    spectrum = initialize_spectrum_struct(Ephot, cfg);

    [qd, ne3d] = debye_wavevector(sp, cfg);
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd*c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;

    E0 = cfg.laser.E0_kVcm * 1e5;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    hw0 = cfg.material.LO_phonon_meV * c.meV;
    N0 = 1/expm1(hw0/kBT);

    % a0 = e E0 / (m* Omega^2), unless a constant value is explicitly
    % requested for a controlled comparison/debug run.
    switch lower(cfg.laser.a0_mode)
        case 'constant'
            a0 = (cfg.laser.a0_nm*c.nm) .* ones(size(Ephot));
        case 'dynamic'
            omega = Ephot/c.hbar;
            a0 = c.e*E0 ./ (mstar*omega.^2);
        otherwise
            error('Unknown laser.a0_mode: %s', cfg.laser.a0_mode);
    end
    spectrum.meta.a0_m = a0;
    spectrum.meta.a0_nm = a0/c.nm;

    % S and V are normalization quantities in Eq. (5)/(7). Do not replace
    % S by Lz^2: Lz is a confinement length, not the in-plane normalization
    % area. The project already exposes normalization_area_m2 for this role.
    S = cfg.oap.normalization_area_m2;
    thickness = cfg.oap.active_thickness_nm * c.nm;
    if ~(isscalar(S) && isfinite(S) && S > 0)
        error('cfg.oap.normalization_area_m2 must be a positive finite scalar.');
    end
    if ~(isscalar(thickness) && isfinite(thickness) && thickness > 0)
        error('cfg.oap.active_thickness_nm must be positive.');
    end
    V = S * thickness;

    % Eq. (5) optical prefactor C and Eq. (7) piezoelectric prefactor C_pi.
    Copt = E0^2*sqrt(cfg.material.eps_static)/8 * ...
        c.e^2*(hw0/c.hbar) / (V*c.eps0) * ...
        (1/cfg.material.eps_high - 1/cfg.material.eps_static);

    Cpiezo = cfg.material.piezo_kappa2*c.e^2*E0^2 * ...
        cfg.material.sound_speed_mps*sqrt(cfg.material.eps_static) / ...
        (8*V*cfg.material.eps_static*c.eps0);

    % Common factor outside braces in Eq. (5)/(7). a0 may be energy-dependent.
    common = a0.^2 .* (S*V*mstar^2/(32*pi^3*c.hbar^4));

    % The geometric expansion used for the Fermi distribution in the note
    % converges only when exp[-(E_n-E_F)/(kBT)] < 1, i.e. E_n > E_F.
    

    supported_orders = cfg.oap.photon_orders(ismember(cfg.oap.photon_orders,[1 2]));
    unsupported_orders = cfg.oap.photon_orders(~ismember(cfg.oap.photon_orders,[1 2]));
    if ~isempty(unsupported_orders)
        warning('analytical_series implements Eq. (5)/(7) only for photon orders 1 and 2; ignoring %s.', ...
                mat2str(unsupported_orders));
    end
    if isempty(supported_orders)
        error('analytical_series requires photon order 1 and/or 2.');
    end

    for it = 1:numel(td)
        i = td(it).initial;
        Q = td(it).Q_half_inv_m;  % Theta_{n'',n} = int_0^inf |I(qz)|^2 dqz
        boltz0 = sp.E_J(i) - sp.EF_J;

        for ell = supported_orders
            Dplus  = td(it).deltaE_J + hw0 - ell*Ephot; % zeta = +1, emission
            Dminus = td(it).deltaE_J - hw0 - ell*Ephot; % zeta = -1, absorption

            for im = 1:numel(cfg.oap.mechanisms)
                mechanism = lower(cfg.oap.mechanisms{im});
                switch mechanism
                    case 'optical'
                        Cmech = Copt;
                        Jtype = 'optical';
                        mk = 'optical';
                    case {'piezoelectric','piezo'}
                        Cmech = Cpiezo;
                        Jtype = 'piezoelectric';
                        mk = 'piezoelectric';
                    otherwise
                        error('Unsupported analytical-series mechanism: %s', mechanism);
                end

                series_em = zeros(size(Ephot));
                series_ab = zeros(size(Ephot));

                for s = 1:cfg.oap.series.s_max
                    boltz = exp(-s*boltz0/kBT);
                    for eta = 0:cfg.oap.series.eta_max
                        sign_series = (-1)^(s+eta) * (eta+1) * boltz;

                        % Evaluate qd^(2eta) * exp[sD/(2kBT)] * J_{ell,zeta}
                        % as one stable product. MATLAB's scaled besselk(...,1)
                        % avoids the exp(+x)*K_nu(x) overflow/cancellation.
                        Jem = weighted_J_exp(Dplus,  eta, ell, Jtype, s, qd, mstar, c.hbar, kBT);
                        Jab = weighted_J_exp(Dminus, eta, ell, Jtype, s, qd, mstar, c.hbar, kBT);

                        term_em = sign_series .* Jem;
                        term_ab = sign_series .* Jab;

                        if any(~isfinite(term_em)) || any(~isfinite(term_ab))
                            error(['Non-finite analytical-series term encountered. This usually means ' ...
                                'the photon-energy grid hits D=0 exactly or the truncated Taylor-Bessel ' ...
                                'series is numerically divergent. Do not silently replace this term by zero.']);
                        end

                        series_em = series_em + term_em;
                        series_ab = series_ab + term_ab;
                    end
                end

                if ell == 1
                    order_factor = ones(size(a0));
                else
                    order_factor = a0.^2/16;
                end

                P_em = Cmech .* common .* Q .* (N0+1) .* order_factor .* series_em;
                P_ab = Cmech .* common .* Q .* N0     .* order_factor .* series_ab;

                % Do not clip negative values here. A negative truncated result
                % is a convergence/validity diagnostic and must remain visible.
                ok = sprintf('order_%d',ell);
                spectrum.mechanism.(mk).(ok).emission_raw = ...
                    spectrum.mechanism.(mk).(ok).emission_raw + P_em;
                spectrum.mechanism.(mk).(ok).absorption_raw = ...
                    spectrum.mechanism.(mk).(ok).absorption_raw + P_ab;
            end
        end
    end

    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        for ell = supported_orders
            ok = sprintf('order_%d',ell);
            P = spectrum.mechanism.(mk).(ok).emission_raw + ...
                spectrum.mechanism.(mk).(ok).absorption_raw;
            spectrum.mechanism.(mk).(ok).total_raw = P;
            spectrum.mechanism.(mk).total_raw = ...
                spectrum.mechanism.(mk).total_raw + P;
        end
        spectrum.total_raw = spectrum.total_raw + ...
            spectrum.mechanism.(mk).total_raw;
    end

    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end

function value = weighted_J_exp(D, eta, ell, kind, s, qd, mstar, hbar, kBT)
%WEIGHTED_J_EXP Compute qd^(2eta)*exp[sD/(2kBT)]*J stably.

    absD = abs(D);
    x = s*absD/(2*kBT);
    base = 2*mstar*absD/hbar^2;

    if strcmpi(kind,'optical')
        power0 = ell;
        nu1 = eta-ell+1;
        nu2 = eta-ell;
    else
        power0 = ell+0.5;
        nu1 = eta-ell+0.5;
        nu2 = eta-ell-0.5;
    end

    value = NaN(size(D));
    nz = absD > 0;
    if ~any(nz)
        return;
    end

    xn = x(nz);
    bn = base(nz);

    % besselk(nu,x,1) = exp(x)*K_nu(x) for real positive x.
    Kdiff_scaled = besselk(abs(nu1),xn,1) - besselk(abs(nu2),xn,1);

    % exp[sD/(2kBT)] K_nu(x):
    % D>0 -> exp(+x)K = scaled K
    % D<0 -> exp(-x)K = exp(-2x)*scaled K
    exp_correction = ones(size(xn));
    localD = D(nz);
    neg = localD < 0;
    exp_correction(neg) = exp(-2*xn(neg));

    % qd^(2eta) * base^(-eta+power0)
    % = base^power0 * (qd^2/base)^eta.
    value(nz) = hbar^2/(2*mstar) .* ...
        bn.^power0 .* (qd^2./bn).^eta .* ...
        exp_correction .* Kdiff_scaled;
end
