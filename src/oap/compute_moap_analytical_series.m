function spectrum = compute_moap_analytical_series(sp, td, cfg)
%COMPUTE_MOAP_ANALYTICAL_SERIES Truncated Eq. (5)/(7) implementation.
%
% This mode follows the Taylor-Maclaurin/Bessel-K expressions in the
% calculation note. The direct_q_integral mode is recommended because the
% series contains alternating, dimensionful terms and can converge slowly.
%
% The implementation is provided for formula-by-formula comparison, not as
% the default production kernel.

    c = cfg.constants;
    Ephot = cfg.oap.photon_energy_meV(:).'*c.meV;
    spectrum = initialize_spectrum_struct(Ephot, cfg);
    [qd, ne3d] = debye_wavevector(sp, cfg);
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.screening_ne_m3 = ne3d;

    E0 = cfg.laser.E0_kVcm*1e5;
    S = (cfg.structure.Lz_nm*c.nm)^2;
    V = S*(cfg.oap.active_thickness_nm*c.nm);
    mstar = cfg.material.mstar_rel*c.m0;
    a0 = cfg.laser.a0_nm*c.nm;

    populations = electron_populations_for_oap(sp, cfg);
    pop_norm = populations/max(sum(populations), realmin);

    N0 = 1/expm1(cfg.material.LO_phonon_meV*c.meV/...
                 (c.kB*cfg.temperature_K));

    Copt = E0^2*sqrt(cfg.material.eps_static)/8 * ...
        c.e^2*(cfg.material.LO_phonon_meV*c.meV/c.hbar) / ...
        (V*c.eps0) * ...
        (1/cfg.material.eps_high-1/cfg.material.eps_static);

    Cpiezo = cfg.material.piezo_kappa2*c.e^2*E0^2*...
        cfg.material.sound_speed_mps*sqrt(cfg.material.eps_static) / ...
        (8*V*cfg.material.eps_static*c.eps0);

    common = a0^2*S*V*mstar^2/(32*pi^3*c.hbar^4);

    for it = 1:numel(td)
        i = td(it).initial;
        Q = td(it).Q_half_inv_m;
        for order = cfg.oap.photon_orders
            if order > 2
                continue; % Eq. (5)/(7) note explicitly lists l=1,2.
            end
            for process = 1:2
                if process == 1
                    sign_ph = +1; bose = N0+1; pname = 'emission_raw';
                else
                    sign_ph = -1; bose = N0; pname = 'absorption_raw';
                end

                Dopt = td(it).deltaE_J + sign_ph*...
                       cfg.material.LO_phonon_meV*c.meV - order*Ephot;

                for im = 1:numel(cfg.oap.mechanisms)
                    mechanism = cfg.oap.mechanisms{im};
                    if strcmpi(mechanism,'optical')
                        Cmech = Copt;
                        Jtype = 'optical';
                        D = Dopt;
                        mk = 'optical';
                    else
                        Cmech = Cpiezo;
                        Jtype = 'piezoelectric';
                        qeff = qd;
                        hweff = c.hbar*cfg.material.sound_speed_mps*qeff;
                        D = td(it).deltaE_J + sign_ph*hweff - order*Ephot;
                        mk = 'piezoelectric';
                    end

                    series_sum = zeros(size(Ephot));
                    for s = 1:cfg.oap.series.s_max
                        boltz = exp(-s*(sp.E_J(i)-sp.EF_J)/...
                                   (c.kB*cfg.temperature_K));
                        for eta = 0:cfg.oap.series.eta_max
                            x = s*abs(D)/(2*c.kB*cfg.temperature_K);
                            J = analytical_J(D, x, eta, order, Jtype, ...
                                             mstar, c.hbar);
                            term = (-1)^(s+eta)*(eta+1)*...
                                   boltz*qd^(2*eta) .* ...
                                   exp(s*D/(2*c.kB*cfg.temperature_K)) .* J;
                            term(~isfinite(term)) = 0;
                            series_sum = series_sum + term;
                        end
                    end

                    order_factor = 1;
                    if order == 2
                        order_factor = a0^2/16;
                    end

                    P = Cmech*common*Q*pop_norm(i)*bose*...
                        order_factor.*series_sum;
                    P(P<0) = 0;

                    ok = sprintf('order_%d',order);
                    spectrum.mechanism.(mk).(ok).(pname) = ...
                        spectrum.mechanism.(mk).(ok).(pname) + P;
                end
            end
        end
    end

    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
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

function J = analytical_J(D, x, eta, ell, kind, mstar, hbar)
    base = 2*mstar*max(abs(D),realmin)/hbar^2;

    if strcmpi(kind,'optical')
        exponent = -eta + ell;
        nu1 = eta-ell+1;
        nu2 = eta-ell;
    else
        exponent = -eta + ell + 0.5;
        nu1 = eta-ell+0.5;
        nu2 = eta-ell-0.5;
    end

    Kdiff = besselk(abs(nu1), max(x,1e-12)) - ...
            besselk(abs(nu2), max(x,1e-12));
    J = hbar^2/(2*mstar) .* base.^exponent .* Kdiff;
end
