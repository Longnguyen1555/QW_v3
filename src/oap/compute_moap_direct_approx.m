function spectrum = compute_moap_direct_approx(sp, td, cfg)
%COMPUTE_MOAP_DIRECT_APPROX Direct q-space quadrature of Eqs. (19)/(25).

    c = cfg.constants;

    hbar = c.hbar;
    e    = c.e;
    kBT  = c.kB * cfg.temperature_K;
    eps0 = c.eps0;
    me = cfg.material.mstar_rel * c.m0;

    kappa0   = cfg.material.eps_static;
    kappa_LO = 1.0 / (1.0/cfg.material.eps_high - 1.0/kappa0);

    hw0 = cfg.material.LO_phonon_meV * c.meV;
    N0 = 1.0 / expm1(hw0 / kBT);
    E0 = cfg.laser.E0_kVcm * 1.0e5;
    S = cfg.oap.normalization_area_m2;

    orders = cfg.oap.photon_orders(:).';
    orders = orders(orders <= 2); % Match Nhan's 1PA/2PA limit.
    use_optical = any(strcmpi(cfg.oap.mechanisms,'optical'));
    use_piezoelectric = any(strcmpi( ...
        cfg.oap.mechanisms,'piezoelectric'));

    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);

    qd = debye_wavevector(sp, cfg);
    qperp = linspace(0,cfg.oap.qperp_max_inv_nm, ...
        cfg.oap.Nqperp)/c.nm;
    qz = td(1).qz_inv_m(:);
    I2 = reshape([td.I2],numel(qz),[]);
    a0 = cfg.laser.a0_nm*c.nm;

    spectrum = initialize_spectrum_struct(Ephot, cfg);

    for iE = 1:nE

        Eph = Ephot(iE);

        Omega = Eph / hbar;
        P0_LO = E0^4 * e^4 * me^2 * (hw0/hbar) * S * sqrt(kappa0) / (256.0*pi^3*hbar^6*eps0*kappa_LO*Omega^2);
        P0_PE = E0^4 * e^4 * me^2 * cfg.material.piezo_kappa2 * kBT * S * sqrt(kappa0) / (256.0 * pi^3 * hbar^7 * kappa0 * eps0 * Omega);

        for it = 1:numel(td)

            DeltaE = td(it).deltaE_J;

            transition_total_LO = 0.0;
            transition_total_PE = 0.0;

            for ell = orders

                order_key = sprintf('order_%d', ell);

                order_factor = a0^(2*ell) / (2^(2*ell) * factorial(ell)^2);
                if use_optical
                    Pplus = P0_LO*abs(td(it).dipole_m)^2*order_factor*(N0+1.0)* ...
                        compute_T_qspace( sp.E_J(td(it).initial),sp.EF_J,DeltaE+hw0-ell*Eph,ell,qd,qperp,qz,I2(:,it),cfg);

                    Pminus = P0_LO*abs(td(it).dipole_m)^2*order_factor*N0* ...
                        compute_T_qspace( ...
                        sp.E_J(td(it).initial),sp.EF_J, ...
                        DeltaE-hw0-ell*Eph,ell,qd, ...
                        qperp,qz,I2(:,it),cfg);

                    spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                    spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                    spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                    transition_total_LO = transition_total_LO+Pplus+Pminus;
                end

                if use_piezoelectric
                    Ppiezo = P0_PE*abs(td(it).dipole_m)^2* ...
                        order_factor*compute_T_qspace( ...
                        sp.E_J(td(it).initial),sp.EF_J, ...
                        DeltaE-ell*Eph,ell,qd, ...
                        qperp,qz,I2(:,it),cfg);
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

function Tval = compute_T_qspace(eps_alpha,EF,D,ell,qd,qperp,qz, ...
        I2_qz,cfg)
% Composite trapz in physical q space: qz first, then qperp.
    c = cfg.constants;
    me = cfg.material.mstar_rel*c.m0;
    E_qperp = c.hbar^2*qperp.^2/(2*me);
    epsilon_star = eps_alpha + D/2 + (me * D^2)./(2*c.hbar^2*qperp.^2);
    occupation = 1./(1+exp((epsilon_star-EF)/(c.kB*cfg.temperature_K)));
    q2 = qperp.^2+qz.^2;
    integrand = q2.^(ell+0.5)./(q2+qd^2).^2 .* abs(E_qperp+D).*occupation.*I2_qz;
    Tval = trapz(qperp,trapz(qz,integrand,1),2);
end
