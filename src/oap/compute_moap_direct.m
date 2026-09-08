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
    [WP, WZ] = meshgrid( ...
        wperp, ...
        trapezoid_weights(td(1).qz_inv_m));

    q = sqrt(QP.^2 + QZ.^2);

    % Cylindrical integral:
    % int_{-inf}^{inf}dqz int_0^inf 2*pi*qperp*dqperp /(2*pi)^3
    measure = (2.0/(2*pi)^2) .* QP .* WP .* WZ;

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
        

        if strcmpi(mechanism, 'optical')
   
            mech_key = 'optical';
        else
           
            mech_key = 'piezoelectric';
        end

        mechanism_total = zeros(1,nE);

        for it = 1:numel(td)

            i = td(it).initial;

            DeltaE = td(it).deltaE_J;
            eps_n  = sp.E_J(i);
            EF     = sp.EF_J;
            Theta  = td(it).Q_half_inv_m;

            for ell = cfg.oap.photon_orders

                ok = sprintf('order_%d', ell);

                for iE = 1:nE

                    Eph = Ephot(iE);

                    P_em = direct_q_channel( ...
                        Eph, eps_n, EF, DeltaE, Theta, ...
                        ell, +1, 'optical', qd, cfg);

                    P_ab = direct_q_channel( ...
                        Eph, eps_n, EF, DeltaE, Theta, ...
                        ell, -1, 'optical', qd, cfg);

                    spectrum.mechanism.optical.(ok).emission_raw(iE) = ...
                        spectrum.mechanism.optical.(ok).emission_raw(iE) + P_em;

                    spectrum.mechanism.optical.(ok).absorption_raw(iE) = ...
                        spectrum.mechanism.optical.(ok).absorption_raw(iE) + P_ab;

                    spectrum.mechanism.optical.(ok).total_raw(iE) = ...
                        spectrum.mechanism.optical.(ok).total_raw(iE) ...
                        + P_em + P_ab;

                end
            end
        end

        spectrum.mechanism.(mech_key).total_raw = mechanism_total;
        spectrum.total_raw = spectrum.total_raw + mechanism_total;
    end

    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
