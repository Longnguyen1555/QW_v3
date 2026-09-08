function spectrum = compute_moap_direct(sp, td, cfg)
%COMPUTE_MOAP_DIRECT
% Direct q_perp integration following the calculation note.
%
% No phenomenological Lorentzian broadening.
% FWHM is extracted afterwards from the calculated spectrum.

    c = cfg.constants;

    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);

    if any(Ephot <= 0)
        error(['direct_q_integral requires hbar*Omega > 0 because ', ...
               'a0 = eE0/(m*Omega^2).']);
    end

    [qd, ne3d] = debye_wavevector(sp, cfg);

    spectrum = initialize_spectrum_struct(Ephot, cfg);

    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd*c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.model = 'direct_q_integral';

    mechanisms = cfg.oap.mechanisms;

    if ischar(mechanisms)
        mechanisms = {mechanisms};
    end

    for im = 1:numel(mechanisms)

        mechanism = lower(mechanisms{im});

        switch mechanism
            case 'optical'
                mech_key = 'optical';

            case {'piezoelectric','piezo'}
                mech_key = 'piezoelectric';

            otherwise
                error('Unsupported mechanism: %s', mechanism);
        end

        mechanism_total = zeros(1,nE);

        for it = 1:numel(td)

            i = td(it).initial;

            DeltaE = td(it).deltaE_J;
            eps_n  = sp.E_J(i);
            EF     = sp.EF_J;

            % Theta_nn' = integral |I_nn'(qz)|^2 dqz
            Theta = td(it).Q_half_inv_m;

            transition_total = zeros(1,nE);

            for ell = cfg.oap.photon_orders

                if ~ismember(ell,[1 2])
                    error('Source derivation is restricted to ell = 1,2.');
                end

                ok = sprintf('order_%d',ell);

                for iE = 1:nE

                    Eph = Ephot(iE);

                    P_em = direct_q_channel( ...
                        Eph, eps_n, EF, DeltaE, Theta, ...
                        ell, +1, mechanism, qd, cfg);

                    P_ab = direct_q_channel( ...
                        Eph, eps_n, EF, DeltaE, Theta, ...
                        ell, -1, mechanism, qd, cfg);

                    spectrum.mechanism.(mech_key).(ok). ...
                        emission_raw(iE) = ...
                        spectrum.mechanism.(mech_key).(ok). ...
                        emission_raw(iE) + P_em;

                    spectrum.mechanism.(mech_key).(ok). ...
                        absorption_raw(iE) = ...
                        spectrum.mechanism.(mech_key).(ok). ...
                        absorption_raw(iE) + P_ab;

                    spectrum.mechanism.(mech_key).(ok). ...
                        total_raw(iE) = ...
                        spectrum.mechanism.(mech_key).(ok). ...
                        total_raw(iE) + P_em + P_ab;

                    transition_total(iE) = ...
                        transition_total(iE) + P_em + P_ab;
                end
            end

            tr_key = sprintf( ...
                'transition_%d_%d', ...
                td(it).initial, td(it).final);

            spectrum.mechanism.(mech_key). ...
                transitions.(tr_key) = transition_total;

            mechanism_total = ...
                mechanism_total + transition_total;
        end

        spectrum.mechanism.(mech_key).total_raw = ...
            mechanism_total;

        spectrum.total_raw = ...
            spectrum.total_raw + mechanism_total;
    end

    spectrum = normalize_spectrum_for_plot(spectrum,cfg);
end