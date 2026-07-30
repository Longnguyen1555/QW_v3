function td = compute_transition_data(sp, cfg)
%COMPUTE_TRANSITION_DATA Dipole and qz form factors for requested transitions.

    c = cfg.constants;
    qz = linspace(0, cfg.oap.qz_max_inv_nm, cfg.oap.Nqz) / c.nm;
    wqz = trapezoid_weights(qz);

    transitions = cfg.oap.transitions;
    ntr = size(transitions, 1);
    td = repmat(struct(), ntr, 1);

    phase = exp(1i * (sp.z_m * qz));

    for t = 1:ntr
        i = transitions(t,1);
        f = transitions(t,2);
        if f <= i
            error('Each transition must satisfy final > initial.');
        end

        product = conj(sp.Psi(:,i)) .* sp.Psi(:,f);
        dipole = trapz(sp.z_m, product .* sp.z_m);
        Iqz = trapz(sp.z_m, bsxfun(@times, product, phase), 1);

        td(t).initial = i;
        td(t).final = f;
        td(t).deltaE_J = sp.E_J(f) - sp.E_J(i);
        td(t).deltaE_meV = td(t).deltaE_J / c.meV;
        td(t).dipole_m = dipole;
        td(t).dipole_nm = dipole / c.nm;
        td(t).qz_inv_m = qz;
        td(t).qz_inv_nm = qz * c.nm;
        td(t).Iqz = Iqz;
        td(t).I2 = abs(Iqz).^2;
        td(t).Q_half_inv_m = sum(td(t).I2 .* wqz);
        td(t).Q_full_inv_m = 2.0 * td(t).Q_half_inv_m;
    end
end
