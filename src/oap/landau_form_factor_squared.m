function J2 = landau_form_factor_squared(qperp, cfg)
%LANDAU_FORM_FACTOR_SQUARED Eq. (12) generalized Laguerre factor.
%
% For N_i=N_f=0, J2=exp(-u), u=lB^2 qperp^2/2.
% At B=0 a Maxwell-Boltzmann recoil cutoff is used instead.

    c = cfg.constants;
    B = cfg.fields.B_T;
    Ni = cfg.oap.landau_initial;
    Nf = cfg.oap.landau_final;

    if B <= 1e-12 || strcmpi(cfg.oap.inplane_factor, 'thermal')
        mstar = cfg.material.mstar_rel*c.m0;
        exponent = -cfg.oap.thermal_cutoff_factor .* ...
            c.hbar^2 .* qperp.^2 ./ (2*mstar*c.kB*cfg.temperature_K);
        J2 = exp(exponent);
        return;
    end

    if strcmpi(cfg.oap.inplane_factor, 'none')
        J2 = ones(size(qperp));
        return;
    end

    lB = sqrt(c.hbar/(c.e*B));
    u = 0.5*lB^2.*qperp.^2;

    n1 = max(Ni, Nf);
    n2 = min(Ni, Nf);
    dN = n1 - n2;

    L = generalized_laguerre_numeric(n2, dN, u);
    log_pref = gammaln(n2+1) - gammaln(n1+1);
    J2 = exp(log_pref - u) .* (u.^dN) .* (L.^2);
end
