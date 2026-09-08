function value = direct_q_channel( ...
    Eph, eps_n, EF, DeltaE, Theta, ...
    ell, zeta, mechanism, qd, cfg)

    c = cfg.constants;

    hbar = c.hbar;
    me   = cfg.material.mstar_rel * c.m0;
    kBT  = c.kB * cfg.temperature_K;

    E0 = cfg.laser.E0_kVcm * 1e5;
    Omega = Eph / hbar;

    % Source definition
    a0 = c.e * E0 / (me * Omega^2);

    hw0 = cfg.material.LO_phonon_meV * c.meV;

    % D_{n,n',ell,zeta}
    D = DeltaE + zeta*hw0 - ell*Eph;

    q = linspace(0, ...
        cfg.oap.qperp_max_inv_nm/c.nm, ...
        cfg.oap.Nqperp);

    % Eq. before the Taylor-Bessel manipulation:
    %
    % B(q) = hbar^2 q^2/(2m) + D
    Bq = hbar^2 .* q.^2 ./ (2*me) + D;

    % Delta-function root:
    % k* = -m B(q)/(hbar^2 q)
    %
    % Only k* >= 0 belongs to the integral k in [0,inf).
    valid = q > 0 & Bq < 0;

    Ik = zeros(size(q));

    kstar = -me .* Bq(valid) ./ ...
        (hbar^2 .* q(valid));

    eps_star = eps_n + ...
        hbar^2 .* kstar.^2 ./ (2*me);

    x = (eps_star - EF) ./ kBT;

    f = stable_fermi_local(x);

    % Integral over k_perp after using delta(g(k))
    Ik(valid) = ...
        (-me^2 .* Bq(valid) ./ ...
        (hbar^4 .* q(valid).^2)) .* f;

    % Screening / electron-phonon coupling q dependence
    switch lower(mechanism)
        case 'optical'
            coupling_q = ...
                q.^2 ./ (q.^2 + qd^2).^2;

        case 'piezoelectric'
            coupling_q = ...
                q.^3 ./ (q.^2 + qd^2).^2;

        otherwise
            error('Unknown mechanism');
    end

    % (a0*q)^(2 ell)/(2^(2 ell) (ell!)^2)
    dressing = ...
        (a0 .* q).^(2*ell) ./ ...
        (2^(2*ell) * factorial(ell)^2);

    % q dq measure after q ~= q_perp
    integrand = ...
        q .* coupling_q .* dressing .* Ik;

    value = Theta .* trapz(q, integrand);

    % Phonon occupation
    N0 = 1/expm1(hw0/kBT);

    if zeta == +1
        value = (N0 + 1) .* value;
    else
        value = N0 .* value;
    end
end


function f = stable_fermi_local(x)

    f = zeros(size(x));

    pos = x >= 0;

    t = exp(-x(pos));
    f(pos) = t ./ (1+t);

    t = exp(x(~pos));
    f(~pos) = 1 ./ (1+t);
end