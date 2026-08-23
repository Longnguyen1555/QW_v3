function EF = solve_fermi_level(Nd_sheet, E, cfg)
%SOLVE_FERMI_LEVEL
% Solve charge neutrality:
%
%   Nd = sum_i m*kBT/(pi*hbar^2)
%        * log(1 + exp((EF-Ei)/kBT))
%
% Root finding is performed in a DIMENSIONLESS energy variable
% to avoid numerical problems caused by solving directly in Joules.

    c = cfg.constants;

    T = cfg.temperature_K;
    if T <= 0
        error('solve_fermi_level requires T > 0 K.');
    end

    kBT = c.kB * T;

    mstar = cfg.material.mstar_rel * c.m0;

    pref = mstar * kBT / (pi * c.hbar^2);

    E = E(:);

    if Nd_sheet <= 0
        EF = min(E) - 50*kBT;
        return;
    end

    % ============================================================
    % Dimensionless formulation
    %
    % mu = (EF - Eref)/kBT
    % ============================================================

    Eref = min(E);

    eps_i = (E - Eref) ./ kBT;

    target = Nd_sheet / pref;

    neutrality = @(mu) ...
        sum(softplus_stable(mu - eps_i)) - target;

    % ============================================================
    % Initial bracket in dimensionless units
    % ============================================================

    lo = -50.0;

    % At low temperature target can become large, therefore
    % construct a sufficiently large upper bound.
    hi = max(50.0, max(eps_i) + target + 20.0);

    flo = neutrality(lo);
    fhi = neutrality(hi);

    % Expand lower side if necessary
    count = 0;
    while flo > 0 && count < 50
        lo = lo - 50;
        flo = neutrality(lo);
        count = count + 1;
    end

    % Expand upper side if necessary
    count = 0;
    while fhi < 0 && count < 50
        hi = hi + 50;
        fhi = neutrality(hi);
        count = count + 1;
    end

    if flo > 0 || fhi < 0
        error('Could not bracket dimensionless Fermi level.');
    end

    % ============================================================
    % Root solve
    % ============================================================

    opts = optimset( ...
        'TolX', 1.0e-12, ...
        'FunValCheck', 'on', ...
        'Display', 'off');

    muF = fzero(neutrality, [lo hi], opts);

    % Return SI energy
    EF = Eref + muF * kBT;

    % ============================================================
    % Mandatory neutrality verification
    % ============================================================

    Ni = pref .* ...
        softplus_stable((EF - E) ./ kBT);

    rel_err = abs(sum(Ni) - Nd_sheet) / Nd_sheet;

    if rel_err > 1.0e-10
        error(['Fermi-level neutrality failure: ', ...
               'relative error = %.3e'], rel_err);
    end
end