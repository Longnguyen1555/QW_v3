function k_parallel = inplane_delta_root(qperp, balance_J, cfg)
%INPLANE_DELTA_ROOT Root selected by the in-plane energy delta function.
%
% With q_perp along the parallel axis,
% Delta epsilon = DeltaE + hbar^2*q_perp^2/(2*m*) +
%                 hbar^2*k_parallel*q_perp/m*.
% balance_J is ell*hbar*Omega - DeltaE -/+ hbar*omega_q, so solving
% Delta epsilon - DeltaE = balance_J gives the root below.

    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    recoil = c.hbar^2 .* qperp.^2 ./ (2*mstar);
    denominator = c.hbar^2 .* qperp ./ mstar;

    k_parallel = NaN(size(qperp));
    nonzero = qperp > 0;
    k_parallel(nonzero) = (balance_J(nonzero) - recoil(nonzero)) ./ ...
        denominator(nonzero);
end
