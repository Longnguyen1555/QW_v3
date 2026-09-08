function [weight, kstar, B, f_initial] = direct_k_delta_weight( ...
    Einitial, EF, deltaE, Eph, ell, zeta, hw_q, qperp, cfg)
%DIRECT_K_DELTA_WEIGHT Exact radial k integral of the energy Dirac delta.
%
% The source uses the collinear convention k_perp dot q_perp = k_perp*q_perp:
%
%   B = hbar^2*q_perp^2/(2*m*) + deltaE + zeta*hw_q - ell*Eph,
%   kstar = -m*B/(hbar^2*q_perp).
%
% Since the radial integral is over k_perp in [0,Inf), B must be <= 0 and
% kstar must be nonnegative.  The exact result is
%
%   m*kstar/(hbar^2*q_perp) * f_n(kstar).
%
% The q_perp=0 endpoint is assigned zero.  In the complete radial q
% integral its measure and coupling make this single endpoint vanish; this
% also avoids evaluating an undefined 0/0 or Inf expression.

    if ~(isscalar(ell) && any(ell == [1 2]))
        error('QW:DirectPhotonOrder', 'ell must be 1 or 2.');
    end
    if ~(isscalar(zeta) && any(zeta == [-1 1]))
        error('QW:DirectPhononProcess', 'zeta must be +1 or -1.');
    end

    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    if ~(isfinite(kBT) && kBT > 0)
        error('QW:DirectTemperature', 'Temperature must be positive.');
    end

    qperp = qperp(:);
    Eph = Eph(:).';
    if isscalar(hw_q)
        hw_q = hw_q .* ones(size(qperp));
    else
        hw_q = hw_q(:);
    end
    if numel(hw_q) ~= numel(qperp)
        error('QW:DirectGridSize', ...
            'hw_q must be scalar or have one value per q_perp point.');
    end

    recoil = c.hbar^2 .* qperp.^2 ./ (2*mstar);
    B = recoil + deltaE + zeta.*hw_q - ell.*Eph;
    qgrid = qperp + zeros(size(B));

    weight = zeros(size(B));
    kstar = zeros(size(B));
    f_initial = zeros(size(B));

    active = (qgrid > 0) & (B <= 0);
    kstar(active) = -mstar .* B(active) ./ ...
        (c.hbar^2 .* qgrid(active));

    kinetic_initial = c.hbar^2 .* kstar(active).^2 ./ (2*mstar);
    x_initial = (Einitial + kinetic_initial - EF) ./ kBT;
    f = stable_fermi_dirac(x_initial);

    if isfield(cfg.oap, 'include_pauli_blocking') && ...
            cfg.oap.include_pauli_blocking
        q_active = qgrid(active);
        Efinal = Einitial + deltaE;
        kinetic_final = c.hbar^2 .* (kstar(active) + q_active).^2 ./ ...
                        (2*mstar);
        f_final = stable_fermi_dirac((Efinal + kinetic_final - EF) ./ kBT);
        f = f .* (1.0 - f_final);
    end

    f_initial(active) = f;
    jacobian_weight = mstar .* kstar(active) ./ ...
                      (c.hbar^2 .* qgrid(active));

    % If f underflows to exactly zero, assign the product's limiting zero
    % without forming Inf*0 at extreme, nonphysical grid points.
    product = zeros(size(jacobian_weight));
    finite_mask = isfinite(jacobian_weight) & (f > 0);
    product(finite_mask) = jacobian_weight(finite_mask) .* f(finite_mask);
    weight(active) = product;
end

function f = stable_fermi_dirac(x)
%STABLE_FERMI_DIRAC Evaluate 1/(1+exp(x)) without exponential overflow.
    f = zeros(size(x));
    positive = (x >= 0);
    t = exp(-x(positive));
    f(positive) = t ./ (1.0 + t);
    t = exp(x(~positive));
    f(~positive) = 1.0 ./ (1.0 + t);
end
