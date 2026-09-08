function density = inplane_delta_fermi_density(qperp, k_parallel, Ei_J, EF_J, cfg)
%INPLANE_DELTA_FERMI_DENSITY Remaining in-plane integral after delta root.
%
% This evaluates
%   int dk_y/(2*pi)^2 f(Ei+hbar^2(k_parallel^2+k_y^2)/(2m*))
%       * m*/(hbar^2*q_perp),
% where the final factor is the exact Jacobian from integrating the energy
% delta function over k_parallel.  The finite k_y interval ends 40 kBT
% (configurable) into the Fermi tail; no collision width is introduced.

    density = zeros(size(qperp));
    if ~isfinite(EF_J)
        return;
    end

    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    A = c.hbar^2 / (2*mstar);

    valid = qperp > 0 & isfinite(k_parallel);
    if ~any(valid(:))
        return;
    end

    qv = qperp(valid);
    kv = k_parallel(valid);
    Eparallel = Ei_J + A .* kv.^2;
    ky_max = sqrt((max(EF_J - Eparallel, 0) + ...
        cfg.oap.fermi_tail_kBT*kBT) ./ A);

    t = linspace(0, 1, cfg.oap.Nkin_perp).';
    ky = t * ky_max(:).';
    x = (Eparallel(:).' + A .* ky.^2 - EF_J) ./ kBT;
    f = stable_fermi(x);
    line_integral = 2.0 .* ky_max(:).' .* trapz(t, f, 1);

    values = (mstar ./ (c.hbar^2 .* qv(:).')) .* ...
        line_integral ./ (2*pi)^2;
    density(valid) = reshape(values, size(qv));
end

function f = stable_fermi(x)
    f = zeros(size(x));
    positive = x >= 0;
    ex = exp(-x(positive));
    f(positive) = ex ./ (1 + ex);
    ex = exp(x(~positive));
    f(~positive) = 1 ./ (1 + ex);
end
