function [T, info] = compute_refined_T_pdf(eps_alpha, EF, D, ell, qd, cfg)
%COMPUTE_REFINED_T_PDF Eq. (28), or the truncated Eq. (32) checked against it.
% auto: use series only when it agrees with quadrature; series: fail otherwise.
% Agreement at finite truncation is not proof of convergence of the infinite sum.

    c = cfg.constants;
    m = cfg.material.mstar_rel*c.m0;
    kBT = c.kB*cfg.temperature_K;
    mode = 'auto';
    if isfield(cfg.oap.series, 'evaluation'), mode = lower(cfg.oap.series.evaluation); end
    assert(ismember(mode, {'auto','series','quadrature'}), 'Unknown series.evaluation.');
    validateattributes(eps_alpha, {'numeric'}, {'scalar','real','finite'});
    validateattributes(EF, {'numeric'}, {'scalar','real','finite'});
    validateattributes(D, {'numeric'}, {'real','finite'});
    validateattributes(qd, {'numeric'}, {'scalar','real','finite','positive'});
    assert(isscalar(ell) && ismember(ell,[1 2]), 'Photon order must be 1 or 2.');
    assert(kBT > 0 && m > 0, 'Temperature and effective mass must be positive.');

    shape = size(D);
    D = D(:).';
    strict = strcmp(mode, 'series');
    if strict && eps_alpha <= EF
        error('QW:Refined:FermiDomain', ...
            ['Eq. (32) requires E_initial > EF (PDF p.7). E_initial=%.9g meV, ' ...
             'EF=%.9g meV, (E_initial-EF)/kBT=%.6g. Use evaluation=''auto'' ' ...
             'for Eq. (28); increasing s_max/v_max cannot fix this condition.'], ...
            eps_alpha/c.meV, EF/c.meV, (eps_alpha-EF)/kBT);
    end
    if strict && any(D == 0)
        error('QW:Refined:ZeroD', ...
            'Eq. (32) is singular at D=0. Use evaluation=''auto'' for Eq. (28).');
    end

    T = parent_integral(eps_alpha-EF, D, ell, qd, cfg.oap.Nqperp, m, c.hbar, kBT);
    candidate = (eps_alpha > EF) & (D ~= 0) & ~strcmp(mode,'quadrature');
    used = false(size(D));
    if any(candidate)
        smax = cfg.oap.series.s_max;
        vmax = cfg.oap.series.v_max;
        validateattributes(smax, {'numeric'}, {'scalar','finite','integer','nonnegative'});
        validateattributes(vmax, {'numeric'}, {'scalar','finite','integer','positive'});
        S = bessel_sum(eps_alpha-EF, D(candidate), ell, qd, smax, vmax, m, c.hbar, kBT);
        reference = T(candidate);
        relative_error = abs(S-reference)./max(abs(reference), realmin);
        ok = isfinite(S) & S >= 0 & relative_error <= 1e-3;
        if strict && any(~ok)
            error('QW:Refined:SeriesMismatch', ...
                ['Eq. (32) at s_max=%d, v_max=%d disagrees with Eq. (28) ' ...
                 '(relative error > 1e-3 or non-finite/negative sum). ' ...
                 'Use evaluation=''auto''; a finite sum need not be converged.'], smax, vmax);
        end
        indices = find(candidate);
        T(indices(ok)) = S(ok);
        used(indices(ok)) = true;
    end
    assert(all(isfinite(T) & T >= 0), 'Invalid refined-PDF integral.');
    info = struct('series_points', nnz(used), 'quadrature_points', nnz(~used));
    T = reshape(T, shape);
end

function T = bessel_sum(gap, D, ell, qd, smax, vmax, m, hbar, kBT)
% Page 7: (s+1), not the misprinted (v+1). Eq. (32): SUM of the two K terms.
    Eqd = hbar^2*qd^2/(2*m);
    T = zeros(size(D));
    for s = 0:smax
        prefactor = qd^(2*ell-2)*Eqd.*(Eqd./abs(D)).^(s-ell);
        for v = 1:vmax
            x = v*abs(D)/(2*kBT);
            K = besselk(abs(s-ell),x,1) + besselk(abs(s-ell+1),x,1);
            % exp[-v(gap+D/2)/kBT] * exp(-x); avoids exp(+x) overflow.
            weight = exp(-v*(gap+max(D,0))/kBT);
            T = T + (-1)^(s+v-1)*(s+1).*prefactor.*weight.*K;
        end
    end
end

function T = parent_integral(gap, D, ell, qd, Nq, m, hbar, kBT)
% Eq. (28): Gauss-Legendre on q=qscale*t/(1-t), 0<t<1.
% Retain the PDF's Eq+abs(D), which differs from abs(Eq+D) when D<0.
    n = max(96,2*Nq);
    beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
    [V,L] = eig(diag(beta,1)+diag(beta,-1), 'vector');
    t = (L+1)/2;
    w = V(1,:).'.^2;
    qscale = max(qd,sqrt(8*m*kBT)/hbar);
    q = qscale*t./(1-t);
    Eq = hbar^2*q.^2/(2*m);
    w = w.*qscale./(1-t).^2 .* q.^(2*ell+1)./(q.^2+qd^2).^2;
    T = zeros(size(D));
    for first = 1:512:numel(D)
        index = first:min(first+511,numel(D));
        d = D(index);
        % Square form of Eq. (30) keeps the kinetic energy nonnegative.
        x = (gap+(Eq+d).^2./(4*Eq))/kBT;
        F = exp(-max(x,0))./(1+exp(-abs(x)));
        T(index) = sum(w.*(Eq+abs(d)).*F,1);
    end
end
