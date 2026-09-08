function [expJ, positive_D] = analytical_exp_bessel_factor( ...
        D, s, eta, ell, kind, mstar, hbar, kBT)
%ANALYTICAL_EXP_BESSEL_FACTOR Stable exp(sD/2kBT)*J for Eq. (5)/(7).
%
% The final Bessel formulas were derived for D > 0.  Points with D <= 0
% are therefore excluded explicitly rather than extrapolated with abs(D)
% or hidden behind an arbitrary numerical floor.  For valid points,
%
%   zeta = s*abs(D)/(2*kBT) = s*D/(2*kBT),
%
% and besselk(nu,zeta,1) evaluates exp(zeta)*K_nu(zeta), avoiding the
% overflow/underflow of forming exp(zeta) and K_nu separately.

    if any(~isfinite(D(:)))
        error('QW:AnalyticalNonfiniteD', ...
            'Analytical-series energy mismatch D must be finite.');
    end
    if ~(isscalar(s) && isfinite(s) && s >= 1 && s == round(s)) || ...
            ~(isscalar(eta) && isfinite(eta) && eta >= 0 && eta == round(eta))
        error('QW:AnalyticalSeriesIndex', ...
            's and eta must be nonnegative integer series indices with s >= 1.');
    end
    if ~(isscalar(ell) && ismember(ell, [1 2]))
        error('QW:AnalyticalPhotonOrder', ...
            'Analytical Bessel factors support ell = 1 or 2 only.');
    end

    switch lower(kind)
        case 'optical'
            power = -eta + ell;
            nu1 = eta - ell + 1;
            nu2 = eta - ell;
        case {'piezoelectric', 'piezo'}
            power = -eta + ell + 0.5;
            nu1 = eta - ell + 0.5;
            nu2 = eta - ell - 0.5;
        otherwise
            error('QW:AnalyticalBesselKind', ...
                'Unknown analytical Bessel kind: %s', kind);
    end

    positive_D = D > 0;
    expJ = zeros(size(D));
    if ~any(positive_D(:))
        return;
    end

    Dvalid = D(positive_D);
    base = 2*mstar*Dvalid/hbar^2; % [m^-2]
    zeta = s*abs(Dvalid)/(2*kBT);
    Kdiff_scaled = besselk(abs(nu1), zeta, 1) - ...
                   besselk(abs(nu2), zeta, 1);
    % Jopt has units J*m^(2*eta-2*ell+2); Jpiezo has one less
    % power of length.  Multiplication by qd^(2*eta) occurs in the caller.
    expJ(positive_D) = hbar^2/(2*mstar) .* ...
        base.^power .* Kdiff_scaled;

    if any(~isfinite(expJ(positive_D)))
        error('QW:AnalyticalNonfiniteBessel', ...
            ['Non-finite scaled Bessel factor in the valid D>0 domain ', ...
             '(kind=%s, ell=%d, s=%d, eta=%d).'], kind, ell, s, eta);
    end
end
