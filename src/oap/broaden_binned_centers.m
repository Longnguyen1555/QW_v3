function shape = broaden_binned_centers( ...
    Egrid, centers, weights, order, gamma)
%BROADEN_BINNED_CENTERS Lorentzian broadening on a uniform energy grid.
%
% If
%
%   X = DeltaE +/- hbar*omega - ell*E
%
% then
%
%   L(X;Gamma)
%     = (1/ell) L(E-Ec; Gamma/ell).
%
% Therefore Gamma/order is the HWHM measured on the photon-energy axis.
%
% NOTE:
% Gamma passed to this routine is phenomenological unless a separate,
% physically derived lifetime model is supplied.

    Egrid = Egrid(:).';

    if numel(Egrid) < 2
        error('Egrid must contain at least two points.');
    end

    if ~(isscalar(order) && ...
         isfinite(order) && ...
         order >= 1 && ...
         order == round(order))
        error('order must be a positive integer.');
    end

    if ~(isscalar(gamma) && ...
         isfinite(gamma) && ...
         gamma >= 0)
        error('gamma must be a finite non-negative scalar.');
    end

    d = diff(Egrid);
    dE = d(1);

    tol = 1e-10 * max(abs(dE), realmin);

    if any(abs(d-dE) > tol)
        error('Egrid must be uniformly spaced.');
    end

    line_density = zeros(size(Egrid));

    centers = centers(:);
    weights = real(weights(:));

    weights(~isfinite(weights)) = 0;

    bin = round((centers-Egrid(1))/dE) + 1;

    valid = ...
        bin >= 1 & ...
        bin <= numel(Egrid) & ...
        weights ~= 0 & ...
        isfinite(centers);

    bin = bin(valid);

    values = weights(valid) ./ (order*dE);

    if ~isempty(bin)
        accumulated = accumarray( ...
            bin, ...
            values, ...
            [numel(Egrid),1], ...
            @sum, ...
            0);

        line_density = accumulated.';
    end

    % Zero-width limit: discrete representation of the delta function.
    if gamma == 0
        shape = line_density;
        return;
    end

    n = numel(Egrid);

    offset = ...
        ((1:n)-ceil((n+1)/2))*dE;

    gamma_E = gamma/order;

    kernel = ...
        gamma_E/pi ./ ...
        (offset.^2 + gamma_E^2);

    shape = ...
        conv(line_density, kernel, 'same') .* dE;
end