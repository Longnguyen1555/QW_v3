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
%
% Centers inside Egrid are binned and broadened by convolution. Centers
% outside Egrid are NOT discarded: their Lorentzian tails are evaluated
% directly on Egrid and added to the result. This is necessary because a
% Lorentzian has infinite support.

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

    if numel(centers) ~= numel(weights)
        error('centers and weights must contain the same number of elements.');
    end

    weights(~isfinite(weights)) = 0;

    finite_nonzero = ...
        isfinite(centers) & ...
        weights ~= 0;

    % Centers physically inside the requested photon-energy window are
    % treated with the original fast bin-and-convolve scheme.
    inside = ...
        finite_nonzero & ...
        centers >= Egrid(1) & ...
        centers <= Egrid(end);

    centers_in = centers(inside);
    weights_in = weights(inside);

    bin = round((centers_in-Egrid(1))/dE) + 1;

    % Guard against floating-point round-off at the two endpoints.
    bin = max(1, min(numel(Egrid), bin));

    % delta(Delta - order*E) = (1/order) delta(E-Ec), and a discrete delta
    % on a uniform grid contributes 1/dE.
    values = weights_in ./ (order*dE);

    if ~isempty(bin)
        accumulated = accumarray( ...
            bin, ...
            values, ...
            [numel(Egrid),1], ...
            @sum, ...
            0);

        line_density = accumulated.';
    end

    % Zero-width limit: centers outside Egrid have no support inside Egrid,
    % so only the discrete in-window delta representation remains.
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

    % Contribution from centers that lie inside Egrid.
    shape = ...
        conv(line_density, kernel, 'same') .* dE;

    % A Lorentzian has infinite support. Therefore a center outside Egrid
    % still contributes a non-zero tail inside the requested window. Add
    % those contributions directly, using the exact (unbinned) center:
    %
    %   weight * (1/pi) * Gamma /
    %       ([order*(E-Ec)]^2 + Gamma^2)
    %
    % This is exactly equivalent to
    %
    %   (weight/order) * L(E-Ec; Gamma/order).
    outside = finite_nonzero & ~inside;

    centers_out = centers(outside);
    weights_out = weights(outside);

    for k = 1:numel(centers_out)
        detuning = order * (Egrid - centers_out(k));
        shape = shape + ...
            weights_out(k) * (gamma/pi) ./ ...
            (detuning.^2 + gamma^2);
    end
end
