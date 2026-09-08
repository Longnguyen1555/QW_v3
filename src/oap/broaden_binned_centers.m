function [shape, binned_shape] = broaden_binned_centers( ...
    Egrid, centers, weights, order, gamma, mode)
%BROADEN_BINNED_CENTERS Bin q-space weights and optionally broaden them.
%
% binned_shape is the direct numerical delta-function representation,
% including the 1/order Jacobian.  In legacy_lorentzian mode, shape is the
% historical convolution
%
%   L(Delta + hw - l*E; gamma) = (1/l) L(E-Ecenter; gamma/l).

    if nargin < 6 || isempty(mode)
        mode = 'legacy_lorentzian';
    end

    Egrid = Egrid(:).';
    dE = Egrid(2)-Egrid(1);
    line_density = zeros(size(Egrid));

    centers = centers(:);
    weights = weights(:);
    bin = round((centers-Egrid(1))/dE) + 1;

    valid = bin >= 1 & bin <= numel(Egrid) & ...
            isfinite(weights) & weights ~= 0 & isfinite(centers);
    bin = bin(valid);
    values = weights(valid) ./ (order*dE);

    if ~isempty(bin)
        accumulated = accumarray(bin, values, [numel(Egrid),1], @sum, 0);
        line_density = accumulated.';
    end

    binned_shape = line_density;

    switch lower(mode)
        case 'none'
            shape = binned_shape;

        case 'legacy_lorentzian'
            if ~isscalar(gamma) || ~isfinite(gamma) || gamma <= 0
                error('Legacy Lorentzian broadening requires gamma > 0.');
            end
            n = numel(Egrid);
            offset = ((1:n)-ceil((n+1)/2))*dE;
            gamma_E = gamma/order;
            kernel = gamma_E/pi ./ (offset.^2 + gamma_E^2);
            shape = conv(binned_shape, kernel, 'same') .* dE;

        otherwise
            error('Unknown oap.broadening.mode: %s', mode);
    end
end
