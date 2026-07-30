function shape = broaden_binned_centers(Egrid, centers, weights, order, gamma)
%BROADEN_BINNED_CENTERS Efficient Lorentzian sum on a uniform energy grid.
%
% L(Delta + hw - l*E; gamma)
% = (1/l) L(E-Ecenter; gamma/l).

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

    n = numel(Egrid);
    offset = ((1:n)-ceil((n+1)/2))*dE;
    gamma_E = gamma/order;
    kernel = gamma_E/pi ./ (offset.^2 + gamma_E^2);

    shape = conv(line_density, kernel, 'same') .* dE;
end
