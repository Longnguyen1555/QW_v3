function w = trapezoid_weights(x)
%TRAPEZOID_WEIGHTS Weights for a uniform or nonuniform 1D trapezoidal rule.
    x = x(:).';
    n = numel(x);
    if n < 2
        error('At least two grid points are required.');
    end
    w = zeros(size(x));
    w(1) = (x(2)-x(1))/2;
    w(end) = (x(end)-x(end-1))/2;
    if n > 2
        w(2:end-1) = (x(3:end)-x(1:end-2))/2;
    end
end
