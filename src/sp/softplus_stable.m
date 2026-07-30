function y = softplus_stable(x)
%SOFTPLUS_STABLE log(1+exp(x)) without overflow.
    y = max(x, 0) + log1p(exp(-abs(x)));
end
