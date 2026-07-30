function L = generalized_laguerre_numeric(n, alpha, x)
%GENERALIZED_LAGUERRE_NUMERIC L_n^alpha(x) for nonnegative integer n.

    if n == 0
        L = ones(size(x));
        return;
    elseif n == 1
        L = 1 + alpha - x;
        return;
    end

    Lnm2 = ones(size(x));
    Lnm1 = 1 + alpha - x;
    for k = 2:n
        L = ((2*k-1+alpha-x).*Lnm1 - (k-1+alpha).*Lnm2) ./ k;
        Lnm2 = Lnm1;
        Lnm1 = L;
    end
    L = Lnm1;
end
