function VH = solve_poisson_dirichlet(z, Nd_z, n_z, cfg)
%SOLVE_POISSON_DIRICHLET Solve d2 VH/dz2 = e^2/(eps eps0)(Nd-n).
% VH is the electron Hartree potential energy in joules.

    c = cfg.constants;
    dz = z(2) - z(1);
    rhs = (c.e^2 / (cfg.material.eps_static*c.eps0)) .* (Nd_z - n_z);

    n_in = numel(z) - 2;
    main = -2.0 * ones(n_in, 1);
    off = ones(n_in-1, 1);
    D2 = spdiags([[off;0], main, [0;off]], [-1 0 1], n_in, n_in);

    VH = zeros(size(z));
    VH(2:end-1) = D2 \ (dz^2 .* rhs(2:end-1));
end
