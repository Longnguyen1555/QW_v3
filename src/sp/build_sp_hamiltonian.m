function H = build_sp_hamiltonian(z_in, dz, Vconf_in, VH_in, cfg)
%BUILD_SP_HAMILTONIAN Finite-difference Eq. (5)-(6) Hamiltonian.
%
% Ueff = (hbar*kx + e*B*z)^2/(2m*) + UA + UH - e*E*z.

    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    B = cfg.fields.B_T;
    Efield = cfg.fields.E_kVcm * 1.0e5;

    alpha_kin = c.hbar^2 / (2.0 * mstar);

   
    magnetic = (c.e^2 * B^2 .* z_in.^2) ./ (2.0*mstar);
    electric = -c.e * Efield .* z_in;
    Veff = Vconf_in + VH_in + magnetic + electric;

    n = numel(z_in);
    off = (-alpha_kin / dz^2) * ones(n-1, 1);
    main = (2.0*alpha_kin/dz^2) + Veff(:);

    H = spdiags([[off; 0], main, [0; off]], [-1 0 1], n, n);
    H = (H + H') / 2;
end
