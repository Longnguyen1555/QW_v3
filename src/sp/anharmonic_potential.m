function U = anharmonic_potential(z, cfg)
%ANHARMONIC_POTENTIAL Manning confinement in PDF Eq. (1).
% U_A(z) = nu*U0*[sech(z/Lz)^4 - sech(z/Lz)^2].

    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    nu = cfg.structure.manning_prefactor;
    u = z ./ Lz;
    U = U0 .* nu .* (sech(u).^4 - sech(u).^2);
end
