function U = anharmonic_potential(z, cfg)
%ANHARMONIC_POTENTIAL Target-paper confinement potential.
% U_A(z) = U0 (z/Lz)^2 [alpha (z/Lz)^6 - 1].

    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    u = z ./ Lz;
    U = U0 .* 6.0 .* (sech(u).^4 - sech(u).^2);
end
