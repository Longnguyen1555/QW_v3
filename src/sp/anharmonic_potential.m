function U = anharmonic_potential(z, cfg)
%ANHARMONIC_POTENTIAL Target-paper confinement potential.
% U_A(z) = U0 (z/Lz)^2 [alpha (z/Lz)^6 - 1].

    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    u = z ./ Lz;
    U = U0 .* (u.^2) .* (cfg.structure.alpha .* u.^6 - 1.0);
end
