function U = anharmonic_potential(z, cfg)
%ANHARMONIC_POTENTIAL Manning/sech confinement potential.
% U_A(z) = U0 * p * (sech(z/Lz)^4 - sech(z/Lz)^2), where p is the
% configured Manning prefactor.  The retired polynomial alpha parameter
% is deliberately not used by this model.

    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    u = z ./ Lz;
    p = cfg.structure.manning_prefactor;
    U = U0 .* p .* (sech(u).^4 - sech(u).^2);
end
