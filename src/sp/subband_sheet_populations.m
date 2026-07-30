function Ni = subband_sheet_populations(EF, E, cfg)
%SUBBAND_SHEET_POPULATIONS 2D Fermi-Dirac population of each subband.
    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    pref = mstar * kBT / (pi * c.hbar^2);
    Ni = pref .* softplus_stable((EF - E(:)) ./ kBT);
end
