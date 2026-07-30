function EF = solve_fermi_level(Nd_sheet, E, cfg)
%SOLVE_FERMI_LEVEL Charge-neutrality solution sum_i N_i(EF)=Nd_sheet.

    c = cfg.constants;
    kBT = c.kB * cfg.temperature_K;

    if Nd_sheet <= 0
        EF = min(E) - 100*kBT;
        return;
    end

    f = @(x) sum(subband_sheet_populations(x, E, cfg)) - Nd_sheet;
    lo = min(E) - 100*kBT;
    hi = max(E) + 100*kBT;

    flo = f(lo);
    fhi = f(hi);
    expand_count = 0;
    while flo*fhi > 0 && expand_count < 20
        lo = lo - 50*kBT;
        hi = hi + 50*kBT;
        flo = f(lo);
        fhi = f(hi);
        expand_count = expand_count + 1;
    end

    if flo*fhi > 0
        error('Could not bracket the Fermi level.');
    end

    EF = fzero(f, [lo hi]);
end
