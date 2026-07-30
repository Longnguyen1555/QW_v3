function spectrum = initialize_spectrum_struct(Ephot, cfg)
%INITIALIZE_SPECTRUM_STRUCT Allocate contribution arrays.

    c = cfg.constants;
    nE = numel(Ephot);

    spectrum.energy_J = Ephot;
    spectrum.energy_meV = Ephot/c.meV;
    spectrum.total_raw = zeros(1,nE);
    spectrum.meta = struct();

    all_mechanisms = {'optical','piezoelectric'};
    for im = 1:numel(all_mechanisms)
        mk = all_mechanisms{im};
        spectrum.mechanism.(mk).total_raw = zeros(1,nE);
        spectrum.mechanism.(mk).transitions = struct();
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d', order);
            spectrum.mechanism.(mk).(ok).emission_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).absorption_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).total_raw = zeros(1,nE);
        end
    end
end
