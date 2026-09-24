function spectrum = compute_moap_analytical_series(sp, td, cfg)
%COMPUTE_MOAP_ANALYTICAL_SERIES Refined PDF Eqs. (22), (26), (28), (32).
% Formula conventions and source inconsistencies: ANALYTICAL_SERIES_REFINED_AUDIT.md.

    c = cfg.constants;
    E = cfg.oap.photon_energy_meV(:).'*c.meV;
    validateattributes(E, {'numeric'}, {'vector','real','finite','positive','nonempty'});
    orders = cfg.oap.photon_orders(:).';
    assert(~isempty(orders) && all(ismember(orders,[1 2])) && ...
        numel(unique(orders)) == numel(orders), 'Use unique photon orders from [1 2].');
    mechs = cellfun(@lower,cfg.oap.mechanisms,'UniformOutput',false);
    mechs(strcmp(mechs,'piezo')) = {'piezoelectric'};
    mechs = unique(mechs,'stable');
    assert(all(ismember(mechs,{'optical','piezoelectric'})), 'Unknown mechanism.');
    spectrum = initialize_spectrum_struct(E,cfg);

    m = cfg.material.mstar_rel*c.m0;
    E0 = cfg.laser.E0_kVcm*1e5;
    kBT = c.kB*cfg.temperature_K;
    hw0 = cfg.material.LO_phonon_meV*c.meV;
    Omega = E/c.hbar;
    N0 = 1/expm1(hw0/kBT);
    [qd, ne] = debye_wavevector(sp,cfg);
    S = cfg.oap.normalization_area_m2;
    validateattributes(S, {'numeric'}, {'scalar','real','finite','positive'});
    switch lower(cfg.laser.a0_mode)
        case 'constant'
            a0 = cfg.laser.a0_nm*c.nm*ones(size(E));
        case 'dynamic'
            a0 = c.e*E0./(m*Omega.^2);
        otherwise
            error('Unknown laser.a0_mode: %s', cfg.laser.a0_mode);
    end

    % Prefactors derived from Eqs. (2),(4),(5),(7)/(23),(11),(18).
    epss = cfg.material.eps_static;
    common = E0^4*c.e^4*m^2*S*sqrt(epss)./(256*pi^3*c.hbar^7*c.eps0*Omega.^2);
    P0.optical = common*hw0*(1/cfg.material.eps_high-1/epss);
    P0.piezoelectric = common*cfg.material.piezo_kappa2*kBT/epss;
    mode = 'auto';
    if isfield(cfg.oap.series,'evaluation'), mode = lower(cfg.oap.series.evaluation); end
    spectrum.meta = struct('model','analytical_series', ...
        'source_pdf','Refined Tinh_toan_chi_tiet.pdf', ...
        'evaluation',mode, 'series_evaluation_requested',mode, ...
        'qd_inv_m',qd, 'qd_inv_nm',qd*c.nm, 'screening_ne_m3',ne, ...
        'normalization_area_m2',S, 'a0_m',a0, 'a0_nm',a0/c.nm, ...
        'Omega_rad_s',Omega, 'P0_optical',P0.optical, ...
        'P0_piezoelectric',P0.piezoelectric, 'piezo_processes_combined',true);
    spectrum.meta.piezo_storage = 'Combined Eq. (26) in emission_raw; absorption_raw=0';
    counts = [0 0];

    for it = 1:numel(td)
        transition = abs(td(it).dipole_m)^2*td(it).Q_half_inv_m;
        for ell = orders
            key = sprintf('order_%d',ell);
            factor = transition*a0.^(2*ell)/(2^(2*ell)*factorial(ell)^2);
            for im = 1:numel(mechs)
                mk = mechs{im};
                shifts = 0; bose = 1;
                if strcmp(mk,'optical'), shifts = [hw0 -hw0]; bose = [N0+1 N0]; end
                part = spectrum.mechanism.(mk).(key);
                processes = {'emission_raw','absorption_raw'};
                for ip = 1:numel(shifts)
                    D = td(it).deltaE_J + shifts(ip) - ell*E;
                    [T,info] = compute_refined_T_pdf(sp.E_J(td(it).initial),sp.EF_J,D,ell,qd,cfg);
                    P = P0.(mk).*factor*bose(ip).*T;
                    name = processes{ip};
                    part.(name) = part.(name) + P;
                    counts = counts + [info.series_points info.quadrature_points];
                end
                part.total_raw = part.emission_raw + part.absorption_raw;
                spectrum.mechanism.(mk).(key) = part;
            end
        end
    end
    for im = 1:numel(mechs)
        mk = mechs{im};
        for ell = orders
            key = sprintf('order_%d',ell);
            spectrum.mechanism.(mk).total_raw = spectrum.mechanism.(mk).total_raw + ...
                spectrum.mechanism.(mk).(key).total_raw;
        end
        spectrum.total_raw = spectrum.total_raw + spectrum.mechanism.(mk).total_raw;
    end
    spectrum.meta.series_diagnostic = struct('series_points',counts(1), ...
        'quadrature_points',counts(2));
    spectrum = normalize_spectrum_for_plot(spectrum,cfg);
end
