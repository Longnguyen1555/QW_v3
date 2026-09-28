function test_refined_analytical_series()
%TEST_REFINED_ANALYTICAL_SERIES Refined PDF series/integral regression.

    cfg = default_config();
    cfg.oap.Nqperp = 180;
    c = cfg.constants;

    % In the stated qd << q and eps_alpha > EF domain, Eq. (32) must
    % reproduce its unexpanded parent integral Eq. (20)/(28).
    eps_alpha = 50*c.meV;
    EF = 0;
    qd = 0.01/c.nm;
    D = [-120 -80 80 120]*c.meV;

    for ell = [1 2]
        cfg.oap.series.evaluation = 'series';
        [Tseries, used] = refined_T_from_spectrum( ...
            eps_alpha, EF, D, ell, qd, cfg);
        assert(used.series_points == numel(D) && used.quadrature_points == 0, ...
            'Forced series silently substituted quadrature.');
        cfg.oap.series.evaluation = 'quadrature';
        Tquadrature = refined_T_from_spectrum( ...
            eps_alpha, EF, D, ell, qd, cfg);
        relative_error = abs(Tseries-Tquadrature)./Tquadrature;
        assert(all(relative_error < 1.0e-10), ...
            'Refined Eq. (32) does not match Eq. (20)/(28).');
    end

    % Specific errors replace the old combined heuristic-domain message.
    cfg.oap.series.evaluation = 'series';
    expect_error(@() refined_T_from_spectrum(-20*c.meV,0,D,1,qd,cfg), ...
        'QW:Refined:FermiDomain');
    expect_error(@() refined_T_from_spectrum(eps_alpha,0,0,1,qd,cfg), ...
        'QW:Refined:ZeroD');
    cfg.oap.series.s_max = 0;
    cfg.oap.series.v_max = 1;
    expect_error(@() refined_T_from_spectrum(5*c.meV,0,-80*c.meV,1,qd,cfg), ...
        'QW:Refined:SeriesMismatch');
    cfg.oap.series.s_max = 8;
    cfg.oap.series.v_max = 20;
    expect_error(@() refined_T_from_spectrum( ...
        eps_alpha,0,1e-3*c.meV,1,0.2/c.nm,cfg), ...
        'QW:Refined:SeriesMismatch');

    % The old implementation diverged when eps_alpha <= EF and at D=0.
    % Auto mode must use the exact parent integral and remain nonnegative.
    cfg.oap.series.evaluation = 'auto';
    grid = [-80 0 1e-3 80]*c.meV;
    [mixed, info] = refined_T_from_spectrum(eps_alpha,0,grid,1,qd,cfg);
    assert(isequal(size(mixed),size(grid)) && all(isfinite(mixed(:))), ...
        'Mixed series/quadrature evaluation lost the spectrum shape.');
    assert(info.series_points == 2 && info.quadrature_points == 2, ...
        'Auto must retain quadrature at D=0 and for the divergent screening sum.');
    [Tfallback, info] = refined_T_from_spectrum( ...
        -20*c.meV, 0, [-40 0 40]*c.meV, 1, 0.2/c.nm, cfg);
    assert(all(isfinite(Tfallback)) && all(Tfallback >= 0), ...
        'Refined parent-integral fallback is invalid.');
    assert(info.series_points == 0 && info.quadrature_points == 3, ...
        'Invalid Fermi-series points were not sent to quadrature.');

    % Independent adaptive quadrature verifies the fixed-node integral,
    % including exact D=0 and both signs of D, for both photon orders.
    for ell = [1 2]
        for d = [-40 0 40]*c.meV
            actual = refined_T_from_spectrum( ...
                -20*c.meV,0,d,ell,0.2/c.nm,cfg);
            expected = reference_integral(-20*c.meV,d,ell,0.2/c.nm,cfg);
            assert(abs(actual-expected)/expected < 1e-8, ...
                'Eq. (28) disagrees with independent adaptive quadrature.');
        end
    end

    % End-to-end synthetic spectrum checks Eq. (22)/(26) assembly,
    % including |Q_ab|^2*K_ab and the combined acoustic contribution.
    cfg.oap.model = 'van';
    cfg.oap.photon_energy_meV = [20 30 33.46 50 66.92 80];
    cfg.oap.photon_orders = [1 2];
    cfg.oap.mechanisms = {'optical','piezoelectric'};
    cfg.oap.Nqperp = 48;

    sp = struct();
    sp.E_J = [50 80]*c.meV;
    sp.EF_J = 0;
    sp.Nd_sheet_m2 = 1.0e17;

    td = struct();
    td.initial = 1;
    td.final = 2;
    td.deltaE_J = 30*c.meV;
    td.deltaE_meV = 30;
    td.dipole_m = 1.0*c.nm;
    td.Q_half_inv_m = 1.0e8;

    spectrum = compute_moap_spectrum(sp, td, cfg);
    assert(all(isfinite(spectrum.total_raw)), ...
        'Refined analytical spectrum contains non-finite values.');
    assert(all(spectrum.total_raw >= 0) && max(spectrum.total_raw) > 0, ...
        'Refined analytical spectrum is negative or zero.');
    assert(all(spectrum.mechanism.piezoelectric.order_1.absorption_raw == 0), ...
        'Eq. (26) acoustic term was incorrectly split into two processes.');
    assert(spectrum.meta.piezo_processes_combined, ...
        'Combined acoustic-process metadata is missing.');

    doubled = td;
    doubled.dipole_m = 2*td.dipole_m;
    scaled = compute_moap_spectrum(sp,doubled,cfg);
    assert(max(abs(scaled.total_raw-4*spectrum.total_raw)) == 0, ...
        'Dipole factor must scale as its absolute square.');
    two = compute_moap_spectrum(sp,[td td],cfg);
    assert(max(abs(two.total_raw-2*spectrum.total_raw)) == 0, ...
        'Multiple transitions must add exactly once.');

    negative_cfg = cfg;
    negative_cfg.oap.model = 'van';
    negative_cfg.oap.photon_energy_meV = [-80 -40 -20];
    negative_cfg.laser.a0_mode = 'dynamic';
    negative = compute_moap_van(sp,td,negative_cfg);
    assert(isequal(negative.energy_meV, ...
        negative_cfg.oap.photon_energy_meV) && ...
        all(isfinite(negative.total_raw)), ...
        'Van lost or invalidated the negative photon-energy branch.');
    zero_cfg = negative_cfg;
    zero_cfg.oap.photon_energy_meV = [-1 0 1];
    expect_error(@() compute_moap_van(sp,td,zero_cfg), ...
        'QW:VanPhotonEnergyZero');

    fprintf('  PASS: refined analytical series\n');
end

function [T, info] = refined_T_from_spectrum(eps_alpha,EF,D,ell,qd,cfg)
% Recover the private refined integral through the public Van spectrum API.
    c = cfg.constants;
    D = D(:).';
    transition = 256*c.meV;
    cfg.oap.photon_energy_meV = (transition-D)/(ell*c.meV);
    cfg.oap.photon_orders = ell;
    cfg.oap.mechanisms = {'piezoelectric'};
    cfg.oap.plot_normalization = 'none';
    cfg.laser.a0_mode = 'constant';

    ne = qd^2*cfg.material.eps_static*c.eps0 * ...
        c.kB*cfg.temperature_K/c.e^2;
    cfg.oap.screening_density_mode = 'fixed';
    cfg.oap.fixed_ne_cm3 = ne/1e6;

    sp = struct('E_J',[eps_alpha eps_alpha+transition], ...
        'EF_J',EF,'Nd_sheet_m2',0);
    td = struct('initial',1,'final',2, ...
        'deltaE_J',transition,'deltaE_meV',transition/c.meV, ...
        'dipole_m',c.nm,'Q_half_inv_m',1.0e8);

    spectrum = compute_moap_van(sp,td,cfg);
    factor = abs(td.dipole_m)^2*td.Q_half_inv_m * ...
        (cfg.laser.a0_nm*c.nm)^(2*ell) / ...
        (2^(2*ell)*factorial(ell)^2);
    key = sprintf('order_%d',ell);
    raw = spectrum.mechanism.piezoelectric.(key).emission_raw;
    T = raw./(spectrum.meta.P0_piezoelectric*factor);
    info = spectrum.meta.series_diagnostic;
end

function expect_error(f, identifier)
    try
        f();
    catch err
        assert(strcmp(err.identifier,identifier), 'Unexpected error: %s',err.message);
        return;
    end
    error('Expected error was not raised: %s',identifier);
end

function value = reference_integral(gap,D,ell,qd,cfg)
% Dimensionless u=q/qthermal prevents the SI-scale integral from underflowing.
    c = cfg.constants;
    kBT = c.kB*cfg.temperature_K;
    qs = sqrt(8*cfg.material.mstar_rel*c.m0*kBT)/c.hbar;
    d = D/kBT;
    f = @(u) u.^(2*ell+1)./(u.^2+(qd/qs)^2).^2 .* ...
        (4*u.^2+abs(d))./(1+exp(gap/kBT+u.^2+d/2+d^2./(16*u.^2)));
    value = qs^(2*ell-2)*kBT*integral(f,0,Inf,'RelTol',1e-10,'AbsTol',0);
end
