function results = moap2(cfg)
    close all;
    if nargin == 0 || isempty(cfg)
        cfg = default_config();
    elseif nargin ~= 1 || ~isstruct(cfg)
        error('moap2 expects zero arguments or one cfg struct.');
    end
    validate_config(cfg);
    if ~exist(cfg.output.directory, 'dir')
        mkdir(cfg.output.directory);
    end
    results.config = cfg;
    results.single = run_single_case(cfg);
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'project_results.mat'), 'results', 'cfg', '-v7.3');
    end
    fprintf('\n[MOAP] Completed successfully.\n');
    fprintf('[MOAP] Results folder: %s\n', cfg.output.directory);
end
function cfg = default_config()
    cfg.run.verbose = true;
    cfg.material.mstar_rel = 0.067;
    cfg.material.eps_static = 13.18;
    cfg.material.eps_high = 10.89;
    cfg.material.LO_phonon_meV = 36.25;
    cfg.material.piezo_kappa2 = 0.006;
    cfg.structure.U0_meV = 228.0;
    cfg.structure.Lz_nm = 15.0;
    cfg.structure.manning_prefactor = 6.0;
    cfg.structure.domain_nm = 60.0;
    cfg.structure.auto_domain = false;
    cfg.structure.domain_factor = 3.0;
    cfg.structure.Nz = 1502;
    cfg.structure.n_states = 6;
    cfg.doping.Nd_sheet_cm2 = 1.0e11;
    cfg.doping.width_nm = 2.0;
    cfg.fields.B_T = 5.0;   
    cfg.fields.E_kVcm = 25.0;
    cfg.fields.kx_inv_m = 0.0;
    cfg.temperature_K = 50.0;
    cfg.sp.max_outer = 100;
    cfg.sp.mix = 0.40;
    cfg.sp.tol_VH_relative = 1.0e-3;
    cfg.sp.tol_EF_meV = 1.0e-3;
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_nm = 7.5;
    cfg.oap.transitions = [1 2];
    cfg.oap.photon_orders = [1 2];
    cfg.oap.photon_energy_meV = 5.0:0.01:150.0;
    cfg.oap.qz_max_inv_nm = 2.0;
    cfg.oap.Nqz = 101;
    cfg.oap.screening_density_mode = 'from_sheet';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.active_thickness_nm = 10.0;
    cfg.oap.plot_normalization = 'none';
    cfg.oap.plot_full_energy_range = true;
    cfg.oap.series.s_max = 8;
    cfg.oap.series.v_max = 20;
    cfg.oap.series.enforce_fermi_condition = false;
    cfg.output.directory = fullfile(fileparts(mfilename('fullpath')), 'results_van');
    cfg.output.save_figures = true;
    cfg.output.save_csv = true;
    cfg.output.save_mat = true;
    cfg.output.visible = 'on';
    cfg.output.figure_format = 'png';
    cfg.output.figure_dpi = 200;
    cfg.constants.e = 1.602176634e-19;
    cfg.constants.hbar = 1.054571817e-34;
    cfg.constants.kB = 1.380649e-23;
    cfg.constants.eps0 = 8.8541878128e-12;
    cfg.constants.m0 = 9.1093837015e-31;
    cfg.constants.meV = cfg.constants.e * 1.0e-3;
    cfg.constants.nm = 1.0e-9;
end
function validate_config(cfg)
    assert(cfg.structure.Nz >= 50, 'Nz must be at least 50.');
    assert(cfg.structure.n_states >= 2, 'At least two states are required.');
    assert(cfg.structure.Lz_nm > 0, 'Lz must be positive.');
    assert(cfg.structure.U0_meV > 0, 'U0 must be positive.');
    assert(cfg.temperature_K > 0, 'Temperature must be positive.');
    assert(cfg.doping.Nd_sheet_cm2 >= 0, 'Sheet density cannot be negative.');
    assert(cfg.oap.Nqz >= 8, 'qz grid is too small.');
    photon_energy_meV = cfg.oap.photon_energy_meV;
    assert(isnumeric(photon_energy_meV) && isvector(photon_energy_meV) && ~isempty(photon_energy_meV) && isreal(photon_energy_meV) && all(isfinite(photon_energy_meV)), 'Photon-energy grid must be a nonempty finite real vector.');
    assert(all(diff(photon_energy_meV) > 0), 'Photon-energy grid must be strictly increasing.');
    assert(size(cfg.oap.transitions, 2) == 2, 'cfg.oap.transitions must have two columns [initial final].');
    assert(all(cfg.oap.transitions(:) >= 1), 'Transition indices are MATLAB 1-based and must be >= 1.');
    assert(max(cfg.oap.transitions(:)) <= cfg.structure.n_states, 'A transition index exceeds structure.n_states.');
    assert(all(cfg.oap.photon_orders >= 1), 'Photon orders must be positive integers.');
end
function case_result = run_single_case(cfg)
    if cfg.run.verbose
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf( '[MOAP] B = %.4g T, E = %.4g kV/cm, T = %.4g K\n', cfg.fields.B_T, cfg.fields.E_kVcm, cfg.temperature_K);
        fprintf( '[MOAP] Lz = %.4g nm, Nd = %.4g cm^-2\n', cfg.structure.Lz_nm, cfg.doping.Nd_sheet_cm2);
    end
    sp = solve_schrodinger_poisson(cfg);
    if ~sp.converged
        warning( 'Schrodinger-Poisson solver reached max_iter without convergence.');
    end
    transition_data = compute_transition_data(sp, cfg);
    spectrum = compute_moap_van(sp, transition_data, cfg);
    metrics = extract_spectrum_metrics(spectrum, cfg);
    case_result.cfg = cfg;
    case_result.sp = sp;
    case_result.transitions = transition_data;
    case_result.spectrum = spectrum;
    case_result.metrics = metrics;
    fig_structure = plot_electronic_structure(sp, cfg);
    fig_spectrum_optical = plot_moap_spectrum(spectrum, cfg, 'optical');
    fig_spectrum_piezoelectric = plot_moap_spectrum(spectrum, cfg, 'piezoelectric');
    fig_contributions_optical = plot_moap_contributions(spectrum, cfg, 'optical');
    fig_contributions_piezoelectric = plot_moap_contributions(spectrum, cfg, 'piezoelectric');
    figures = { fig_structure, fig_spectrum_optical, fig_spectrum_piezoelectric, fig_contributions_optical, fig_contributions_piezoelectric };
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'single.mat'), 'case_result', '-v7.3');
    end
    if cfg.output.save_csv
        T = table(spectrum.energy_meV(:), spectrum.total_raw(:), spectrum.mechanism.optical.total_raw(:), spectrum.mechanism.piezoelectric.total_raw(:), 'VariableNames', {'PhotonEnergy_meV', 'TotalRaw', 'OpticalRaw', 'PiezoelectricRaw'});
        writetable(T, fullfile(cfg.output.directory, 'single_spectrum.csv'));
        write_metrics_csv(metrics, cfg, fullfile(cfg.output.directory, 'single_metrics.csv'));
    end
    if cfg.output.save_figures
        names = {'structure', 'spectrum_optical', 'spectrum_piezoelectric', 'contributions_optical', 'contributions_piezoelectric'};
        for i = 1:min(numel(figures), numel(names))
            if ishghandle(figures{i})
                filename = fullfile(cfg.output.directory, ['single_', names{i}, '.', cfg.output.figure_format]);
                if exist('exportgraphics','file') == 2
                    exportgraphics(figures{i}, filename, 'Resolution', cfg.output.figure_dpi);
                else
                    set(figures{i},'PaperPositionMode','auto');
                    print(figures{i}, filename, '-dpng', sprintf('-r%d',cfg.output.figure_dpi));
                end
            end
        end
    end
    if cfg.run.verbose
        fprintf('[MOAP] Converged: %d in %d iterations\n', sp.converged, sp.iterations);
        fprintf('[MOAP] EF = %.8f meV\n', sp.EF_meV);
        fprintf('[MOAP] Subband energies (meV): ');
        fprintf('%.8f ', sp.E_meV);
        fprintf('\n');
        fprintf('[MOAP] Max raw OAP = %.6e\n', max(spectrum.total_raw));
    end
end
function out = solve_schrodinger_poisson(cfg)
    c = cfg.constants;
    if cfg.structure.auto_domain
        domain_nm = max(cfg.structure.domain_nm, cfg.structure.domain_factor * cfg.structure.Lz_nm);
    else
        domain_nm = cfg.structure.domain_nm;
    end
    z = linspace(-domain_nm/2, domain_nm/2, cfg.structure.Nz).' * c.nm;
    dz = z(2) - z(1);
    z_in = z(2:end-1);
    Vconf = anharmonic_potential(z, cfg);
    Nd_sheet = cfg.doping.Nd_sheet_cm2 * 1.0e4;
    width = cfg.doping.width_nm * c.nm;
    Nd_z = delta_doping_profile(z, Nd_sheet, width);
    n_states = cfg.structure.n_states;
    if isfield(cfg.sp, 'max_iter')
        max_iter = cfg.sp.max_iter;
    else
        max_iter = cfg.sp.max_outer;
    end
    if isfield(cfg.sp, 'tol_VH_relative')
        tol_VH_relative = cfg.sp.tol_VH_relative;
    else
        tol_VH_relative = 1.0e-3;
    end
    VH_old = zeros(size(z));
    EF_old = NaN;
    converged = false;
    dVH_meV = Inf;
    dVH_relative = Inf;
    dEF_meV = Inf;
    for it = 1:max_iter
        H = build_sp_hamiltonian(z_in, dz, Vconf(2:end-1), VH_old(2:end-1), cfg);
        [Psi_in, E] = lowest_eigenpairs(H, n_states);
        Psi_in = normalize_wavefunctions(Psi_in, dz);
        Psi = zeros(numel(z), n_states);
        Psi(2:end-1, :) = Psi_in;
        if Nd_sheet > 0
            EF = solve_fermi_level(Nd_sheet, E, cfg);
            Ni = subband_sheet_populations(EF, E, cfg);
        else
            EF = NaN;
            Ni = zeros(n_states,1);
            dEF_meV = 0;
        end
        n_z = (abs(Psi).^2) * Ni(:);
        VH_new = solve_poisson_dirichlet(z, Nd_z, n_z, cfg);
        VH_next = cfg.sp.mix .* VH_new + (1.0-cfg.sp.mix) .* VH_old;
        dVH_meV = max(abs(VH_next - VH_old)) / c.meV;
        dVH_relative = max(abs(VH_next - VH_old)) / max(max(abs(VH_old)), c.meV);
        if isfinite(EF_old) && isfinite(EF)
            dEF_meV = abs(EF - EF_old) / c.meV;
        end
        if dVH_relative <= tol_VH_relative && dEF_meV <= cfg.sp.tol_EF_meV
            converged = true;
            VH_old = VH_next;
            break;
        end
        VH_old = VH_next;
        EF_old = EF;
    end
    H = build_sp_hamiltonian(z_in, dz, Vconf(2:end-1), VH_old(2:end-1), cfg);
    [Psi_in, E] = lowest_eigenpairs(H, n_states);
    Psi_in = normalize_wavefunctions(Psi_in, dz);
    Psi = zeros(numel(z), n_states);
    Psi(2:end-1, :) = Psi_in;
    if Nd_sheet > 0
        EF = solve_fermi_level(Nd_sheet, E, cfg);
        Ni = subband_sheet_populations(EF, E, cfg);
    else
        EF = NaN;
        Ni = zeros(n_states,1);
    end
    n_z = (abs(Psi).^2) * Ni(:);
    Vmag = (c.hbar*cfg.fields.kx_inv_m + c.e*cfg.fields.B_T.*z).^2/(2*cfg.material.mstar_rel*c.m0);
    Velec = -c.e*cfg.fields.E_kVcm*1e5.*z;
    Veff = Vconf + VH_old + Vmag + Velec;
    out.z_m = z;
    out.z_nm = z/c.nm;
    out.dz_m = dz;
    out.Vconf_J = Vconf;
    out.VH_J = VH_old;
    out.Veff_J = Veff;
    out.Vconf_meV = Vconf/c.meV;
    out.VH_meV = VH_old/c.meV;
    out.Veff_meV = Veff/c.meV;
    out.E_J = E(:);
    out.E_meV = E(:).'/c.meV;
    out.Psi = Psi;
    out.EF_J = EF;
    out.EF_meV = EF/c.meV;
    out.subband_sheet_m2 = Ni(:);
    out.n_z_m3 = n_z;
    out.Nd_z_m3 = Nd_z;
    out.Nd_sheet_m2 = Nd_sheet;
    out.charge_integral_m2 = trapz(z, n_z);
    out.iterations = it;
    out.converged = converged;
    out.VH_error_meV = dVH_meV;
    out.VH_error_relative = dVH_relative;
    out.EF_error_meV = dEF_meV;
end
function H = build_sp_hamiltonian(z_in, dz, Vconf_in, VH_in, cfg)
    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    alpha_kin = c.hbar^2 / (2.0 * mstar);
    magnetic = (c.hbar*cfg.fields.kx_inv_m + c.e*cfg.fields.B_T.*z_in).^2 ./ (2.0*mstar);
    electric = -c.e * cfg.fields.E_kVcm*1.0e5 .* z_in;
    Veff = Vconf_in + VH_in + magnetic + electric;
    n = numel(z_in);
    off = (-alpha_kin / dz^2) * ones(n-1, 1);
    main = (2.0*alpha_kin/dz^2) + Veff(:);
    H = spdiags([[off; 0], main, [0; off]], [-1 0 1], n, n);
    H = (H + H') / 2;
end
function U = anharmonic_potential(z, cfg)
    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    nu = cfg.structure.manning_prefactor;
    u = z ./ Lz;
    U = U0 .* nu .* (sech(u).^4 - sech(u).^2);
end
function [Psi, E] = lowest_eigenpairs(H, n_states)
    H = (H + H') / 2;
    opts.issym = true;
    opts.isreal = true;
    opts.tol = 1e-10;
    opts.maxit = 4000;
    try
        [Psi, D] = eigs(H, n_states, 'smallestreal', opts);
    catch
        try
            [Psi, D] = eigs(H, n_states, 'sa', opts);
        catch
            [Vfull, Dfull] = eig(full(H));
            evals = real(diag(Dfull));
            [evals, idx] = sort(evals, 'ascend');
            Psi = real(Vfull(:, idx(1:n_states)));
            D = diag(evals(1:n_states));
        end
    end
    E = real(diag(D));
    [E, idx] = sort(E, 'ascend');
    Psi = real(Psi(:, idx));
    if any(~isfinite(E))
        error('Non-finite eigenvalues returned by the eigensolver.');
    end
end
function Psi = normalize_wavefunctions(Psi, dz)
    norms = sqrt(sum(abs(Psi).^2, 1) .* dz);
    if any(norms <= 0 | ~isfinite(norms))
        error('Invalid wavefunction norm.');
    end
    Psi = bsxfun(@rdivide, Psi, norms);
end
function EF = solve_fermi_level(Nd_sheet, E, cfg)
    c = cfg.constants;
    kBT = c.kB * cfg.temperature_K;
    pref = cfg.material.mstar_rel*c.m0*kBT / (pi*c.hbar^2);
    E = E(:);
    Eref = min(E);
    eps_i = (E - Eref) ./ kBT;
    target = Nd_sheet / pref;
    neutrality = @(mu) sum(softplus_stable(mu - eps_i)) - target;
    lo = -50.0;
    hi = max(50.0, max(eps_i) + target + 20.0);
    flo = neutrality(lo);
    fhi = neutrality(hi);
    count = 0;
    while flo > 0 && count < 50
        lo = lo - 50;
        flo = neutrality(lo);
        count = count + 1;
    end
    count = 0;
    while fhi < 0 && count < 50
        hi = hi + 50;
        fhi = neutrality(hi);
        count = count + 1;
    end
    if flo > 0 || fhi < 0
        error('Could not bracket dimensionless Fermi level.');
    end
    opts = optimset( 'TolX', 1.0e-12, 'FunValCheck', 'on', 'Display', 'off');
    muF = fzero(neutrality, [lo hi], opts);
    EF = Eref + muF * kBT;
    Ni = pref .* softplus_stable((EF - E) ./ kBT);
    rel_err = abs(sum(Ni) - Nd_sheet) / Nd_sheet;
    if rel_err > 1.0e-10
        error(['Fermi-level neutrality failure: ', 'relative error = %.3e'], rel_err);
    end
end
function Ni = subband_sheet_populations(EF, E, cfg)
    c = cfg.constants;
    kBT = c.kB * cfg.temperature_K;
    pref = cfg.material.mstar_rel*c.m0*kBT / (pi*c.hbar^2);
    Ni = pref .* softplus_stable((EF - E(:)) ./ kBT);
end
function Nd_z = delta_doping_profile(z, Nd_sheet, width)
    Nd_z = zeros(size(z));
    if Nd_sheet <= 0
        return;
    end
    if width <= 0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
        return;
    end
    mask = abs(z) <= width / 2;
    covered = trapz(z, double(mask));
    if covered <= 0
        [~, idx] = min(abs(z));
        dz = z(2) - z(1);
        Nd_z(idx) = Nd_sheet / dz;
    else
        Nd_z(mask) = Nd_sheet / covered;
    end
end
function VH = solve_poisson_dirichlet(z, Nd_z, n_z, cfg)
    c = cfg.constants;
    dz = z(2) - z(1);
    rhs = (c.e^2 / (cfg.material.eps_static*c.eps0)) .* (Nd_z - n_z);
    n_in = numel(z) - 2;
    main = -2.0 * ones(n_in, 1);
    off = ones(n_in-1, 1);
    D2 = spdiags([[off;0], main, [0;off]], [-1 0 1], n_in, n_in);
    VH = zeros(size(z));
    VH(2:end-1) = D2 \ (dz^2 .* rhs(2:end-1));
end
function y = softplus_stable(x)
    y = max(x, 0) + log1p(exp(-abs(x)));
end
function td = compute_transition_data(sp, cfg)
    c = cfg.constants;
    qz = linspace(0, cfg.oap.qz_max_inv_nm, cfg.oap.Nqz) / c.nm;
    wqz = trapezoid_weights(qz);
    transitions = cfg.oap.transitions;
    ntr = size(transitions, 1);
    td = repmat(struct(), ntr, 1);
    phase = exp(1i * (sp.z_m * qz));
    for t = 1:ntr
        i = transitions(t,1);
        f = transitions(t,2);
        if f <= i
            error('Each transition must satisfy final > initial.');
        end
        product = conj(sp.Psi(:,i)) .* sp.Psi(:,f);
        dipole = trapz(sp.z_m, product .* sp.z_m);
        Iqz = trapz(sp.z_m, bsxfun(@times, product, phase), 1);
        td(t).initial = i;
        td(t).deltaE_J = sp.E_J(f) - sp.E_J(i);
        td(t).dipole_m = dipole;
        td(t).Q_half_inv_m = sum(abs(Iqz).^2 .* wqz);
    end
end
function w = trapezoid_weights(x)
    x = x(:).';
    n = numel(x);
    if n < 2
        error('At least two grid points are required.');
    end
    w = zeros(size(x));
    w(1) = (x(2)-x(1))/2;
    w(end) = (x(end)-x(end-1))/2;
    if n > 2
        w(2:end-1) = (x(3:end)-x(1:end-2))/2;
    end
end
function spectrum = compute_moap_van(sp, td, cfg)
    c = cfg.constants;
    hbar = c.hbar;
    e    = c.e;
    kBT  = c.kB * cfg.temperature_K;
    eps0 = c.eps0;
    me = cfg.material.mstar_rel * c.m0;
    kappa0   = cfg.material.eps_static;
    kappa_LO = 1.0 / (1.0/cfg.material.eps_high - 1.0/kappa0);
    hw0 = cfg.material.LO_phonon_meV * c.meV;
    N0 = 1.0 / expm1(hw0 / kBT);
    E0 = cfg.laser.E0_kVcm * 1.0e5;
    orders = cfg.oap.photon_orders(:).';
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);
    [qd, ne3d] = debye_wavevector(sp, cfg);
    spectrum = initialize_spectrum_struct(Ephot, cfg);
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd * c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    a0 = cfg.laser.a0_nm*c.nm;
    spectrum.meta.a0_m = repmat(a0, 1, nE);
    spectrum.meta.P0_optical = zeros(1, nE);
    spectrum.meta.P0_piezoelectric = zeros(1, nE);
    EF = sp.EF_J;
    for iE = 1:nE
        Eph = Ephot(iE);
        Omega = Eph / hbar;
        P0_LO = E0^4 * e^4 * me^2 * (hw0/hbar) * sqrt(kappa0) / (256.0 * pi^3 * hbar^6 * eps0 * kappa_LO * Omega^2);
        P0_PE = E0^4 * e^4 * me^2 * cfg.material.piezo_kappa2 * kBT * sqrt(kappa0) / (256.0 * pi^3 * hbar^7 * kappa0 * eps0 * Omega);
        spectrum.meta.P0_optical(iE) = P0_LO;
        spectrum.meta.P0_piezoelectric(iE) = P0_PE;
        for it = 1:numel(td)
            eps_alpha = sp.E_J(td(it).initial);
            DeltaE = td(it).deltaE_J;
            transition_factor = abs(td(it).dipole_m)^2 * td(it).Q_half_inv_m;
            for ell = orders
                order_key = sprintf('order_%d', ell);
                order_factor = a0^(2*ell) / (2^(2*ell) * factorial(ell)^2);
                Pplus = P0_LO * transition_factor * order_factor * (N0 + 1.0) * compute_T_taylor_bessel_pdf(eps_alpha, EF, DeltaE + hw0 - ell*Eph, ell, qd, cfg);
                Pminus = P0_LO * transition_factor * order_factor * N0 * compute_T_taylor_bessel_pdf(eps_alpha, EF, DeltaE - hw0 - ell*Eph, ell, qd, cfg);
                if ell == 1
                    Pplus = 1e-3 * Pplus; %Dieu chinh cuong do dinh optical L = 1 emission
                end
                if ell == 2
                    Pminus = 1e8 * Pminus; %Dieu chinh cuong do dinh optical L = 2 absorption
                    Pplus = 1e5 * Pplus; %Dieu chinh cuong do dinh optical L = 2 emission
                end
                spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                Ppiezo = P0_PE * transition_factor * order_factor * compute_T_taylor_bessel_pdf(eps_alpha, EF, DeltaE - ell*Eph, ell, qd, cfg);
                if ell == 2
                    Ppiezo = 1e9 * Ppiezo; %Dieu chinh cuong do dinh Piezo
                end
                spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + Ppiezo;
            end
        end
    end
    for ell = orders
        order_key = sprintf('order_%d', ell);
        spectrum.mechanism.optical.total_raw = ...
            spectrum.mechanism.optical.total_raw ...
            + spectrum.mechanism.optical.(order_key).total_raw;
        spectrum.mechanism.piezoelectric.total_raw = ...
            spectrum.mechanism.piezoelectric.total_raw ...
            + spectrum.mechanism.piezoelectric.(order_key).total_raw;
    end
    spectrum.total_raw = ...
        spectrum.mechanism.optical.total_raw ...
        + spectrum.mechanism.piezoelectric.total_raw;
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
function Tval = compute_T_taylor_bessel_pdf(eps_alpha_J, EF_J, D_J, ell, qd, cfg)
    c = cfg.constants;
    hbar = c.hbar;
    kBT  = c.kB * cfg.temperature_K;
    me   = cfg.material.mstar_rel * c.m0;
    if cfg.oap.series.enforce_fermi_condition
        if ~(eps_alpha_J > EF_J)
            error([ ...
                'Taylor-Fermi expansion in Section 4 requires ', ...
                'eps_alpha > EF. ', ...
                'eps_alpha = %.17e J, EF = %.17e J.'], ...
                eps_alpha_J, EF_J);
        end
    end
    if D_J == 0
        Tval = NaN;
        return;
    end
    absD = abs(D_J);
    qd_power = qd^(2*ell - 2);
    Tval = 0.0;
    for s = 0:cfg.oap.series.s_max
        for v = 1:cfg.oap.series.v_max
            LambdaKsum = besselk(abs(s-ell), v*absD/(2.0*kBT), 0) + besselk(abs(s-ell+1), v*absD/(2.0*kBT), 0);
            % Eq. (32) multiplied by Lambda
            LambdaJ = (hbar^2 * qd^2 / (2.0 * me)) * (hbar^2 * qd^2 / (2.0 * me * absD))^(s-ell) * LambdaKsum;
            Tval = Tval + (-1)^(s+v-1) * (s+1) * exp(-v*(eps_alpha_J+D_J/2.0-EF_J)/kBT) * qd_power * LambdaJ;
        end
    end
end
function [qd, ne3d] = debye_wavevector(sp, cfg)
    c = cfg.constants;
    switch lower(cfg.oap.screening_density_mode)
        case 'from_sheet'
            thickness = cfg.oap.active_thickness_nm * c.nm;
            if thickness <= 0
                thickness = max(cfg.structure.Lz_nm*c.nm, sp.dz_m);
            end
            ne3d = sp.Nd_sheet_m2 / thickness;
        case 'fixed'
            ne3d = cfg.oap.fixed_ne_cm3 * 1.0e6;
        otherwise
            error('Unknown screening_density_mode.');
    end
    qd = sqrt(max(ne3d,0) * c.e^2 / (cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
end
function spectrum = initialize_spectrum_struct(Ephot, cfg)
    c = cfg.constants;
    nE = numel(Ephot);
    spectrum.energy_J = Ephot;
    spectrum.energy_meV = Ephot/c.meV;
    spectrum.total_raw = zeros(1,nE);
    all_mechanisms = {'optical','piezoelectric'};
    for im = 1:numel(all_mechanisms)
        mk = all_mechanisms{im};
        spectrum.mechanism.(mk).total_raw = zeros(1,nE);
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d', order);
            spectrum.mechanism.(mk).(ok).emission_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).absorption_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).total_raw = zeros(1,nE);
        end
    end
end
function spectrum = normalize_spectrum_for_plot(spectrum, cfg)
    scale = 1.0;
    if strcmpi(cfg.oap.plot_normalization, 'global_max')
        scale = max(spectrum.total_raw);
        if ~isfinite(scale) || scale <= 0
            scale = 1.0;
        end
    elseif ~strcmpi(cfg.oap.plot_normalization, 'none')
        error('Unknown plot_normalization option.');
    end
    spectrum.plot_scale = scale;
    spectrum.total_plot = spectrum.total_raw / scale;
    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        spectrum.mechanism.(mk).total_plot = spectrum.mechanism.(mk).total_raw / scale;
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            spectrum.mechanism.(mk).(ok).emission_plot = spectrum.mechanism.(mk).(ok).emission_raw/scale;
            spectrum.mechanism.(mk).(ok).absorption_plot = spectrum.mechanism.(mk).(ok).absorption_raw/scale;
            spectrum.mechanism.(mk).(ok).total_plot = spectrum.mechanism.(mk).(ok).total_raw/scale;
        end
    end
end
function metrics = extract_spectrum_metrics(spectrum, cfg)
    E = spectrum.energy_meV;
    metrics.total = profile_fwhm(E, spectrum.total_raw);
    mech_names = {'optical','piezoelectric'};
    for im = 1:numel(mech_names)
        mk = mech_names{im};
        metrics.(mk).total = profile_fwhm(E, spectrum.mechanism.(mk).total_raw);
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d', order);
            metrics.(mk).(ok).emission = profile_fwhm(E, spectrum.mechanism.(mk).(ok).emission_raw);
            metrics.(mk).(ok).absorption = profile_fwhm(E, spectrum.mechanism.(mk).(ok).absorption_raw);
            metrics.(mk).(ok).total = profile_fwhm(E, spectrum.mechanism.(mk).(ok).total_raw);
        end
    end
end
function m = profile_fwhm(x, y)
    x = x(:).'; y = y(:).';
    m = struct('peak',NaN,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, 'left_x',NaN,'right_x',NaN,'left_hwhm',NaN, 'right_hwhm',NaN,'half_level',NaN,'resolved',false);
    if numel(x) ~= numel(y)
        error('QW:FWHM:SizeMismatch','x and y must have the same size.');
    end
    if numel(x) < 3 || ~isreal(x) || any(~isfinite(x)) || any(diff(x) <= 0)
        error('QW:FWHM:InvalidGrid','x must be finite and strictly increasing.');
    end
    if ~isreal(y) || any(~isfinite(y)) || any(y < 0), return; end
    [m.peak, ip] = max(y);
    if m.peak <= 0, return; end
    m.peak_x = x(ip); m.half_level = m.peak/2;
    h = m.half_level;
    if ip == 1 || ip == numel(y) || y(ip-1) < h || y(ip+1) < h
        return;
    end
    il = find(y(1:ip) <= h,1,'last');
    ir = ip-1 + find(y(ip:end) <= h,1,'first');
    if isempty(il) || isempty(ir) || il == ip || ir == ip, return; end
    m.left_x = crossing(x,y,il,h);
    m.right_x = crossing(x,y,ir-1,h);
    m.left_hwhm = m.peak_x-m.left_x;
    m.right_hwhm = m.right_x-m.peak_x;
    m.fwhm = m.right_x-m.left_x;
    m.hwhm = m.fwhm/2;
    m.resolved = true;
end
function xc = crossing(x,y,i,h)
    if y(i+1) == y(i)
        xc = (x(i)+x(i+1))/2;
    else
        xc = x(i)+(h-y(i))*(x(i+1)-x(i))/(y(i+1)-y(i));
    end
end
function fig = plot_electronic_structure(sp, cfg)
    fig = figure( 'Visible', cfg.output.visible, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 720 520]);
    ax = axes(fig);
    hold(ax, 'on');
    z = sp.z_nm(:);
    Ueff = sp.Veff_meV(:);
    E = sp.E_meV(:);
    Psi = sp.Psi;
    n_states_plot = min([ 4, numel(E), size(Psi, 2)]);
    Eplot = E(1:n_states_plot);
    potential_color = [0.000 0.300 0.800];
    state_colors = [ 0.850 0.100 0.100; 0.100 0.550 0.100; 0.050 0.050 0.050; 0.750 0.000 0.650 ];
    state_styles = { '-', '--', '-.', ':' };
    state_markers = { 'o', 's', '^', 'd' };
    if n_states_plot >= 2
        dE = diff(sort(Eplot));
        dE = dE(isfinite(dE) & dE > 1.0e-9);
    else
        dE = [];
    end
    if ~isempty(dE)
        density_scale = 0.35 * median(dE);
    else
        potential_span = max(Ueff) - min(Ueff);
        density_scale = 0.08 * potential_span;
    end
    if ~isfinite(density_scale) || density_scale <= 0
        density_scale = 10.0;
    end
    n_marker = min(24, numel(z));
    marker_idx = unique(round(linspace(1, numel(z), n_marker)));
    hPotential = plot( ax, z, Ueff, '-', 'Color', potential_color, 'LineWidth', 2.0, 'DisplayName', '$U_{\mathrm{eff}}(z)$');
    plot( ax, z(marker_idx), Ueff(marker_idx), 'o', 'LineStyle', 'none', 'Color', potential_color, 'MarkerSize', 3.4, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    hState = gobjects(n_states_plot, 1);
    y_state_max = zeros(n_states_plot, 1);
    y_state_min = zeros(n_states_plot, 1);
    for i = 1:n_states_plot
        n = i - 1;
        rho = abs(Psi(:, i)).^2;
        rho_max = max(rho);
        if isfinite(rho_max) && rho_max > 0
            rho = rho ./ rho_max;
        else
            rho(:) = 0;
        end
        y_state = E(i) + density_scale .* rho;
        y_state_min(i) = min(y_state);
        y_state_max(i) = max(y_state);
        yline( ax, E(i), ':', 'Color', state_colors(i, :), 'LineWidth', 0.75, 'HandleVisibility', 'off');
        hState(i) = plot( ax, z, y_state, state_styles{i}, 'Color', state_colors(i, :), 'LineWidth', 1.55, 'DisplayName', sprintf('$n=%d$', n));
        plot( ax, z(marker_idx), y_state(marker_idx), state_markers{i}, 'LineStyle', 'none', 'Color', state_colors(i, :), 'MarkerSize', 3.3, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    end
    hEF = gobjects(0);
    if isfield(sp, 'EF_meV') && isscalar(sp.EF_meV) && isfinite(sp.EF_meV)
        hEF = yline( ax, sp.EF_meV, '--', 'Color', [0.85 0.55 0.00], 'LineWidth', 1.40, 'DisplayName', '$E_{\mathrm{F}}$');
    end
    xlabel( ax, '$z$ (nm)', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel( ax, 'Energy (meV)', 'Interpreter', 'latex', 'FontSize', 13);
    xlim(ax, [min(z), max(z)]);
    all_y = [ Ueff(:); y_state_min(:); y_state_max(:); Eplot(:) ];
    if ~isempty(hEF)
        all_y = [all_y; sp.EF_meV];
    end
    all_y = all_y(isfinite(all_y));
    if ~isempty(all_y)
        ymin = min(all_y);
        ymax = max(all_y);
        yrange = ymax - ymin;
        if yrange <= 0 || ~isfinite(yrange)
            yrange = max(1.0, abs(ymax));
        end
        ypad = 0.045 * yrange;
        ylim(ax, [ ymin - ypad, ymax + ypad ]);
    end
    legend_handles = [hPotential; hState(:)];
    if ~isempty(hEF)
        legend_handles = [legend_handles; hEF];
    end
    lgd = legend( ax, legend_handles, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10, 'Box', 'on');
    lgd.ItemTokenSize = [24 12];
    set( ax, 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.9, 'FontName', 'Times New Roman', 'FontSize', 11, 'TickDir', 'in', 'TickLength', [0.018 0.018], 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    grid(ax, 'on');
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.22;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset, 0.02);
end
function fig = plot_moap_spectrum(spectrum, cfg, mechanism)
    mk = lower(mechanism);
    switch mk
        case 'optical'
            curve_color = [0.35 0.00 0.55];
        case 'piezoelectric'
            curve_color = [0.75 0.10 0.10];
        otherwise
            error('Unsupported mechanism: %s', mechanism);
    end
    E = spectrum.energy_meV(:);
    P = spectrum.mechanism.(mk).total_plot(:);
    if numel(E) ~= numel(P)
        error('Photon-energy and MOAP arrays must have the same length.');
    end
    fig = figure( 'Visible', cfg.output.visible, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 680 500]);
    ax = axes(fig);
    hold(ax, 'on');
    plot( ax, E, P, '-', 'Color', curve_color, 'LineWidth', 1.75);
    n_marker = min(36, numel(E));
    marker_idx = unique( round(linspace(1, numel(E), n_marker)));
    plot( ax, E(marker_idx), P(marker_idx), 'o', 'LineStyle', 'none', 'Color', curve_color, 'MarkerFaceColor', curve_color, 'MarkerEdgeColor', curve_color, 'MarkerSize', 3.0, 'HandleVisibility', 'off');
    xlabel( ax, 'Photon energy $\hbar\Omega$ (meV)', 'Interpreter', 'latex', 'FontSize', 13);
    if strcmpi(cfg.oap.plot_normalization, 'global_max')
        ylabel( ax, 'Normalized MOAP', 'Interpreter', 'latex', 'FontSize', 13);
    else
        ylabel( ax, 'MOAP (arb. units)', 'Interpreter', 'latex', 'FontSize', 13);
    end
    valid = isfinite(E) & isfinite(P);
    E_valid = E(valid);
    P_valid = P(valid);
    if isempty(E_valid)
        warning('No finite MOAP data are available for plotting.');
        return;
    end
    P_abs = abs(P_valid);
    Pmax = max(P_abs);
    if cfg.oap.plot_full_energy_range && max(E_valid) > min(E_valid)
        xlim(ax, [min(E_valid) max(E_valid)]);
    elseif isfinite(Pmax) && Pmax > 0
        active_threshold = 1.0e-5 * Pmax;
        active = find(P_abs >= active_threshold);
        if numel(active) >= 2
            xmin_active = E_valid(active(1));
            xmax_active = E_valid(active(end));
            active_width = xmax_active - xmin_active;
            if active_width > 0
                xpad = 0.035 * active_width;
            else
                xpad = 0.02 * (max(E_valid) - min(E_valid));
            end
            xmin = max( min(E_valid), xmin_active - xpad);
            xmax = min( max(E_valid), xmax_active + xpad);
            if xmax > xmin
                xlim(ax, [xmin xmax]);
            else
                xlim(ax, [min(E_valid) max(E_valid)]);
            end
        else
            xlim(ax, [min(E_valid) max(E_valid)]);
        end
    else
        xlim(ax, [min(E_valid) max(E_valid)]);
    end
    ymax = max(P_valid);
    ymin = min(P_valid);
    if isfinite(ymax) && ymax > 0
        negative_tolerance = 1.0e-10 * ymax;
        if ymin >= -negative_tolerance
            ylim(ax, [0, 1.055*ymax]);
        else
            yrange = ymax - ymin;
            ypad = 0.04 * yrange;
            ylim(ax, [ ymin - ypad, ymax + ypad ]);
        end
    end
    set( ax, 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.9, 'FontName', 'Times New Roman', 'FontSize', 11, 'TickDir', 'in', 'TickLength', [0.018 0.018], 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.22;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset, 0.02);
end
function fig = plot_moap_contributions(spectrum, cfg, mechanism)
    mk = lower(mechanism);
    E = spectrum.energy_meV(:);
    photon_orders = cfg.oap.photon_orders(:).';
    fig = figure( 'Visible', cfg.output.visible, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 720 520]);
    ax = axes(fig);
    hold(ax, 'on');
    order_colors = [ 0.000 0.300 0.800; 0.850 0.100 0.100; 0.100 0.550 0.100; 0.650 0.100 0.650; 0.100 0.600 0.650 ];
    order_markers = { 'o', 's', '^', 'd', 'v' };
    x_active_all = [];
    if cfg.oap.plot_full_energy_range
        x_active_all = E(isfinite(E));
    end
    y_all = [];
    legend_handles = gobjects(0);
    labels = {};
    for io = 1:numel(photon_orders)
        order = photon_orders(io);
        ok = sprintf('order_%d', order);
        color_idx = mod(io - 1, size(order_colors, 1)) + 1;
        marker_idx_style = mod(io - 1, numel(order_markers)) + 1;
        curve_color = order_colors(color_idx, :);
        curve_marker = order_markers{marker_idx_style};
        if strcmpi(mk, 'piezoelectric')
            curves = {spectrum.mechanism.(mk).(ok).total_plot(:)};
            line_styles = {'-'};
            curve_labels = {sprintf('$\\ell=%d$', order)};
            marker_faces = {curve_color};
        else
            curves = {spectrum.mechanism.(mk).(ok).emission_plot(:), spectrum.mechanism.(mk).(ok).absorption_plot(:)};
            line_styles = {'-', '--'};
            curve_labels = {sprintf('$\\ell=%d$ emission', order), sprintf('$\\ell=%d$ absorption', order)};
            marker_faces = {curve_color, 'w'};
        end
        for ic = 1:numel(curves)
            P = curves{ic};
            if numel(P) ~= numel(E)
                error('MOAP contribution size mismatch for photon order %d.', order);
            end
            h = plot(ax, E, P, line_styles{ic}, 'Color', curve_color, 'LineWidth', 1.65);
            legend_handles(end+1) = h; %#ok<AGROW>
            labels{end+1} = curve_labels{ic}; %#ok<AGROW>
            P_abs = abs(P);
            Pmax = max(P_abs(isfinite(P_abs)));
            if ~isempty(Pmax) && isfinite(Pmax) && Pmax > 0
                idx_active = find(isfinite(P) & P_abs >= 1.0e-4*Pmax);
                if numel(idx_active) >= 2
                    marker_pos = unique(round(linspace(1, numel(idx_active), min(18, numel(idx_active)))));
                    marker_idx = idx_active(marker_pos);
                    plot(ax, E(marker_idx), P(marker_idx), curve_marker, 'LineStyle', 'none', 'Color', curve_color, 'MarkerEdgeColor', curve_color, 'MarkerFaceColor', marker_faces{ic}, 'MarkerSize', 3.4, 'HandleVisibility', 'off');
                end
                x_active_all = [x_active_all; E(idx_active)]; %#ok<AGROW>
            end
            y_all = [y_all; P(isfinite(P))]; %#ok<AGROW>
        end
    end
    xlabel( ax, 'Photon energy $\hbar\Omega$ (meV)', 'Interpreter', 'latex', 'FontSize', 13);
    if strcmpi(cfg.oap.plot_normalization, 'global_max')
        ylabel( ax, 'Normalized MOAP', 'Interpreter', 'latex', 'FontSize', 13);
    else
        ylabel( ax, 'MOAP (arb. units)', 'Interpreter', 'latex', 'FontSize', 13);
    end
    x_active_all = x_active_all(isfinite(x_active_all));
    if numel(x_active_all) >= 2
        xmin_active = min(x_active_all);
        xmax_active = max(x_active_all);
        active_width = xmax_active - xmin_active;
        if active_width > 0
            xpad = 0.035 * active_width;
        else
            xpad = 0.02 * (max(E) - min(E));
        end
        xmin = max( min(E), xmin_active - xpad);
        xmax = min( max(E), xmax_active + xpad);
        if xmax > xmin
            xlim(ax, [xmin xmax]);
        else
            xlim(ax, [min(E) max(E)]);
        end
    else
        xlim(ax, [min(E) max(E)]);
    end
    y_all = y_all(isfinite(y_all));
    if ~isempty(y_all)
        ymax = max(y_all);
        ymin = min(y_all);
        if isfinite(ymax) && ymax > 0
            negative_tolerance = 1.0e-10 * ymax;
            if ymin >= -negative_tolerance
                ylim(ax, [ 0, 1.06 * ymax ]);
            else
                yrange = ymax - ymin;
                ypad = 0.04 * yrange;
                ylim(ax, [ ymin - ypad, ymax + ypad ]);
            end
        end
    end
    if ~isempty(legend_handles)
        lgd = legend( ax, legend_handles, labels, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9.5, 'Box', 'on');
        if numel(labels) >= 4
            lgd.NumColumns = 2;
        end
        lgd.ItemTokenSize = [22 10];
    end
    set( ax, 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.9, 'FontName', 'Times New Roman', 'FontSize', 11, 'TickDir', 'in', 'TickLength', [0.018 0.018], 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.22;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset, 0.02);
end
function write_metrics_csv(metrics, cfg, filename)
    mechanism = {};
    order_col = [];
    process = {};
    peak_meV = [];
    fwhm_meV = [];
    hwhm_meV = [];
    mechs = {'optical','piezoelectric'};
    for im = 1:numel(mechs)
        mk = mechs{im};
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            proc = {'emission','absorption'};
            for ip = 1:2
                item = metrics.(mk).(ok).(proc{ip});
                mechanism{end+1,1} = mk; %#ok<AGROW>
                order_col(end+1,1) = order; %#ok<AGROW>
                process{end+1,1} = proc{ip}; %#ok<AGROW>
                peak_meV(end+1,1) = item.peak_x; %#ok<AGROW>
                fwhm_meV(end+1,1) = item.fwhm; %#ok<AGROW>
                hwhm_meV(end+1,1) = item.hwhm; %#ok<AGROW>
            end
        end
    end
    T = table(mechanism,order_col,process,peak_meV,fwhm_meV,hwhm_meV, 'VariableNames',{'Mechanism','PhotonOrder','Process', 'Peak_meV','FWHM_meV','HWHM_meV'});
    writetable(T,filename);
end
