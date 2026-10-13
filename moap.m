function results = moap(cfg)
    close all;
    project_root = fileparts(mfilename('fullpath'));
    if nargin == 0 || isempty(cfg)
        cfg = default_config();
        cfg.run.task = 'single';
        cfg.structure.Lz_nm       = 15.0;
        cfg.structure.U0_meV      = 228.0;
        cfg.doping.Nd_sheet_cm2   = 1.0e11;
        cfg.fields.B_T            = 5.0;
        cfg.fields.E_kVcm         = 25.0;
        cfg.temperature_K         = 50.0;
        cfg.oap.model             = 'direct_approx';
        cfg.oap.transitions       = [1 2];
        cfg.oap.photon_orders     = [1 2];
        cfg.oap.mechanisms        = {'optical', 'piezoelectric'};
        cfg.oap.photon_energy_meV = 1.0:0.01:400.0;
        cfg.output.directory      = fullfile(project_root, 'results_approx');
        cfg.output.save_figures   = true;
        cfg.output.save_csv       = true;
        cfg.output.save_mat       = true;
        cfg.output.visible        = 'on';
    elseif nargin ~= 1 || ~isstruct(cfg)
        error('direct_approx_all_in_one expects zero arguments or one cfg struct.');
    end
 
    validate_config(cfg);
    ensure_directory(cfg.output.directory);
    results = struct();
    results.config = cfg;
    results.single = run_single_case(cfg, 'single');
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'project_results.mat'), 'results', 'cfg', '-v7.3');
    end
    fprintf('\n[MOAP] Completed successfully.\n');
    fprintf('[MOAP] Results folder: %s\n', cfg.output.directory);
end
function cfg = default_config()
%DEFAULT_CONFIG Central configuration for the complete QW-MOAP workflow.
    c = physical_constants();
    cfg.run.task = 'single';
    cfg.run.verbose = true;
    cfg.run.debug_moap_direct = false;
    cfg.material.name = 'GaAs';
    cfg.material.mstar_rel = 0.067;
    cfg.material.eps_static = 13.18;
    cfg.material.eps_high = 10.89;
    cfg.material.LO_phonon_meV = 36.25;
    cfg.material.sound_speed_mps = 5.22e3;
    cfg.material.mass_density_kgm3 = 5317;
    cfg.material.piezo_kappa2 = 0.006;
    cfg.material.e14_Cm2 = 0.16;
    cfg.material.eps_poisson = 13.18;
    cfg.structure.potential_model = 'manning_sech';
    cfg.structure.U0_meV = 220.0;
    cfg.structure.Lz_nm = 12.0;
    cfg.structure.manning_prefactor = 6.0;
    cfg.structure.alpha = 0.30;
    cfg.structure.domain_nm = 60.0;
    cfg.structure.auto_domain = false;
    cfg.structure.domain_factor = 3.0;
    cfg.structure.Nz = 1502;
    cfg.structure.n_states = 6;
    cfg.doping.Nd_sheet_cm2 = 1.0e13;
    cfg.doping.width_nm = 2.0;
    cfg.fields.B_T = 0.0;
    cfg.fields.E_kVcm = 0.0;
    cfg.fields.kx_inv_m = 0.0;
    cfg.temperature_K = 90.0;
    cfg.sp.max_outer = 100;
    cfg.sp.max_inner = 500;
    cfg.sp.mix = 0.40; % stable PDF fixed-point mixing for configured subbands
    cfg.sp.tol_VH_meV = 1.0e-3;
    cfg.sp.tol_VH_relative = 1.0e-3; % PDF relative U_H iteration criterion
    cfg.sp.tol_EF_meV = 1.0e-3;
    cfg.sp.poisson_boundary = 'dirichlet_zero';
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_mode = 'constant'; %'dynamic'
    cfg.laser.a0_nm = 7.5;
    cfg.oap.model = 'direct_q_integral';
    cfg.oap.transitions = [1 2];
    cfg.oap.photon_orders = [1 2 3];
    cfg.oap.mechanisms = {'optical', 'piezoelectric'};
    cfg.oap.photon_energy_meV = 7.5:0.01:100.0;
    cfg.oap.population_mode = 'maxwell_boltzmann';
    cfg.oap.include_pauli_blocking = false;
    cfg.oap.qz_max_inv_nm = 50;
    cfg.oap.qperp_max_inv_nm = 50;
    cfg.oap.Nqz = 100;
    cfg.oap.Nqperp = 100;
    cfg.oap.inplane_factor = 'none';
    cfg.oap.landau_initial = 0;
    cfg.oap.landau_final = 0;
    cfg.oap.thermal_cutoff_factor = 1.0;
    cfg.oap.gamma_optical_meV = 0.80;
    cfg.oap.gamma_piezo_meV = 0.60;
    cfg.oap.screening_density_mode = 'from_sheet';
    cfg.oap.fixed_ne_cm3 = 1.0e18;
    cfg.oap.active_thickness_nm = 10.0;
    cfg.oap.plot_normalization = 'none';
    cfg.oap.plot_full_energy_range = true;
    cfg.oap.series.s_max = 99;   % screening index s = 0:s_max (van), 1:s_max (Nhan)
    cfg.oap.series.eta_max = 5;  % HÆ¯á»šNG DáºªN screening-series index eta
    cfg.oap.series.v_max = 99;   % Refined PDF Fermi index v = 1:v_max
    cfg.oap.series.evaluation = 'auto'; % series = checked Eq. (32), quadrature = Eq. (28)
    cfg.oap.series.enforce_fermi_condition = false;
    cfg.oap.normalization_area_m2 = 1.0;
    cfg.sweep.parameter = 'B_T';
    cfg.sweep.values = [5 10 15];
    cfg.output.directory = fullfile(pwd, 'results');
    cfg.output.save_figures = true;
    cfg.output.save_csv = true;
    cfg.output.save_mat = true;
    cfg.output.visible = 'on';
    cfg.output.figure_format = 'png';
    cfg.output.figure_dpi = 200;
    cfg.constants = c;
end

function c = physical_constants()
%PHYSICAL_CONSTANTS CODATA/SI constants used throughout the project.
    c.e = 1.602176634e-19;
    c.hbar = 1.054571817e-34;
    c.kB = 1.380649e-23;
    c.eps0 = 8.8541878128e-12;
    c.m0 = 9.1093837015e-31;
    c.c0 = 299792458;
    c.meV = c.e * 1.0e-3;
    c.nm = 1.0e-9;
end
function validate_config(cfg)
%VALIDATE_CONFIG Fail early on nonphysical or unsupported inputs.
    assert(cfg.structure.Nz >= 50, 'Nz must be at least 50.');
    assert(cfg.structure.n_states >= 2, 'At least two states are required.');
    assert(cfg.structure.Lz_nm > 0, 'Lz must be positive.');
    assert(cfg.structure.U0_meV > 0, 'U0 must be positive.');
    assert(cfg.temperature_K > 0, 'Temperature must be positive.');
    assert(cfg.doping.Nd_sheet_cm2 >= 0, 'Sheet density cannot be negative.');
    assert(cfg.oap.Nqz >= 8 && cfg.oap.Nqperp >= 8, 'q grids are too small.');
    photon_energy_meV = cfg.oap.photon_energy_meV;
    assert(isnumeric(photon_energy_meV) && isvector(photon_energy_meV) && ~isempty(photon_energy_meV) && isreal(photon_energy_meV) && all(isfinite(photon_energy_meV)), 'Photon-energy grid must be a nonempty finite real vector.');
    assert(all(diff(photon_energy_meV) > 0), 'Photon-energy grid must be strictly increasing.');
    assert(size(cfg.oap.transitions, 2) == 2, 'cfg.oap.transitions must have two columns [initial final].');
    assert(all(cfg.oap.transitions(:) >= 1), 'Transition indices are MATLAB 1-based and must be >= 1.');
    assert(max(cfg.oap.transitions(:)) <= cfg.structure.n_states, 'A transition index exceeds structure.n_states.');
    assert(all(cfg.oap.photon_orders >= 1), 'Photon orders must be positive integers.');
end
function ensure_directory(path_name)
%ENSURE_DIRECTORY Create directory if needed.
    if ~exist(path_name, 'dir')
        mkdir(path_name);
    end
end
function case_result = run_single_case(cfg, tag)
    if nargin < 2
        tag = 'case';
    end
    if cfg.run.verbose
        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('[MOAP] Running case: %s\n', tag);
        fprintf( '[MOAP] B = %.4g T, E = %.4g kV/cm, T = %.4g K\n', cfg.fields.B_T, cfg.fields.E_kVcm, cfg.temperature_K);
        fprintf( '[MOAP] Lz = %.4g nm, alpha = %.4g, Nd = %.4g cm^-2\n', cfg.structure.Lz_nm, cfg.structure.alpha, cfg.doping.Nd_sheet_cm2);
    end
    sp = solve_schrodinger_poisson(cfg);
    if ~sp.converged
        warning( 'Schrodinger-Poisson solver reached max_iter without convergence.');
    end
    transition_data = compute_transition_data(sp, cfg);
    spectrum = compute_moap_direct_approx(sp, transition_data, cfg);
    metrics = extract_spectrum_metrics(spectrum, cfg);
    case_result = struct();
    case_result.cfg = cfg;
    case_result.sp = sp;
    case_result.transitions = transition_data;
    case_result.spectrum = spectrum;
    case_result.metrics = metrics;
    fig_structure = plot_electronic_structure(sp, cfg);
    fig_spectrum_optical = plot_moap_spectrum(spectrum, cfg, tag, 'optical');
    fig_spectrum_piezoelectric = plot_moap_spectrum(spectrum, cfg,tag,'piezoelectric');
    fig_contributions_optical = plot_moap_contributions(spectrum, cfg, tag, 'optical');
    fig_contributions_piezoelectric = plot_moap_contributions(spectrum, cfg, tag, 'piezoelectric');
    figures = { fig_structure, fig_spectrum_optical, fig_spectrum_piezoelectric, fig_contributions_optical, fig_contributions_piezoelectric };
    save_case_outputs(case_result, cfg, tag, figures);
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
%SOLVE_SCHRODINGER_POISSON Self-consistent Eqs. (1)-(6).
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
        tol_VH_relative = 1.0e-3; % PDF self-consistency threshold.
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
        n_z = electron_density_from_subbands(Psi, Ni);
        VH_new = solve_poisson_dirichlet(z, Nd_z, n_z, cfg);
        VH_next = cfg.sp.mix .* VH_new + (1.0-cfg.sp.mix) .* VH_old;
        dVH_meV = max(abs(VH_next - VH_old)) / c.meV;
        VH_scale = max(max(abs(VH_old)), c.meV);
        dVH_relative = max(abs(VH_next - VH_old)) / VH_scale;
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
    n_z = electron_density_from_subbands(Psi, Ni);
    B = cfg.fields.B_T;
    Efield = cfg.fields.E_kVcm*1e5;
    mstar = cfg.material.mstar_rel*c.m0;
    kx = cfg.fields.kx_inv_m;
    Vmag = (c.hbar*kx + c.e*B.*z).^2/(2*mstar);
    Velec = -c.e*Efield.*z;
    Veff = Vconf + VH_old + Vmag + Velec;
    out = struct();
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
%BUILD_SP_HAMILTONIAN Finite-difference Eq. (5)-(6) Hamiltonian.
    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    B = cfg.fields.B_T;
    Efield = cfg.fields.E_kVcm * 1.0e5;
    kx = cfg.fields.kx_inv_m;
    alpha_kin = c.hbar^2 / (2.0 * mstar);
    magnetic = (c.hbar*kx + c.e*B.*z_in).^2 ./ (2.0*mstar);
    electric = -c.e * Efield .* z_in;
    Veff = Vconf_in + VH_in + magnetic + electric;
    n = numel(z_in);
    off = (-alpha_kin / dz^2) * ones(n-1, 1);
    main = (2.0*alpha_kin/dz^2) + Veff(:);
    H = spdiags([[off; 0], main, [0; off]], [-1 0 1], n, n);
    H = (H + H') / 2;
end
function U = anharmonic_potential(z, cfg)
%ANHARMONIC_POTENTIAL Manning confinement in PDF Eq. (1).
    c = cfg.constants;
    Lz = cfg.structure.Lz_nm * c.nm;
    U0 = cfg.structure.U0_meV * c.meV;
    nu = cfg.structure.manning_prefactor;
    u = z ./ Lz;
    U = U0 .* nu .* (sech(u).^4 - sech(u).^2);
end
function [Psi, E] = lowest_eigenpairs(H, n_states)
%LOWEST_EIGENPAIRS Robust lowest eigenpairs for a symmetric sparse matrix.
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
%NORMALIZE_WAVEFUNCTIONS Enforce integral |psi_i|^2 dz = 1.
    norms = sqrt(sum(abs(Psi).^2, 1) .* dz);
    if any(norms <= 0 | ~isfinite(norms))
        error('Invalid wavefunction norm.');
    end
    Psi = bsxfun(@rdivide, Psi, norms);
end
function EF = solve_fermi_level(Nd_sheet, E, cfg)
%SOLVE_FERMI_LEVEL
    c = cfg.constants;
    T = cfg.temperature_K;
    if T <= 0
        error('solve_fermi_level requires T > 0 K.');
    end
    kBT = c.kB * T;
    mstar = cfg.material.mstar_rel * c.m0;
    pref = mstar * kBT / (pi * c.hbar^2);
    E = E(:);
    if Nd_sheet <= 0
        EF = min(E) - 50*kBT;
        return;
    end
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
%SUBBAND_SHEET_POPULATIONS 2D Fermi-Dirac population of each subband.
    c = cfg.constants;
    mstar = cfg.material.mstar_rel * c.m0;
    kBT = c.kB * cfg.temperature_K;
    pref = mstar * kBT / (pi * c.hbar^2);
    Ni = pref .* softplus_stable((EF - E(:)) ./ kBT);
end
function n_z = electron_density_from_subbands(Psi, Ni)
%ELECTRON_DENSITY_FROM_SUBBANDS n(z)=sum_i Ni |psi_i(z)|^2.
    n_z = (abs(Psi).^2) * Ni(:);
end
function Nd_z = delta_doping_profile(z, Nd_sheet, width)
%DELTA_DOPING_PROFILE Finite-width representation of a centered delta layer.
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
%SOLVE_POISSON_DIRICHLET Solve d2 VH/dz2 = e^2/(eps eps0)(Nd-n).
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
%SOFTPLUS_STABLE log(1+exp(x)) without overflow.
    y = max(x, 0) + log1p(exp(-abs(x)));
end
function td = compute_transition_data(sp, cfg)
%COMPUTE_TRANSITION_DATA Dipole and qz form factors for requested transitions.
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
        td(t).final = f;
        td(t).deltaE_J = sp.E_J(f) - sp.E_J(i);
        td(t).deltaE_meV = td(t).deltaE_J / c.meV;
        td(t).dipole_m = dipole;
        td(t).dipole_nm = dipole / c.nm;
        td(t).qz_inv_m = qz;
        td(t).qz_inv_nm = qz * c.nm;
        td(t).Iqz = Iqz;
        td(t).I2 = abs(Iqz).^2;
        td(t).Q_half_inv_m = sum(td(t).I2 .* wqz);
        td(t).Q_full_inv_m = 2.0 * td(t).Q_half_inv_m;
    end
end
function w = trapezoid_weights(x)
%TRAPEZOID_WEIGHTS Weights for a uniform or nonuniform 1D trapezoidal rule.
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
function spectrum = compute_moap_direct_approx(sp, td, cfg)
%COMPUTE_MOAP_DIRECT_APPROX Direct q-space quadrature of Eqs. (19)/(25).
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
    S = cfg.oap.normalization_area_m2;
    orders = cfg.oap.photon_orders(:).';
    orders = orders(orders <= 2); % Match Nhan's 1PA/2PA limit.
    use_optical = any(strcmpi(cfg.oap.mechanisms,'optical'));
    use_piezoelectric = any(strcmpi( cfg.oap.mechanisms,'piezoelectric'));
    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);
    qd = debye_wavevector(sp, cfg);
    qperp = linspace(0,cfg.oap.qperp_max_inv_nm, cfg.oap.Nqperp)/c.nm;
    qz = td(1).qz_inv_m(:);
    I2 = reshape([td.I2],numel(qz),[]);
    a0 = cfg.laser.a0_nm*c.nm;
    spectrum = initialize_spectrum_struct(Ephot, cfg);
    for iE = 1:nE
        Eph = Ephot(iE);
        Omega = Eph / hbar;
        P0_LO = E0^4 * e^4 * me^2 * (hw0/hbar) * S * sqrt(kappa0) / (256.0*pi^3*hbar^6*eps0*kappa_LO*Omega^2);
        P0_PE = E0^4 * e^4 * me^2 * cfg.material.piezo_kappa2 * kBT * S * sqrt(kappa0) / (256.0 * pi^3 * hbar^7 * kappa0 * eps0 * Omega);
        for it = 1:numel(td)
            DeltaE = td(it).deltaE_J;
            transition_total_LO = 0.0;
            transition_total_PE = 0.0;
            for ell = orders
                order_key = sprintf('order_%d', ell);
                order_factor = a0^(2*ell) / (2^(2*ell) * factorial(ell)^2);
                if use_optical
                    Pplus = P0_LO*abs(td(it).dipole_m)^2*order_factor*(N0+1.0)* compute_T_qspace( sp.E_J(td(it).initial),sp.EF_J,DeltaE+hw0-ell*Eph,ell,qd,qperp,qz,I2(:,it),cfg);
                    Pminus = P0_LO*abs(td(it).dipole_m)^2*order_factor*N0* compute_T_qspace( sp.E_J(td(it).initial),sp.EF_J, DeltaE-hw0-ell*Eph,ell,qd, qperp,qz,I2(:,it),cfg);
                    spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                    spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                    spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                    transition_total_LO = transition_total_LO+Pplus+Pminus;
                end
                if use_piezoelectric
                    Ppiezo = P0_PE*abs(td(it).dipole_m)^2* order_factor*compute_T_qspace( sp.E_J(td(it).initial),sp.EF_J, DeltaE-ell*Eph,ell,qd, qperp,qz,I2(:,it),cfg);
                    spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + Ppiezo;
                    transition_total_PE = transition_total_PE+Ppiezo;
                end
            end
            tr_key = sprintf('transition_%d_%d', td(it).initial, td(it).final);
            if ~isfield(spectrum.mechanism.optical.transitions, tr_key)
                spectrum.mechanism.optical.transitions.(tr_key) = zeros(1, nE);
                spectrum.mechanism.piezoelectric.transitions.(tr_key) = zeros(1, nE);
            end
            spectrum.mechanism.optical.transitions.(tr_key)(iE) = transition_total_LO;
            spectrum.mechanism.piezoelectric.transitions.(tr_key)(iE) = transition_total_PE;
        end
    end
    for ell = orders
        order_key = sprintf('order_%d', ell);
        spectrum.mechanism.optical.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.optical.(order_key).total_raw;
        spectrum.mechanism.piezoelectric.total_raw = spectrum.mechanism.piezoelectric.total_raw + spectrum.mechanism.piezoelectric.(order_key).total_raw;
    end
    spectrum.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.piezoelectric.total_raw;
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end
function Tval = compute_T_qspace(eps_alpha,EF,D,ell,qd,qperp,qz, I2_qz,cfg)
% Composite trapz in physical q space: qz first, then qperp.
    c = cfg.constants;
    me = cfg.material.mstar_rel*c.m0;
    E_qperp = c.hbar^2*qperp.^2/(2*me);
    epsilon_star = eps_alpha + D/2 + (me * D^2)./(2*c.hbar^2*qperp.^2);
    occupation = 1./(1+exp((epsilon_star-EF)/(c.kB*cfg.temperature_K)));
    q2 = qperp.^2+qz.^2;
    integrand = q2.^(ell+0.5)./(q2+qd^2).^2 .* abs(E_qperp+D).*occupation.*I2_qz;
    Tval = trapz(qperp,trapz(qz,integrand,1),2);
end
function [qd, ne3d] = debye_wavevector(sp, cfg)
%DEBYE_WAVEVECTOR qd^2 = ne e^2/(eps_s eps0 kBT).
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
function spectrum = normalize_spectrum_for_plot(spectrum, cfg)
%NORMALIZE_SPECTRUM_FOR_PLOT Keep raw results and add plotting arrays.
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
%EXTRACT_SPECTRUM_METRICS Peak position, FWHM and HWHM for each component.
    E = spectrum.energy_meV;
    metrics = struct();
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
%PROFILE_FWHM Width of a resolved peak at half its maximum.
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
%PLOT_ELECTRONIC_STRUCTURE
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
function fig = plot_moap_spectrum(spectrum, cfg, ~, mechanism)
%PLOT_MOAP_SPECTRUM
    if nargin < 4
        error('plot_moap_spectrum requires a mechanism.');
    end
    mk = lower(mechanism);
    if ~isfield(spectrum.mechanism, mk)
        error('Unknown MOAP mechanism: %s', mechanism);
    end
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
    full_energy_range = isfield(cfg.oap,'plot_full_energy_range') && cfg.oap.plot_full_energy_range;
    P_abs = abs(P_valid);
    Pmax = max(P_abs);
    if full_energy_range && max(E_valid) > min(E_valid)
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
    Pfinite = P_valid(isfinite(P_valid));
    if ~isempty(Pfinite)
        ymax = max(Pfinite);
        ymin = min(Pfinite);
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
    end
    set( ax, 'Box', 'on', 'Layer', 'top', 'LineWidth', 0.9, 'FontName', 'Times New Roman', 'FontSize', 11, 'TickDir', 'in', 'TickLength', [0.018 0.018], 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    grid(ax, 'on');
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.22;
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset, 0.02);
end
function fig = plot_moap_contributions(spectrum, cfg, ~, mechanism)
%PLOT_MOAP_CONTRIBUTIONS
    if nargin < 4
        error('plot_moap_contributions requires a mechanism.');
    end
    mk = lower(mechanism);
    if ~isfield(spectrum.mechanism, mk)
        error('Unknown MOAP mechanism: %s', mechanism);
    end
    switch mk
        case 'optical'
        case 'piezoelectric'
        otherwise
            error('Unsupported mechanism: %s', mechanism);
    end
    E = spectrum.energy_meV(:);
    photon_orders = cfg.oap.photon_orders(:).';
    n_orders = numel(photon_orders);
    fig = figure( 'Visible', cfg.output.visible, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 720 520]);
    ax = axes(fig);
    hold(ax, 'on');
    order_colors = [ 0.000 0.300 0.800; 0.850 0.100 0.100; 0.100 0.550 0.100; 0.650 0.100 0.650; 0.100 0.600 0.650 ];
    order_markers = { 'o', 's', '^', 'd', 'v' };
    x_active_all = [];
    if isfield(cfg.oap,'plot_full_energy_range') && cfg.oap.plot_full_energy_range
        x_active_all = E(isfinite(E));
    end
    y_all = [];
    legend_handles = gobjects(0);
    labels = {};
    for io = 1:n_orders
        order = photon_orders(io);
        ok = sprintf('order_%d', order);
        if ~isfield(spectrum.mechanism.(mk), ok)
            continue;
        end
        P_em = spectrum.mechanism.(mk).(ok).emission_plot(:);
        P_ab = spectrum.mechanism.(mk).(ok).absorption_plot(:);
        if numel(P_em) ~= numel(E) || numel(P_ab) ~= numel(E)
            error( 'MOAP contribution size mismatch for photon order %d.', order);
        end
        color_idx = mod(io - 1, size(order_colors, 1)) + 1;
        marker_idx_style = mod(io - 1, numel(order_markers)) + 1;
        curve_color = order_colors(color_idx, :);
        curve_marker = order_markers{marker_idx_style};
        h_em = plot( ax, E, P_em, '-', 'Color', curve_color, 'LineWidth', 1.65);
        legend_handles(end+1) = h_em; %#ok<AGROW>
        labels{end+1} = sprintf( '$%d\\mathrm{PA}$ emission', order); %#ok<AGROW>
        Pem_abs = abs(P_em);
        Pem_max = max(Pem_abs(isfinite(Pem_abs)));
        if ~isempty(Pem_max) && isfinite(Pem_max) && Pem_max > 0
            threshold_em = 1.0e-4 * Pem_max;
            idx_em_active = find( isfinite(P_em) & Pem_abs >= threshold_em);
            if numel(idx_em_active) >= 2
                n_marker_em = min(18, numel(idx_em_active));
                marker_pos = unique(round(linspace( 1, numel(idx_em_active), n_marker_em)));
                marker_idx_em = idx_em_active(marker_pos);
                plot( ax, E(marker_idx_em), P_em(marker_idx_em), curve_marker, 'LineStyle', 'none', 'Color', curve_color, 'MarkerEdgeColor', curve_color, 'MarkerFaceColor', curve_color, 'MarkerSize', 3.4, 'HandleVisibility', 'off');
            end
            x_active_all = [ x_active_all; E(idx_em_active) ]; %#ok<AGROW>
        end
        h_ab = plot( ax, E, P_ab, '--', 'Color', curve_color, 'LineWidth', 1.65);
        legend_handles(end+1) = h_ab; %#ok<AGROW>
        labels{end+1} = sprintf( '$%d\\mathrm{PA}$ absorption', order); %#ok<AGROW>
        Pab_abs = abs(P_ab);
        Pab_max = max(Pab_abs(isfinite(Pab_abs)));
        if ~isempty(Pab_max) && isfinite(Pab_max) && Pab_max > 0
            threshold_ab = 1.0e-4 * Pab_max;
            idx_ab_active = find( isfinite(P_ab) & Pab_abs >= threshold_ab);
            if numel(idx_ab_active) >= 2
                n_marker_ab = min(18, numel(idx_ab_active));
                marker_pos = unique(round(linspace( 1, numel(idx_ab_active), n_marker_ab)));
                marker_idx_ab = idx_ab_active(marker_pos);
                plot( ax, E(marker_idx_ab), P_ab(marker_idx_ab), curve_marker, 'LineStyle', 'none', 'Color', curve_color, 'MarkerEdgeColor', curve_color, 'MarkerFaceColor', 'w', 'MarkerSize', 3.4, 'HandleVisibility', 'off');
            end
            x_active_all = [ x_active_all; E(idx_ab_active) ]; %#ok<AGROW>
        end
        y_all = [ y_all; P_em(isfinite(P_em)); P_ab(isfinite(P_ab)) ]; %#ok<AGROW>
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
function save_case_outputs(case_result, cfg, tag, figures)
%SAVE_CASE_OUTPUTS Save MAT, CSV and figures for one case.
    ensure_directory(cfg.output.directory);
    safe_tag = regexprep( tag, '[^A-Za-z0-9_-]', '_');
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, [safe_tag '.mat']), 'case_result', '-v7.3');
    end
    if cfg.output.save_csv
        E = case_result.spectrum.energy_meV(:);
        total = case_result.spectrum.total_raw(:);
        optical = case_result.spectrum.mechanism.optical.total_raw(:);
        piezo = case_result.spectrum.mechanism.piezoelectric.total_raw(:);
        T = table( E, total, optical, piezo, 'VariableNames', { 'PhotonEnergy_meV', 'TotalRaw', 'OpticalRaw', 'PiezoelectricRaw' });
        writetable( T, fullfile( cfg.output.directory, [safe_tag '_spectrum.csv']));
        write_metrics_csv( case_result.metrics, cfg, fullfile( cfg.output.directory, [safe_tag '_metrics.csv']));
    end
    if cfg.output.save_figures
        names = { 'structure', 'spectrum_optical', 'spectrum_piezoelectric', 'contributions_optical', 'contributions_piezoelectric' };
        nfig = min( numel(figures), numel(names));
        for i = 1:nfig
            if ishghandle(figures{i})
                filename = fullfile( cfg.output.directory, [ safe_tag, '_', names{i}, '.', cfg.output.figure_format ]);
                export_figure_compat( figures{i}, filename, cfg.output.figure_dpi);
            end
        end
    end
end
function export_figure_compat(fig, filename, dpi)
%EXPORT_FIGURE_COMPAT Save figure in old and new MATLAB releases.
    [folder,~,~] = fileparts(filename);
    ensure_directory(folder);
    if exist('exportgraphics','file') == 2
        exportgraphics(fig, filename, 'Resolution', dpi);
    else
        set(fig,'PaperPositionMode','auto');
        print(fig, filename, '-dpng', sprintf('-r%d',dpi));
    end
end
function write_metrics_csv(metrics, cfg, filename)
%WRITE_METRICS_CSV Flatten component peak metrics to CSV.
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
