function results = moap2(cfg)
% MOAP numerical model with local raw-profile FWHM measurement.
    if nargin == 0 || isempty(cfg), cfg = struct(); end
    % Recursive merge permits e.g. moap2(struct('run',struct('task','linewidth_sweep'))).
    cfg = merge_moap2_config(default_moap2_config(),cfg);
    cfg = apply_moap2_defaults(cfg);
    if ~exist(cfg.output.directory, 'dir')
        mkdir(cfg.output.directory);
    end
    task = lower(cfg.run.task);
    results = struct('config',cfg);
    switch task
        case 'single'
            results.single = run_moap2_case(cfg, true, true);
        case 'sweep_spectra'
            results.sweep = run_moap2_sweep(cfg, false, true, false);
        case 'linewidth_sweep'
            results.linewidth = run_moap2_sweep(cfg, true, false, true);
        case 'all_demo'
            results.single = run_moap2_case(cfg, true, true);
            results.sweep = run_moap2_sweep(cfg, true, true, true);
            results.linewidth = struct('parameter',results.sweep.parameter, ...
                'values',results.sweep.values, ...
                'linewidth',results.sweep.linewidth);
        otherwise
            error('Unsupported cfg.run.task: %s', cfg.run.task);
    end
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'project_results.mat'), 'results', 'cfg', '-v7.3');
    end
    fprintf('\n[MOAP] Completed successfully.\n[MOAP] Results folder: %s\n',cfg.output.directory);
end

function cfg = default_moap2_config()
    cfg = struct();
    cfg.run.task = 'all_demo'; % single | sweep_spectra | linewidth_sweep | all_demo
    cfg.run.verbose = true;
    cfg.material.mstar_rel = 0.067;
    cfg.material.eps_static = 13.18;
    cfg.material.eps_high = 10.89;
    cfg.material.LO_phonon_meV = 36.25;
    cfg.material.piezo_kappa2 = 0.006;
    cfg.structure.U0_meV = 228;
    cfg.structure.Lz_nm = 15;
    cfg.structure.manning_prefactor = 6;
    cfg.structure.domain_nm = 60;
    cfg.structure.domain_factor = 3;
    cfg.structure.Nz = 1502;
    cfg.structure.n_states = 6;
    cfg.doping.Nd_sheet_cm2 = 1e11;
    cfg.doping.width_nm = 2;
    cfg.fields.B_T = 5;
    cfg.fields.E_kVcm = 25;
    cfg.fields.kx_inv_m = 0;
    cfg.temperature_K = 100;
    cfg.sp.max_outer = 100;
    cfg.sp.mix = 0.5;
    cfg.sp.tol_VH_relative = 1e-3;
    cfg.sp.tol_EF_meV = 1e-3;
    cfg.laser.E0_kVcm = 4.5;
    cfg.laser.a0_nm = 7.5;
    cfg.oap.transitions = [1 2];
    cfg.oap.photon_orders = [1 2];
    cfg.oap.mechanisms = {'optical','piezoelectric'};
    cfg.oap.photon_energy_meV = 5:0.01:150;
    cfg.oap.qz_max_inv_nm = 2;
    cfg.oap.Nqz = 101;
    cfg.oap.plot_full_energy_range = true;
    cfg.oap.series.s_max = 2;
    cfg.oap.series.v_max = 2;
    cfg.sweep.parameter = 'T_K';
    cfg.sweep.values = 5:50:300;
    cfg.output.directory = fullfile(fileparts(mfilename('fullpath')), 'results_van');
    cfg.output.save_figures = true;
    cfg.output.save_csv = true;
    cfg.output.save_mat = true;
    cfg.output.visible = 'on';
    cfg.output.figure_format = 'png';
    cfg.output.figure_dpi = 200;
    cfg.output.show_experimental_fwhm = true;
    cfg.constants.e = 1.602176634e-19;
    cfg.constants.hbar = 1.054571817e-34;
    cfg.constants.kB = 1.380649e-23;
    cfg.constants.eps0 = 8.8541878128e-12;
    cfg.constants.m0 = 9.1093837015e-31;
    cfg.constants.meV = cfg.constants.e*1e-3;
    cfg.constants.nm = 1e-9;
end

function out = merge_moap2_config(out, overrides)
    assert(isstruct(overrides) && isscalar(overrides), ...
        'Configuration overrides must be a scalar struct.');
    keys = fieldnames(overrides);
    for i = 1:numel(keys)
        k = keys{i};
        if isfield(out,k) && isstruct(out.(k)) && isstruct(overrides.(k))
            out.(k) = merge_moap2_config(out.(k),overrides.(k));
        else
            out.(k) = overrides.(k);
        end
    end
end

function cfg = apply_moap2_defaults(cfg)
    if ~isfield(cfg,'run'), cfg.run = struct(); end
    if ~isfield(cfg.run,'task'), cfg.run.task = 'single'; end
    if ~isfield(cfg.run,'verbose'), cfg.run.verbose = true; end
    if ~isfield(cfg,'oap'), cfg.oap = struct(); end
    if ~isfield(cfg.oap,'mechanisms')
        cfg.oap.mechanisms = {'optical','piezoelectric'};
    end
    if ~isfield(cfg,'sweep'), cfg.sweep = struct(); end
    if ~isfield(cfg.sweep,'parameter'), cfg.sweep.parameter = 'T_K'; end
    if ~isfield(cfg.sweep,'values'), cfg.sweep.values = cfg.temperature_K; end
end

function case_result = run_moap2_case(cfg, make_plots, save_outputs, compute_metrics)
    if nargin < 4,compute_metrics = true; end
    assert(cfg.temperature_K > 0 && isfinite(cfg.temperature_K), ...
        'Temperature must be strictly positive (use e.g. 1 K instead of 0 K).');
    if cfg.run.verbose
        fprintf('\n============================================================\n');
        fprintf('[MOAP] B = %.4g T, E = %.4g kV/cm, T = %.4g K\n', ...
            cfg.fields.B_T,cfg.fields.E_kVcm,cfg.temperature_K);
        fprintf('[MOAP] Lz = %.4g nm, Nd = %.4g cm^-2\n', ...
            cfg.structure.Lz_nm,cfg.doping.Nd_sheet_cm2);
    end
    sp = solve_schrodinger_poisson(cfg);
    if ~sp.converged
        warning('Schrodinger-Poisson solver reached max_iter without convergence.');
    end
    transition_data = compute_transition_data(sp,cfg);
    spectrum = compute_moap(sp,transition_data,cfg);
    if compute_metrics
        metrics = extract_spectrum_metrics(spectrum,cfg);
    else
        metrics = struct();
    end
    case_result = struct('cfg',cfg,'sp',sp,'transitions',transition_data, ...
        'spectrum',spectrum,'metrics',metrics);
    figures = {};
    names = {};
    if make_plots
        figures = {plot_electronic_structure(sp,cfg), ...
            plot_moap_contributions(spectrum,cfg,'optical'), ...
            plot_moap_contributions(spectrum,cfg,'piezoelectric')};
        names = {'structure','contributions_optical','contributions_piezoelectric'};
    end
    if save_outputs
        save_moap2_single(case_result,figures,names,cfg);
    end
    if cfg.run.verbose
        fprintf('[MOAP] Converged: %d in %d iterations\n',sp.converged,sp.iterations);
        fprintf('[MOAP] EF = %.8f meV\n',sp.EF_meV);
        fprintf('[MOAP] Subband energies (meV): '); fprintf('%.8f ',sp.E_meV);
        fprintf('\n[MOAP] Max raw OAP = %.6e\n',max(spectrum.total_raw));
    end
end

function save_moap2_single(case_result, figures, names, cfg)
    spectrum = case_result.spectrum;
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory,'single.mat'),'case_result','-v7.3');
    end
    if cfg.output.save_csv
        T = table(spectrum.energy_meV(:),spectrum.total_raw(:), ...
            spectrum.mechanism.optical.total_raw(:), ...
            spectrum.mechanism.piezoelectric.total_raw(:), ...
            'VariableNames',{'PhotonEnergy_meV','TotalRaw','OpticalRaw','PiezoelectricRaw'});
        writetable(T,fullfile(cfg.output.directory,'single_spectrum.csv'));
        write_fwhm_metrics_csv(case_result.metrics,cfg, ...
            fullfile(cfg.output.directory,'single_metrics.csv'),NaN);
    end
    save_moap2_figures(figures,names,'single',cfg);
end
function out = solve_schrodinger_poisson(cfg)
    c = cfg.constants;
    domain_nm = cfg.structure.domain_nm;
    z = linspace(-domain_nm/2, domain_nm/2, cfg.structure.Nz).' * c.nm;
    dz = z(2) - z(1);
    z_in = z(2:end-1);
    u = z/(cfg.structure.Lz_nm*c.nm);
    Vconf = cfg.structure.U0_meV*c.meV*cfg.structure.manning_prefactor.*(sech(u).^4-sech(u).^2);
    Nd_sheet = cfg.doping.Nd_sheet_cm2 * 1.0e4;
    width = cfg.doping.width_nm * c.nm;
    Nd_z = zeros(size(z));
    if Nd_sheet > 0
        if width <= 0
            [~,idx] = min(abs(z)); Nd_z(idx) = Nd_sheet/dz;
        else
            mask = abs(z) <= width/2; covered = trapz(z,double(mask));
            if covered <= 0
                [~,idx] = min(abs(z)); Nd_z(idx) = Nd_sheet/dz;
            else
                Nd_z(mask) = Nd_sheet/covered;
            end
        end
    end
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
    mstar = cfg.material.mstar_rel*c.m0;
    alpha_kin = c.hbar^2/(2*mstar);
    magnetic = (c.hbar*cfg.fields.kx_inv_m+c.e*cfg.fields.B_T.*z_in).^2/(2*mstar);
    electric = -c.e*cfg.fields.E_kVcm*1e5.*z_in;
    n = numel(z_in);
    offH = (-alpha_kin/dz^2)*ones(n-1,1);
    offD = ones(n-1,1);
    D2 = spdiags([[offD;0],-2*ones(n,1),[0;offD]],[-1 0 1],n,n);
    for it = 1:max_iter
        Veff = Vconf(2:end-1)+VH_old(2:end-1)+magnetic+electric;
        H = spdiags([[offH;0],2*alpha_kin/dz^2+Veff(:),[0;offH]],[-1 0 1],n,n);
        H = (H+H')/2;
        [Psi_in, E] = lowest_eigenpairs(H, n_states);
        norms = sqrt(sum(abs(Psi_in).^2,1)*dz);
        Psi_in = bsxfun(@rdivide,Psi_in,norms);
        Psi = zeros(numel(z), n_states);
        Psi(2:end-1, :) = Psi_in;
        if Nd_sheet > 0
            [EF,Ni] = solve_fermi_level(Nd_sheet,E,cfg);
        else
            EF = NaN;
            Ni = zeros(n_states,1);
            dEF_meV = 0;
        end
        n_z = (abs(Psi).^2) * Ni(:);
        rhs = (c.e^2/(cfg.material.eps_static*c.eps0)).*(Nd_z-n_z);
        VH_new = zeros(size(z));
        VH_new(2:end-1) = D2\(dz^2.*rhs(2:end-1));
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
    Veff = Vconf(2:end-1)+VH_old(2:end-1)+magnetic+electric;
    H = spdiags([[offH;0],2*alpha_kin/dz^2+Veff(:),[0;offH]],[-1 0 1],n,n);
    H = (H+H')/2;
    [Psi_in, E] = lowest_eigenpairs(H, n_states);
    norms = sqrt(sum(abs(Psi_in).^2,1)*dz);
    Psi_in = bsxfun(@rdivide,Psi_in,norms);
    Psi = zeros(numel(z), n_states);
    Psi(2:end-1, :) = Psi_in;
    if Nd_sheet > 0
        [EF,Ni] = solve_fermi_level(Nd_sheet,E,cfg);
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
end
function [EF,Ni] = solve_fermi_level(Nd_sheet, E, cfg)
    c = cfg.constants;
    kBT = c.kB * cfg.temperature_K;
    pref = cfg.material.mstar_rel*c.m0*kBT / (pi*c.hbar^2);
    E = E(:);
    Eref = min(E);
    eps_i = (E - Eref) ./ kBT;
    target = Nd_sheet / pref;
    neutrality = @(mu) sum(max(mu-eps_i,0)+log1p(exp(-abs(mu-eps_i))))-target;
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
    opts = optimset('TolX',1e-12,'Display','off');
    muF = fzero(neutrality, [lo hi], opts);
    EF = Eref + muF * kBT;
    x = (EF-E)/kBT;
    Ni = pref.*(max(x,0)+log1p(exp(-abs(x))));
end
function td = compute_transition_data(sp, cfg)
    c = cfg.constants;
    qz = linspace(0, cfg.oap.qz_max_inv_nm, cfg.oap.Nqz) / c.nm;
    wqz = zeros(size(qz));
    wqz([1 end]) = [qz(2)-qz(1),qz(end)-qz(end-1)]/2;
    wqz(2:end-1) = (qz(3:end)-qz(1:end-2))/2;
    transitions = cfg.oap.transitions;
    ntr = size(transitions, 1);
    td = repmat(struct(), ntr, 1);
    phase = exp(1i * (sp.z_m * qz));
    for t = 1:ntr
        i = transitions(t,1);
        f = transitions(t,2);
        product = conj(sp.Psi(:,i)) .* sp.Psi(:,f);
        dipole = trapz(sp.z_m, product .* sp.z_m);
        Iqz = trapz(sp.z_m, bsxfun(@times, product, phase), 1);
        td(t).initial = i;
        td(t).deltaE_J = sp.E_J(f) - sp.E_J(i);
        td(t).dipole_m = dipole;
        td(t).K = sum(abs(Iqz).^2 .* wqz);
    end
end
function spectrum = compute_moap(sp, td, cfg)
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

    thickness = cfg.structure.Lz_nm*c.nm;
    qd = sqrt(sp.Nd_sheet_m2*c.e^2/(thickness*cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
    

    spectrum.energy_meV = Ephot/c.meV;
    for mechanism = {'optical','piezoelectric'}
        mk = mechanism{1};
        spectrum.mechanism.(mk).total_raw = zeros(1,nE);
        for ell = orders
            ok = sprintf('order_%d',ell);
            spectrum.mechanism.(mk).(ok).emission_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).absorption_raw = zeros(1,nE);
            spectrum.mechanism.(mk).(ok).total_raw = zeros(1,nE);
        end
    end
    a0 = cfg.laser.a0_nm*c.nm;
    for iE = 1:nE
        Eph = Ephot(iE);
        Omega = Eph / hbar;
        P0_LO = E0^4*e^4*me^2*(hw0/hbar)*sqrt(kappa0)/(256*pi^3*hbar^6*eps0*kappa_LO*Omega^2); %Ma hoa dung Omega^2
        P0_PE = E0^4*e^4*me^2*cfg.material.piezo_kappa2*kBT*sqrt(kappa0)/(256*pi^3*hbar^7*kappa0*eps0*Omega^2);
        for it = 1:numel(td)
            eps_alpha = sp.E_J(td(it).initial);
            DeltaE = td(it).deltaE_J;
            transition_factor = abs(td(it).dipole_m)^2 * td(it).K;
            for ell = orders
                order_key = sprintf('order_%d', ell);
                order_factor = a0^(2*ell) / (2^(2*ell) * factorial(ell)^2); %Ma hoa dung a0^(2*ell)
                Pplus = P0_LO * transition_factor * order_factor * (N0 + 1.0) * compute_taylor_bessel(eps_alpha, sp.EF_J, DeltaE + hw0 - ell*Eph, ell, qd, cfg);
                Pminus = P0_LO * transition_factor * order_factor * N0 * compute_taylor_bessel(eps_alpha, sp.EF_J, DeltaE - hw0 - ell*Eph, ell, qd, cfg);
                spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + P0_PE*transition_factor*order_factor*compute_taylor_bessel(eps_alpha,sp.EF_J,DeltaE-ell*Eph,ell,qd,cfg);
                peak(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE);
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
    spectrum.total_raw = spectrum.mechanism.optical.total_raw+spectrum.mechanism.piezoelectric.total_raw;
    spectrum.mechanism.optical.total_plot = spectrum.mechanism.optical.total_raw;
    spectrum.mechanism.piezoelectric.total_plot = spectrum.mechanism.piezoelectric.total_raw;
    for ell = orders
        ok = sprintf('order_%d',ell);
        spectrum.mechanism.optical.(ok).emission_plot = spectrum.mechanism.optical.(ok).emission_raw;
        spectrum.mechanism.optical.(ok).absorption_plot = spectrum.mechanism.optical.(ok).absorption_raw;
        spectrum.mechanism.optical.(ok).total_plot = spectrum.mechanism.optical.(ok).total_raw;
        spectrum.mechanism.piezoelectric.(ok).total_plot = spectrum.mechanism.piezoelectric.(ok).total_raw;
    end
end
function Tval = compute_taylor_bessel(eps_alpha_J, EF_J, D_J, ell, qd, cfg)
    c = cfg.constants;
    hbar = c.hbar;
    kBT  = c.kB * cfg.temperature_K;
    me   = cfg.material.mstar_rel * c.m0;
    if D_J == 0
        Tval = NaN;
        return;
    end
    absD = abs(D_J);
    qd_power = qd^(2*ell - 2);
    Tval = 0.0;
    for s = ell:ell
        for v = 1:cfg.oap.series.v_max
            LambdaKsum = besselk(s-ell, v*absD/(2.0*kBT), 0);
            LambdaJ = (hbar^2 * qd^2 / (2.0 * me)) * (hbar^2 * qd^2 / (2.0 * me * absD))^(s-ell) * LambdaKsum;
            Tval = Tval + (1)^(s+v-1) * (s+1) * exp(-v*(eps_alpha_J+D_J/2.0-EF_J)/kBT) * qd_power * LambdaJ;
        end
    end
end
function metrics = extract_spectrum_metrics(spectrum, cfg)
% Peak position, FWHM and HWHM for each available component.
    E = spectrum.energy_meV;
    metrics.definition = 'profile_fwhm';
    for order = cfg.oap.photon_orders
        ok = sprintf('order_%d',order);
        optical = spectrum.mechanism.optical.(ok);
        metrics.optical.(ok).emission = profile_fwhm(E,optical.emission_raw);
        metrics.optical.(ok).absorption = profile_fwhm(E,optical.absorption_raw);
        metrics.optical.(ok).total = profile_fwhm(E,optical.total_raw);
        metrics.piezoelectric.(ok).total = profile_fwhm(E, ...
            spectrum.mechanism.piezoelectric.(ok).total_raw);
    end
end

function m = profile_fwhm(x, y)
% Width of a resolved peak at half its maximum.
    x = x(:).'; y = y(:).';
    % Keep resolved=false as the diagnostic flag, but report zero linewidth
    % when the profile does not contain a resolvable peak.
    m = struct('peak',NaN,'peak_x',NaN,'fwhm',0,'hwhm',0, ...
        'left_x',NaN,'right_x',NaN,'left_hwhm',0, ...
        'right_hwhm',0,'half_level',NaN,'resolved',false);
    if numel(x) ~= numel(y)
        error('QW:FWHM:SizeMismatch','x and y must have the same size.');
    end
    if numel(x) < 3 || ~isreal(x) || any(~isfinite(x)) || any(diff(x) <= 0)
        error('QW:FWHM:InvalidGrid','x must be finite and strictly increasing.');
    end
    if ~isreal(y) || any(~isfinite(y)) || any(y < 0), return; end

    [m.peak,ip] = max(y);
    if m.peak <= 0, return; end
    m.peak_x = x(ip); m.half_level = m.peak/2;
    h = m.half_level;
    if ip == 1 || ip == numel(y) || y(ip-1) < h || y(ip+1) < h, return; end
    il = find(y(1:ip) <= h,1,'last');
    ir = ip-1+find(y(ip:end) <= h,1,'first');
    if isempty(il) || isempty(ir) || il == ip || ir == ip, return; end
    m.left_x = x(il);%crossing(x,y,il,h);
    m.right_x = x(ir);%crossing(x,y,ir-1,h);
    m.left_hwhm = m.peak_x-m.left_x;
    m.right_hwhm = m.right_x-m.peak_x;
    m.fwhm = m.right_x-m.left_x;
    m.hwhm = m.fwhm/2;
    m.resolved = true;
end

function xc = crossing(x, y, i, h)
    if y(i+1) == y(i)
        xc = (x(i)+x(i+1))/2;
    else
        xc = x(i)+(h-y(i))*(x(i+1)-x(i))/(y(i+1)-y(i));
    end
end

function sweep_result = run_moap2_sweep(cfg, need_metrics, make_spectra, make_linewidth)
    values = cfg.sweep.values(:).';
    assert(~isempty(values) && all(isfinite(values)), ...
        'cfg.sweep.values must contain finite values.');
    cases = cell(size(values));
    labels = cell(size(values));
    for i = 1:numel(values)
        cfg_i = set_moap2_sweep_parameter(cfg,cfg.sweep.parameter,values(i));
        cases{i} = run_moap2_case(cfg_i,false,false,need_metrics);
        labels{i} = moap2_sweep_label(cfg.sweep.parameter,values(i));
        if cfg.run.verbose
            fprintf('[SWEEP] %s: %d/%d complete\n', ...
                cfg.sweep.parameter,i,numel(values));
        end
    end
    sweep_result = struct('parameter',cfg.sweep.parameter,'values',values, ...
        'labels',{labels},'cases',{cases});
    if need_metrics
        sweep_result.linewidth = collect_linewidths(cases,values,cfg);
        sweep_result.linewidth.reference = experimental_fwhm_references();
    end
    mechanisms = cfg.oap.mechanisms;
    if ischar(mechanisms), mechanisms = {mechanisms}; end
    if make_spectra
        figures = {}; names = {};
        for i = 1:numel(mechanisms)
            mk = lower(mechanisms{i});
            if ~ismember(mk,{'optical','piezoelectric'}), continue; end
            figures{end+1} = plot_moap2_sweep_spectra(cases,labels,cfg,mk); %#ok<AGROW>
            names{end+1} = ['sweep_spectra_' mk]; %#ok<AGROW>
        end
        save_moap2_sweep_spectra(sweep_result,figures,names,cfg);
    end
    if make_linewidth
        figures = {}; names = {};
        for i = 1:numel(mechanisms)
            mk = lower(mechanisms{i});
            if ~ismember(mk,{'optical','piezoelectric'}), continue; end
            figures{end+1} = plot_moap2_linewidth(sweep_result.linewidth,cfg,mk); %#ok<AGROW>
            names{end+1} = ['linewidth_sweep_' mk]; %#ok<AGROW>
        end
        save_moap2_linewidth(sweep_result,figures,names,cfg);
    end
end

function cfg = set_moap2_sweep_parameter(cfg, parameter, value)
    switch lower(parameter)
        case 'b_t', cfg.fields.B_T = value;
        case 't_k', cfg.temperature_K = value;
        case 'lz_nm', cfg.structure.Lz_nm = value;
        case 'e_kvcm', cfg.fields.E_kVcm = value;
        case 'nd_sheet_cm2', cfg.doping.Nd_sheet_cm2 = value;
        case 'manning_prefactor', cfg.structure.manning_prefactor = value;
        case 'u0_mev', cfg.structure.U0_meV = value;
        otherwise, error('Unsupported sweep parameter: %s',parameter);
    end
end

function label = moap2_sweep_label(parameter, value)
    switch lower(parameter)
        case 'b_t', label = sprintf('$B=%.3g\\,\\mathrm{T}$',value);
        case 't_k', label = sprintf('$T=%.3g\\,\\mathrm{K}$',value);
        case 'lz_nm', label = sprintf('$L_z=%.3g\\,\\mathrm{nm}$',value);
        case 'e_kvcm', label = sprintf('$E=%.3g\\,\\mathrm{kV/cm}$',value);
        case 'nd_sheet_cm2', label = sprintf('$N_D=%.3g\\,\\mathrm{cm}^{-2}$',value);
        case 'manning_prefactor', label = sprintf('$\\nu=%.3g$',value);
        case 'u0_mev', label = sprintf('$U_0=%.3g\\,\\mathrm{meV}$',value);
        otherwise, label = sprintf('%s = %.3g',strrep(parameter,'_','\_'),value);
    end
end

function linewidth = collect_linewidths(cases, values, cfg)
    linewidth.values = values(:).';
    for order = cfg.oap.photon_orders
        ok = sprintf('order_%d',order);
        for process = {'emission','absorption'}
            pk = process{1};
            linewidth.optical.(ok).([pk '_width']) = nan(size(values));
        end
        linewidth.piezoelectric.(ok).total_width = nan(size(values));
        for i = 1:numel(cases)
            for process = {'emission','absorption'}
                pk = process{1};
                linewidth.optical.(ok).([pk '_width'])(i) = ...
                    cases{i}.metrics.optical.(ok).(pk).fwhm;
            end
            linewidth.piezoelectric.(ok).total_width(i) = ...
                cases{i}.metrics.piezoelectric.(ok).total.fwhm;
        end
    end
end

function fig = plot_moap2_sweep_spectra(cases, labels, cfg, mechanism)
    idx = unique([1 floor((numel(cases)+1)/2) numel(cases)],'stable');
    colors = lines(numel(idx));
    [fig,ax] = new_moap2_axes(cfg,[100 100 710 520]);
    handles = gobjects(size(idx)); y_all = [];
    for j = 1:numel(idx)
        i = idx(j);
        E = cases{i}.spectrum.energy_meV(:);
        P = cases{i}.spectrum.mechanism.(mechanism).total_raw(:);
        handles(j) = plot(ax,E,P,'Color',colors(j,:),'LineWidth',1.65);
        y_all = [y_all; P(isfinite(P))]; %#ok<AGROW>
    end
    xlabel(ax,'Photon energy $\hbar\Omega$ (meV)','Interpreter','latex','FontSize',13);
    ylabel(ax,'MOAP (arb. units)','Interpreter','latex','FontSize',13);
    legend(ax,handles,labels(idx),'Location','best','Interpreter','latex');
    if ~isempty(y_all) && max(y_all) > 0, ylim(ax,[0 1.055*max(y_all)]); end
    format_moap2_axes(ax);
    suppress_moap_intensity_exponent(ax);
end

function fig = plot_moap2_linewidth(linewidth, cfg, mechanism)
    x = linewidth.values(:);
    colors = [0 0.3 0.8; 0.85 0.1 0.1; 0.7 0 0.65];
    markers = {'o','s','^'};
    [fig,ax] = new_moap2_axes(cfg,[100 100 700 510]);
    handles = gobjects(0); labels = {}; y_all = [];
    for j = 1:numel(cfg.oap.photon_orders)
        order = cfg.oap.photon_orders(j); ok = sprintf('order_%d',order);
        color = colors(mod(j-1,size(colors,1))+1,:);
        marker = markers{mod(j-1,numel(markers))+1};
        if strcmp(mechanism,'optical')
            processes = {'emission','absorption'};
        else
            processes = {'total'};
        end
        for p = 1:numel(processes)
            process = processes{p};
            y = linewidth.(mechanism).(ok).([process '_width'])(:);
            valid = isfinite(x) & isfinite(y);
            if ~any(valid), continue; end
            style = '-'; face = color; suffix = '';
            if strcmp(process,'absorption'), style = '--'; face = 'w'; suffix = ' absorption'; end
            if strcmp(process,'emission'), suffix = ' emission'; end
            handles(end+1) = plot(ax,x(valid),y(valid),style,'Color',color, ...
                'LineWidth',1.6,'Marker',marker,'MarkerSize',5, ...
                'MarkerEdgeColor',color,'MarkerFaceColor',face); %#ok<AGROW>
            labels{end+1} = sprintf('$\\ell=%d$%s',order,suffix); %#ok<AGROW>
            y_all = [y_all; y(valid)]; %#ok<AGROW>
        end
    end
    show_refs = cfg.output.show_experimental_fwhm && ...
        strcmp(mechanism,'optical') && strcmpi(cfg.sweep.parameter,'T_K');
    if show_refs
        refs = linewidth.reference;
        handles(end+1) = plot(ax,refs.unuma.T_K,refs.unuma.W_meV,':*', ...
            'Color',[0 .48 .1],'LineWidth',1.2);
        handles(end+1) = plot(ax,refs.von_allmen.T_K,refs.von_allmen.W_meV,':*', ...
            'Color',[.9 .38 .05],'LineWidth',1.2);
        labels(end+1:end+2) = {'T. Unuma et al. [54,55] (digitized)', ...
            'P. Von Allmen et al. [56] (digitized)'};
        y_all = [y_all; refs.unuma.W_meV(:); refs.von_allmen.W_meV(:)];
    end
    xlabel(ax,moap2_sweep_axis_label(cfg.sweep.parameter), ...
        'Interpreter','latex','FontSize',13);
    ylabel(ax,'FWHM (meV)','Interpreter','latex','FontSize',13);
    if ~isempty(handles)
        legend(ax,handles,labels,'Location','best','Interpreter','latex', ...
            'FontSize',9.5,'Box','on');
    end
    if ~isempty(y_all) && max(y_all) > 0, ylim(ax,[0 1.06*max(y_all)]); end
    if show_refs
        xlim(ax,[min([x;refs.unuma.T_K(:);refs.von_allmen.T_K(:)]), ...
            max([x;refs.unuma.T_K(:);refs.von_allmen.T_K(:)])]);
    elseif isscalar(x)
        dx = max(0.05*abs(x),1); xlim(ax,[x-dx x+dx]);
    elseif numel(x) > 1
        xlim(ax,[min(x) max(x)]);
    end
    format_moap2_axes(ax);
end

function refs = experimental_fwhm_references()
    refs.unuma.T_K = [5 30 75 125 167 217 257 296];
    refs.unuma.W_meV = [11.03 11.03 11.43 11.93 13.12 13.67 14.13 14.83];
    refs.von_allmen.T_K = [5 15 30 40 50 60 75 97 116 135 155 174 194 213 233 252 272 292];
    refs.von_allmen.W_meV = [10.2 10.2 10.2 10.1 10.1 10.1 10.21 ...
        10.35 10.71 10.88 11.40 11.93 12.70 14.36 16.04 15.65 16.98 18.59];
end

function [fig,ax] = new_moap2_axes(cfg, position)
    fig = figure('Visible',cfg.output.visible,'Color','w', ...
        'Units','pixels','Position',position);
    ax = axes(fig); hold(ax,'on');
end

function label = moap2_sweep_axis_label(parameter)
    switch lower(parameter)
        case 'b_t', label = 'Magnetic field $B$ (T)';
        case 't_k', label = 'Temperature $T$ (K)';
        case 'lz_nm', label = 'Well width $L_z$ (nm)';
        case 'e_kvcm', label = 'Electric field $E$ (kV/cm)';
        case 'nd_sheet_cm2', label = 'Sheet density $N_D$ (cm$^{-2}$)';
        case 'manning_prefactor', label = 'Manning parameter $\nu$';
        case 'u0_mev', label = 'Potential depth $U_0$ (meV)';
        otherwise, label = strrep(parameter,'_','\_');
    end
end

function format_moap2_axes(ax)
    set(ax,'Box','on','Layer','top','LineWidth',0.9, ...
        'FontName','Times New Roman','FontSize',11,'TickDir','in', ...
        'TickLength',[0.018 0.018],'XMinorTick','on','YMinorTick','on', ...
        'TickLabelInterpreter','latex');
    grid(ax,'on'); ax.GridLineStyle = '--'; ax.GridAlpha = 0.2;
    ax.XMinorGrid = 'off'; ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset,0.02);
end

function suppress_moap_intensity_exponent(ax)
    tick_labels = ax.YTickLabel;
    ax.YAxis.Exponent = 0;
    ax.YTickLabel = tick_labels;
end

function save_moap2_sweep_spectra(sweep_result, figures, names, cfg)
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory,'sweep_spectra.mat'), ...
            'sweep_result','-v7.3');
    end
    if cfg.output.save_csv
        E = sweep_result.cases{1}.spectrum.energy_meV(:);
        T = table(E,'VariableNames',{'PhotonEnergy_meV'});
        for i = 1:numel(sweep_result.cases)
            for mechanism = {'optical','piezoelectric'}
                mk = mechanism{1};
                name = matlab.lang.makeValidName(sprintf('%s_%s_%g', ...
                    mk,cfg.sweep.parameter,sweep_result.values(i)));
                T.(name) = sweep_result.cases{i}.spectrum.mechanism.(mk).total_raw(:);
            end
        end
        writetable(T,fullfile(cfg.output.directory,'sweep_spectra.csv'));
    end
    save_moap2_figures(figures,names,'',cfg);
end

function save_moap2_linewidth(sweep_result, figures, names, cfg)
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory,'linewidth_sweep.mat'), ...
            'sweep_result','-v7.3');
    end
    if cfg.output.save_csv
        tables = cell(1,numel(sweep_result.cases));
        for i = 1:numel(sweep_result.cases)
            tables{i} = linewidth_metrics_table(sweep_result.cases{i}.metrics, ...
                cfg,sweep_result.values(i));
        end
        writetable(vertcat(tables{:}), ...
            fullfile(cfg.output.directory,'linewidth_sweep.csv'));
        if strcmpi(cfg.sweep.parameter,'T_K')
            export_experimental_fwhm_csv(sweep_result.linewidth.reference, ...
                fullfile(cfg.output.directory,'experimental_fwhm_digitized.csv'));
        end
    end
    save_moap2_figures(figures,names,'',cfg);
end

function export_experimental_fwhm_csv(refs, filename)
    temperature = [refs.unuma.T_K(:);refs.von_allmen.T_K(:)];
    fwhm = [refs.unuma.W_meV(:);refs.von_allmen.W_meV(:)];
    source = [repmat({'T. Unuma et al. [54,55]'},numel(refs.unuma.T_K),1); ...
        repmat({'P. Von Allmen et al. [56]'},numel(refs.von_allmen.T_K),1)];
    writetable(table(source,temperature,fwhm,'VariableNames', ...
        {'ExperimentalReference','Temperature_K','FWHM_meV'}),filename);
end

function write_fwhm_metrics_csv(metrics, cfg, filename, sweep_value)
    writetable(linewidth_metrics_table(metrics,cfg,sweep_value),filename);
end

function T = linewidth_metrics_table(metrics, cfg, sweep_value)
    parameter_value = []; mechanism = {}; photon_order = [];
    process = {}; peak = []; left_x = []; right_x = [];
    fwhm = []; hwhm = []; resolved = [];
    for order = cfg.oap.photon_orders
        ok = sprintf('order_%d',order);
        entries = {{'optical','emission',metrics.optical.(ok).emission}, ...
            {'optical','absorption',metrics.optical.(ok).absorption}, ...
            {'piezoelectric','total',metrics.piezoelectric.(ok).total}};
        for i = 1:numel(entries)
            item = entries{i}{3};
            parameter_value(end+1,1) = sweep_value; %#ok<AGROW>
            mechanism{end+1,1} = entries{i}{1}; %#ok<AGROW>
            photon_order(end+1,1) = order; %#ok<AGROW>
            process{end+1,1} = entries{i}{2}; %#ok<AGROW>
            peak(end+1,1) = item.peak_x; %#ok<AGROW>
            left_x(end+1,1) = item.left_x; %#ok<AGROW>
            right_x(end+1,1) = item.right_x; %#ok<AGROW>
            fwhm(end+1,1) = item.fwhm; %#ok<AGROW>
            hwhm(end+1,1) = item.hwhm; %#ok<AGROW>
            resolved(end+1,1) = item.resolved; %#ok<AGROW>
        end
    end
    T = table(parameter_value,mechanism,photon_order,process,peak, ...
        left_x,right_x,fwhm,hwhm,logical(resolved), ...
        'VariableNames',{'SweepValue','Mechanism','PhotonOrder','Process', ...
        'Peak_meV','LeftHalfMaximum_meV','RightHalfMaximum_meV', ...
        'FWHM_meV','HWHM_meV','Resolved'});
end

function save_moap2_figures(figures, names, prefix, cfg)
    if ~cfg.output.save_figures, return; end
    for i = 1:min(numel(figures),numel(names))
        if ~ishghandle(figures{i}), continue; end
        stem = names{i};
        if ~isempty(prefix), stem = sprintf('%s_%s',prefix,stem); end
        filename = fullfile(cfg.output.directory, ...
            [stem '.' cfg.output.figure_format]);
        if exist('exportgraphics','file') == 2
            exportgraphics(figures{i},filename,'Resolution',cfg.output.figure_dpi);
        else
            set(figures{i},'PaperPositionMode','auto');
            print(figures{i},filename,'-dpng',sprintf('-r%d',cfg.output.figure_dpi));
        end
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
    hPotential = plot( ax, z, Ueff, '-', 'Color', potential_color, 'LineWidth', 2.0, 'DisplayName', '$U_{\mathrm{eff}}(z)$');
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
    format_moap2_axes(ax);
    ax.GridLineStyle = ':';
end
function fig = plot_moap_contributions(spectrum, cfg, mechanism)
    mk = lower(mechanism);
    E = spectrum.energy_meV(:);
    photon_orders = cfg.oap.photon_orders(:).';
    fig = figure( 'Visible', cfg.output.visible, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 720 520]);
    ax = axes(fig);
    hold(ax, 'on');
    order_colors = [ 0.000 0.300 0.800; 0.850 0.100 0.100; 0.100 0.550 0.100; 0.650 0.100 0.650; 0.100 0.600 0.650 ];
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
        curve_color = order_colors(color_idx, :);
        if strcmpi(mk, 'piezoelectric')
            curves = {spectrum.mechanism.(mk).(ok).total_plot(:)};
            line_styles = {'-'};
            curve_labels = {sprintf('$\\ell=%d$', order)};
        else
            curves = {spectrum.mechanism.(mk).(ok).emission_plot(:), spectrum.mechanism.(mk).(ok).absorption_plot(:)};
            line_styles = {'-', '--'};
            curve_labels = {sprintf('$\\ell=%d$ emission', order), sprintf('$\\ell=%d$ absorption', order)};
        end
        for ic = 1:numel(curves)
            P = curves{ic};
            h = plot(ax, E, P, line_styles{ic}, 'Color', curve_color, 'LineWidth', 1.65);
            legend_handles(end+1) = h; %#ok<AGROW>
            labels{end+1} = curve_labels{ic}; %#ok<AGROW>
            P_abs = abs(P);
            Pmax = max(P_abs(isfinite(P_abs)));
            if ~isempty(Pmax) && isfinite(Pmax) && Pmax > 0
                idx_active = isfinite(P) & P_abs >= 1.0e-4*Pmax;
                x_active_all = [x_active_all; E(idx_active)]; %#ok<AGROW>
            end
            y_all = [y_all; P(isfinite(P))]; %#ok<AGROW>
        end
    end
    xlabel( ax, 'Photon energy $\hbar\Omega$ (meV)', 'Interpreter', 'latex', 'FontSize', 13);

    ylabel( ax, 'MOAP (arb. units)', 'Interpreter', 'latex', 'FontSize', 13);

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
    format_moap2_axes(ax);
    suppress_moap_intensity_exponent(ax);
end
