function results = moap3(cfg)
    close all;
    if nargin == 0 || isempty(cfg)
        cfg.run.task = 'single'; % single | sweep_spectra | linewidth_sweep | all_demo
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
        cfg.fields.B_T = 0;
        cfg.fields.E_kVcm = 25;
        cfg.fields.kx_inv_m = 0;
        cfg.temperature_K = 100;
        cfg.sp.max_outer = 100;
        cfg.sp.mix = 0.5;
        cfg.sp.tol_VH_relative = 1e-3;
        cfg.sp.tol_EF_meV = 1e-3;
        cfg.laser.E0_kVcm = 4.5;
        cfg.laser.a0_nm = 7.5;  %
        cfg.oap.transitions = [1 2];
        cfg.oap.photon_orders = [1 2];
        cfg.oap.mechanisms = {'optical','piezoelectric'};
        cfg.oap.photon_energy_meV = 5:.001:100;
        cfg.oap.qz_max_inv_nm = 2;
        cfg.oap.Nqz = 101;
        cfg.oap.plot_full_energy_range = true;
        cfg.oap.series.s_max = 2;
        cfg.oap.series.v_max = 5;
        cfg.sweep.parameter = 'T_K';
        cfg.sweep.values = [50 100 150];
        cfg.linewidth.core_detuning_meV = 1e-4;
        cfg.linewidth.foot_fraction = 1e-3;
        cfg.linewidth.max_search_meV = 40;
        cfg.linewidth.root_tolerance_meV = 1e-8;
        cfg.output.directory = fullfile(fileparts(mfilename('fullpath')), 'results_moap3');
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
        cfg.constants.meV = cfg.constants.e*1e-3;
        cfg.constants.nm = 1e-9;
    end
    cfg = apply_defaults(cfg);
    if ~exist(cfg.output.directory, 'dir')
        mkdir(cfg.output.directory);
    end
    task = lower(cfg.run.task);
    results = struct('config',cfg);
    switch task
        case 'single'
            results.single = run_case(cfg, true, true);
        case 'sweep_spectra'
            results.sweep = run_sweep(cfg, false, true, false);
        case 'linewidth_sweep'
            results.linewidth = run_sweep(cfg, true, false, true);
        case 'all_demo'
            results.single = run_case(cfg, true, true);
            results.sweep = run_sweep(cfg, true, true, true);
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

function cfg = apply_defaults(cfg)
    if ~isfield(cfg,'run'), cfg.run = struct(); end
    if ~isfield(cfg.run,'task'), cfg.run.task = 'single'; end
    if ~isfield(cfg.run,'verbose'), cfg.run.verbose = true; end
    if ~isfield(cfg.oap,'mechanisms')
        cfg.oap.mechanisms = {'optical','piezoelectric'};
    end
    if ~isfield(cfg,'sweep'), cfg.sweep = struct(); end
    if ~isfield(cfg.sweep,'parameter'), cfg.sweep.parameter = 'T_K'; end
    if ~isfield(cfg.sweep,'values'), cfg.sweep.values = cfg.temperature_K; end
    if ~isfield(cfg,'linewidth'), cfg.linewidth = struct(); end
    defaults = struct('core_detuning_meV',1e-4,'foot_fraction',1e-3, ...
        'max_search_meV',40,'root_tolerance_meV',1e-8);
    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(cfg.linewidth,name)
            cfg.linewidth.(name) = defaults.(name);
        end
    end
    assert(cfg.linewidth.core_detuning_meV > 0, ...
        'linewidth.core_detuning_meV must be positive.');
    assert(cfg.linewidth.foot_fraction > 0 && ...
        cfg.linewidth.foot_fraction < 1, ...
        'linewidth.foot_fraction must lie strictly between zero and one.');
    assert(cfg.linewidth.max_search_meV > cfg.linewidth.core_detuning_meV, ...
        'linewidth.max_search_meV must exceed core_detuning_meV.');
end

function case_result = run_case(cfg, make_plots, save_outputs, compute_metrics)
    if nargin < 4, compute_metrics = true; end
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
        metrics = extract_foot_metrics(sp,transition_data,cfg);
    else
        metrics = struct();
    end
    case_result = struct('cfg',cfg,'sp',sp,'transitions',transition_data, ...
        'spectrum',spectrum,'metrics',metrics);
    figures = {};
    names = {};
    if make_plots
        figures = {plot_electronic_structure(sp,cfg), ...
            plot_moap_spectrum(spectrum,cfg,'optical'), ...
            plot_moap_spectrum(spectrum,cfg,'piezoelectric'), ...
            plot_moap_contributions(spectrum,cfg,'optical'), ...
            plot_moap_contributions(spectrum,cfg,'piezoelectric')};
        names = {'structure','spectrum_optical','spectrum_piezoelectric', ...
            'contributions_optical','contributions_piezoelectric'};
    end
    if save_outputs
        if cfg.output.save_mat
            save(fullfile(cfg.output.directory,'single.mat'),'case_result','-v7.3');
        end
        if cfg.output.save_csv
            T = table(spectrum.energy_meV(:),spectrum.total_raw(:), ...
                spectrum.mechanism.optical.total_raw(:), ...
                spectrum.mechanism.piezoelectric.total_raw(:), ...
                'VariableNames',{'PhotonEnergy_meV','TotalRaw','OpticalRaw','PiezoelectricRaw'});
            writetable(T,fullfile(cfg.output.directory,'single_spectrum.csv'));
            writetable(foot_metrics_table(metrics,cfg,NaN), ...
                fullfile(cfg.output.directory,'single_metrics.csv'));
        end
        save_figures(figures,names,'single',cfg);
    end
    if cfg.run.verbose
        fprintf('[MOAP] Converged: %d in %d iterations\n',sp.converged,sp.iterations);
        fprintf('[MOAP] EF = %.8f meV\n',sp.EF_meV);
        fprintf('[MOAP] Subband energies (meV): '); fprintf('%.8f ',sp.E_meV);
        fprintf('\n[MOAP] Max raw OAP = %.6e\n',max(spectrum.total_raw));
    end
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

    qd = sqrt(sp.Nd_sheet_m2*c.e^2/(cfg.structure.Lz_nm*c.nm* ...
        cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
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
        Omega = Eph/hbar;
        P0_LO = E0^4*e^4*me^2*(hw0/hbar)*sqrt(kappa0)/ ...
            (256*pi^3*hbar^6*eps0*kappa_LO*Omega^2);
        P0_PE = E0^4*e^4*me^2*cfg.material.piezo_kappa2*kBT*sqrt(kappa0)/ ...
            (256*pi^3*hbar^7*kappa0*eps0*Omega^2);
        for it = 1:numel(td)
            eps_alpha = sp.E_J(td(it).initial);
            DeltaE = td(it).deltaE_J;
            transition_factor = abs(td(it).dipole_m)^2 * td(it).K;
            for ell = 2
                order_key = sprintf('order_%d',ell);
                order_factor = a0^(ell)/(2^(2*ell)*factorial(ell)^2);
                Pplus = P0_LO * transition_factor * order_factor * (N0 + 1.0) * compute_taylor_bessel(eps_alpha, sp.EF_J, DeltaE + hw0 - ell*Eph, ell, qd, cfg);
                Pminus = P0_LO * transition_factor * order_factor * N0 * compute_taylor_bessel(eps_alpha, sp.EF_J, DeltaE - hw0 - ell*Eph, ell, qd, cfg);
                spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + P0_PE*transition_factor*order_factor*compute_taylor_bessel(eps_alpha,sp.EF_J,DeltaE-ell*Eph,ell,qd,cfg);
                em(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE);
                ab(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE);
                pie(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE);
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
    spectrum.total_raw = spectrum.mechanism.optical.total_raw + ...
        spectrum.mechanism.piezoelectric.total_raw;
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
    for s = 0:cfg.oap.series.s_max
        for v = 1:cfg.oap.series.v_max
            LambdaKsum = besselk(abs(s-ell), v*absD/(2.0*kBT), 0)*0 + besselk(abs(s-ell+1), v*absD/(2.0*kBT), 0);
            LambdaJ = (hbar^2 * qd^2 / (2.0 * me)) * (hbar^2 * qd^2 / (2.0 * me * absD))^(s-ell) * LambdaKsum;
            Tval = Tval + (-1)^(s+v-1) * (s+1) * exp(-v*(eps_alpha_J+D_J/2.0-EF_J)/kBT) * qd_power * LambdaJ;
        end
    end
end
function metrics = extract_foot_metrics(sp, td, cfg)
%EXTRACT_FOOT_METRICS Near-base linewidth on the continuous Van profile.
% The exact Van resonance is singular, so no sampled peak height is used.
% A finite reference level is evaluated at a fixed detuning from the
% analytical resonance center.  The two foot crossings are then obtained
% with root finding.  Consequently the result does not depend on the
% spacing or phase of cfg.oap.photon_energy_meV.
    metrics = struct();
    metrics.definition = 'foot_to_foot';
    metrics.core_detuning_meV = cfg.linewidth.core_detuning_meV;
    metrics.foot_fraction = cfg.linewidth.foot_fraction;
    for order = cfg.oap.photon_orders
        ok = sprintf('order_%d',order);
        metrics.optical.(ok).emission = resonance_foot_metric( ...
            sp,td,cfg,'optical',order,'emission');
        metrics.optical.(ok).absorption = resonance_foot_metric( ...
            sp,td,cfg,'optical',order,'absorption');
        metrics.optical.(ok).total = empty_foot_metric(cfg);
        metrics.piezoelectric.(ok).emission = empty_foot_metric(cfg);
        metrics.piezoelectric.(ok).absorption = empty_foot_metric(cfg);
        metrics.piezoelectric.(ok).total = resonance_foot_metric( ...
            sp,td,cfg,'piezoelectric',order,'total');
    end
end

function m = resonance_foot_metric(sp, td, cfg, mechanism, order, process)
    m = empty_foot_metric(cfg);
    deltaE = [td.deltaE_J]/cfg.constants.meV;
    if strcmp(mechanism,'piezoelectric')
        centers = deltaE/order;
    elseif strcmp(process,'emission')
        centers = (deltaE+cfg.material.LO_phonon_meV)/order;
    else
        centers = (deltaE-cfg.material.LO_phonon_meV)/order;
    end
    centers = centers(isfinite(centers) & centers > 0);
    d0 = cfg.linewidth.core_detuning_meV;
    strengths = nan(size(centers));
    for i = 1:numel(centers)
        if centers(i) <= d0, continue; end
        y_left = component_value(sp,td,cfg,mechanism,order,process,centers(i)-d0);
        y_right = component_value(sp,td,cfg,mechanism,order,process,centers(i)+d0);
        candidates = [y_left y_right];
        candidates = candidates(isfinite(candidates) & candidates > 0);
        if ~isempty(candidates), strengths(i) = max(candidates); end
    end
    [reference_peak,idx] = max(strengths);
    if isempty(idx) || ~isfinite(reference_peak) || reference_peak <= 0
        return;
    end
    center = centers(idx);
    foot_level = cfg.linewidth.foot_fraction*reference_peak;
    left_delta = find_foot_delta(sp,td,cfg,mechanism,order,process, ...
        center,-1,foot_level);
    right_delta = find_foot_delta(sp,td,cfg,mechanism,order,process, ...
        center,+1,foot_level);
    m.reference_peak = reference_peak;
    m.peak_x = center;
    m.foot_level = foot_level;
    if ~isfinite(left_delta) || ~isfinite(right_delta), return; end
    m.left_x = center-left_delta;
    m.right_x = center+right_delta;
    m.left_width = left_delta;
    m.right_width = right_delta;
    m.width = left_delta+right_delta;
    m.resolved = isfinite(m.width) && m.width > 0;
end

function delta = find_foot_delta(sp,td,cfg,mechanism,order,process,center,side,level)
    delta = NaN;
    d0 = cfg.linewidth.core_detuning_meV;
    dmax = cfg.linewidth.max_search_meV;
    if side < 0, dmax = min(dmax,center*(1-1e-10)); end
    if dmax <= d0, return; end
    objective = @(d) foot_objective(d,sp,td,cfg,mechanism,order,process, ...
        center,side,level);
    a = d0;
    fa = objective(a);
    if ~isfinite(fa) || fa <= 0, return; end
    % A modest geometric expansion finds the first outward crossing.  A
    % factor of two can jump across the shallow foot of an asymmetric Van
    % profile and land on its far tail after the profile has risen again.
    expansion = 1.25;
    b = min(expansion*a,dmax);
    fb = objective(b);
    while isfinite(fb) && fb > 0 && b < dmax
        a = b;
        b = min(expansion*b,dmax);
        fb = objective(b);
    end
    if ~isfinite(fb) || fb > 0, return; end
    options = optimset('TolX',cfg.linewidth.root_tolerance_meV,'Display','off');
    delta = fzero(objective,[a b],options);
end

function value = foot_objective(delta,sp,td,cfg,mechanism,order,process,center,side,level)
    energy = center+side*delta;
    y = component_value(sp,td,cfg,mechanism,order,process,energy);
    if ~isfinite(y) || y <= 0
        % A sign change in the truncated series lies below any positive
        % foot level.  Keep the objective finite so fzero can bracket the
        % positive crossing immediately before that point.
        value = -1e6;
    else
        value = log(y)-log(level);
    end
end

function value = component_value(sp,td,cfg,mechanism,order,process,energy_meV)
    cfg_eval = cfg;
    cfg_eval.oap.photon_energy_meV = energy_meV;
    s = compute_moap(sp,td,cfg_eval);
    ok = sprintf('order_%d',order);
    if strcmp(mechanism,'piezoelectric')
        value = s.mechanism.piezoelectric.(ok).total_raw;
    else
        value = s.mechanism.optical.(ok).([process '_raw']);
    end
end

function m = empty_foot_metric(cfg)
    m = struct('reference_peak',NaN,'peak_x',NaN,'width',NaN, ...
        'left_x',NaN,'right_x',NaN,'left_width',NaN,'right_width',NaN, ...
        'foot_level',NaN,'core_detuning_meV', ...
        cfg.linewidth.core_detuning_meV,'foot_fraction', ...
        cfg.linewidth.foot_fraction,'resolved',false);
end

function sweep_result = run_sweep(cfg, need_metrics, make_spectra, make_linewidth)
    values = cfg.sweep.values(:).';
    assert(~isempty(values) && all(isfinite(values)), ...
        'cfg.sweep.values must contain finite values.');
    cases = cell(size(values));
    labels = cell(size(values));
    parameter = lower(cfg.sweep.parameter);
    for i = 1:numel(values)
        cfg_i = cfg;
        value = values(i);
        switch parameter
            case 'b_t', cfg_i.fields.B_T = value; labels{i} = sprintf('$B=%.3g\\,\\mathrm{T}$',value);
            case 't_k', cfg_i.temperature_K = value; labels{i} = sprintf('$T=%.3g\\,\\mathrm{K}$',value);
            case 'lz_nm', cfg_i.structure.Lz_nm = value; labels{i} = sprintf('$L_z=%.3g\\,\\mathrm{nm}$',value);
            case 'e_kvcm', cfg_i.fields.E_kVcm = value; labels{i} = sprintf('$E=%.3g\\,\\mathrm{kV/cm}$',value);
            case 'nd_sheet_cm2', cfg_i.doping.Nd_sheet_cm2 = value; labels{i} = sprintf('$N_D=%.3g\\,\\mathrm{cm}^{-2}$',value);
            case 'manning_prefactor', cfg_i.structure.manning_prefactor = value; labels{i} = sprintf('$\\nu=%.3g$',value);
            case 'u0_mev', cfg_i.structure.U0_meV = value; labels{i} = sprintf('$U_0=%.3g\\,\\mathrm{meV}$',value);
            otherwise, error('Unsupported sweep parameter: %s',cfg.sweep.parameter);
        end
        cases{i} = run_case(cfg_i,false,false,need_metrics);
        if cfg.run.verbose
            fprintf('[SWEEP] %s: %d/%d complete\n', ...
                cfg.sweep.parameter,i,numel(values));
        end
    end
    sweep_result = struct('parameter',cfg.sweep.parameter,'values',values, ...
        'labels',{labels},'cases',{cases});
    if need_metrics
        linewidth.values = values;
        for order = cfg.oap.photon_orders
            key = sprintf('order_%d',order);
            for process = {'emission','absorption'}
                linewidth.optical.(key).([process{1} '_width']) = nan(size(values));
            end
            linewidth.piezoelectric.(key).total_width = nan(size(values));
            for i = 1:numel(cases)
                linewidth.optical.(key).emission_width(i) = ...
                    cases{i}.metrics.optical.(key).emission.width;
                linewidth.optical.(key).absorption_width(i) = ...
                    cases{i}.metrics.optical.(key).absorption.width;
                linewidth.piezoelectric.(key).total_width(i) = ...
                    cases{i}.metrics.piezoelectric.(key).total.width;
            end
        end
        sweep_result.linewidth = linewidth;
    end
    mechanisms = cfg.oap.mechanisms;
    if ischar(mechanisms), mechanisms = {mechanisms}; end
    if make_spectra
        figures = {}; names = {};
        for i = 1:numel(mechanisms)
            mk = lower(mechanisms{i});
            if ~ismember(mk,{'optical','piezoelectric'}), continue; end
            figures{end+1} = plot_sweep_spectra(cases,labels,cfg,mk); %#ok<AGROW>
            names{end+1} = ['sweep_spectra_' mk]; %#ok<AGROW>
        end
        if cfg.output.save_mat
            save(fullfile(cfg.output.directory,'sweep_spectra.mat'), ...
                'sweep_result','-v7.3');
        end
        if cfg.output.save_csv
            T = table(cases{1}.spectrum.energy_meV(:), ...
                'VariableNames',{'PhotonEnergy_meV'});
            for i = 1:numel(cases)
                for mechanism = {'optical','piezoelectric'}
                    mk = mechanism{1};
                    name = matlab.lang.makeValidName(sprintf('%s_%s_%g', ...
                        mk,cfg.sweep.parameter,values(i)));
                    T.(name) = cases{i}.spectrum.mechanism.(mk).total_raw(:);
                end
            end
            writetable(T,fullfile(cfg.output.directory,'sweep_spectra.csv'));
        end
        save_figures(figures,names,'',cfg);
    end
    if make_linewidth
        figures = {}; names = {};
        for i = 1:numel(mechanisms)
            mk = lower(mechanisms{i});
            if ~ismember(mk,{'optical','piezoelectric'}), continue; end
            figures{end+1} = plot_linewidth(sweep_result.linewidth,cfg,mk); %#ok<AGROW>
            names{end+1} = ['linewidth_sweep_' mk]; %#ok<AGROW>
        end
        if cfg.output.save_mat
            save(fullfile(cfg.output.directory,'linewidth_sweep.mat'), ...
                'sweep_result','-v7.3');
        end
        if cfg.output.save_csv
            tables = cell(1,numel(cases));
            for i = 1:numel(cases)
                tables{i} = foot_metrics_table(cases{i}.metrics,cfg,values(i));
            end
            writetable(vertcat(tables{:}), ...
                fullfile(cfg.output.directory,'linewidth_sweep.csv'));
        end
        save_figures(figures,names,'',cfg);
    end
end

function fig = plot_sweep_spectra(cases, labels, cfg, mechanism)
    selected = unique([1,floor((numel(cases)+1)/2),numel(cases)],'stable');
    colors = [0.85 0 0.85; 0.05 0.15 0.65; 0.60 0.05 0];
    fig = figure('Visible',cfg.output.visible,'Color','w','Units','pixels', ...
        'Position',[100 100 710 520]);
    ax = axes(fig); hold(ax,'on');
    handles = gobjects(size(selected));
    ymax = 0;
    for j = 1:numel(selected)
        i = selected(j);
        E = cases{i}.spectrum.energy_meV;
        P = cases{i}.spectrum.mechanism.(mechanism).total_raw;
        handles(j) = plot(ax,E,P,'Color',colors(j,:),'LineWidth',1.65);
        ymax = max(ymax,max(P(isfinite(P))));
    end
    xlabel(ax,'Photon energy $\hbar\Omega$ (meV)','Interpreter','latex');
    ylabel(ax,'MOAP (arb. units)','Interpreter','latex');
    legend(ax,handles,labels(selected),'Location','best','Interpreter','latex');
    if ymax > 0, ylim(ax,[0 1.055*ymax]); end
    format_axes(ax);
end
function fig = plot_linewidth(linewidth, cfg, mechanism)
    x = linewidth.values(:);
    colors = [0 0.3 0.8; 0.85 0.1 0.1; 0.7 0 0.65];
    markers = {'o','s','^'};
    fig = figure('Visible',cfg.output.visible,'Color','w','Units','pixels', ...
        'Position',[100 100 700 510]);
    ax = axes(fig); hold(ax,'on');
    handles = gobjects(2*numel(cfg.oap.photon_orders),1);
    labels = cell(size(handles));
    count = 0;
    for j = 1:numel(cfg.oap.photon_orders)
        ell = cfg.oap.photon_orders(j);
        key = sprintf('order_%d',ell);
        if strcmp(mechanism,'optical'), processes = {'emission','absorption'};
        else, processes = {'total'}; end
        for k = 1:numel(processes)
            process = processes{k};
            y = linewidth.(mechanism).(key).([process '_width'])(:);
            valid = isfinite(x) & isfinite(y);
            if ~any(valid), continue; end
            color = colors(mod(j-1,size(colors,1))+1,:);
            marker = markers{mod(j-1,numel(markers))+1};
            style = '-'; face = color; suffix = '';
            if strcmp(process,'absorption')
                style = '--'; face = 'w'; suffix = ' absorption';
            elseif strcmp(process,'emission')
                suffix = ' emission';
            end
            count = count+1;
            handles(count) = plot(ax,x(valid),y(valid),style, ...
                'Color',color,'LineWidth',1.6,'Marker',marker, ...
                'MarkerFaceColor',face);
            labels{count} = sprintf('$\\ell=%d$%s',ell,suffix);
        end
    end
    switch lower(cfg.sweep.parameter)
        case 'b_t', axis_label = 'Magnetic field $B$ (T)';
        case 't_k', axis_label = 'Temperature $T$ (K)';
        case 'lz_nm', axis_label = 'Well width $L_z$ (nm)';
        case 'e_kvcm', axis_label = 'Electric field $E$ (kV/cm)';
        case 'nd_sheet_cm2', axis_label = 'Sheet density $N_D$ (cm$^{-2}$)';
        case 'manning_prefactor', axis_label = 'Manning parameter $\nu$';
        case 'u0_mev', axis_label = 'Potential depth $U_0$ (meV)';
        otherwise, axis_label = strrep(cfg.sweep.parameter,'_','\_');
    end
    xlabel(ax,axis_label,'Interpreter','latex');
    ylabel(ax,'Foot-to-foot linewidth (meV)','Interpreter','latex');
    if count > 0
        legend(ax,handles(1:count),labels(1:count),'Location','best','Interpreter','latex');
    end
    if isscalar(x), xlim(ax,x+[-1 1]*max(0.05*abs(x),1));
    elseif numel(x)>1, xlim(ax,[min(x) max(x)]); end
    format_axes(ax);
end
function format_axes(ax)
    set(ax,'Box','on','Layer','top','LineWidth',0.9, ...
        'FontName','Times New Roman','FontSize',11,'TickDir','in', ...
        'TickLength',[0.018 0.018],'XMinorTick','on','YMinorTick','on', ...
        'TickLabelInterpreter','latex');
    grid(ax,'on'); ax.GridLineStyle = '--'; ax.GridAlpha = 0.2;
    ax.XMinorGrid = 'off'; ax.YMinorGrid = 'off';
    ax.LooseInset = max(ax.TightInset,0.02);
end

function T = foot_metrics_table(metrics, cfg, sweep_value)
    parameter_value = []; mechanism = {}; photon_order = [];
    process = {}; resonance = []; left_foot = []; right_foot = [];
    linewidth = []; resolved = [];
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
            resonance(end+1,1) = item.peak_x; %#ok<AGROW>
            left_foot(end+1,1) = item.left_x; %#ok<AGROW>
            right_foot(end+1,1) = item.right_x; %#ok<AGROW>
            linewidth(end+1,1) = item.width; %#ok<AGROW>
            resolved(end+1,1) = item.resolved; %#ok<AGROW>
        end
    end
    T = table(parameter_value,mechanism,photon_order,process,resonance, ...
        left_foot,right_foot,linewidth,logical(resolved), ...
        'VariableNames',{'SweepValue','Mechanism','PhotonOrder','Process', ...
        'Resonance_meV','LeftFoot_meV','RightFoot_meV','Linewidth_meV','Resolved'});
end

function save_figures(figures, names, prefix, cfg)
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
    fig = figure('Visible',cfg.output.visible,'Color','w','Units','pixels', ...
        'Position',[100 100 720 520]);
    ax = axes(fig); hold(ax,'on');
    z = sp.z_nm(:);
    E = sp.E_meV(:);
    nplot = min([4,numel(E),size(sp.Psi,2)]);
    if nplot > 1
        spacing = diff(sort(E(1:nplot)));
        spacing = spacing(isfinite(spacing) & spacing > 1e-9);
    else
        spacing = [];
    end
    if isempty(spacing)
        scale = 0.08*(max(sp.Veff_meV)-min(sp.Veff_meV));
    else
        scale = 0.35*median(spacing);
    end
    if ~isfinite(scale) || scale <= 0, scale = 10; end
    colors = [0.85 0.10 0.10; 0.10 0.55 0.10; 0.05 0.05 0.05; 0.75 0 0.65];
    styles = {'-','--','-.',':'};
    handles = gobjects(nplot+2,1);
    count = 1;
    handles(count) = plot(ax,z,sp.Veff_meV,'Color',[0 0.3 0.8],'LineWidth',2, ...
        'DisplayName','$U_{\mathrm{eff}}(z)$');
    for i = 1:nplot
        density = abs(sp.Psi(:,i)).^2;
        if max(density) > 0, density = density/max(density); end
        yline(ax,E(i),':','Color',colors(i,:),'HandleVisibility','off');
        count = count+1;
        handles(count) = plot(ax,z,E(i)+scale*density,styles{i}, ...
            'Color',colors(i,:),'LineWidth',1.55, ...
            'DisplayName',sprintf('$n=%d$',i-1));
    end
    if isfinite(sp.EF_meV)
        count = count+1;
        handles(count) = yline(ax,sp.EF_meV,'--','Color',[0.85 0.55 0], ...
            'LineWidth',1.4,'DisplayName','$E_{\mathrm F}$');
    end
    xlabel(ax,'$z$ (nm)','Interpreter','latex');
    ylabel(ax,'Energy (meV)','Interpreter','latex');
    xlim(ax,[min(z) max(z)]);
    legend(ax,handles(1:count),'Location','best','Interpreter','latex');
    format_axes(ax);
end

function fig = plot_moap_spectrum(spectrum, cfg, mechanism)
    E = spectrum.energy_meV(:);
    P = spectrum.mechanism.(lower(mechanism)).total_plot(:);
    if strcmpi(mechanism,'optical'), color = [0.35 0 0.55];
    else, color = [0.75 0.1 0.1]; end
    fig = figure('Visible',cfg.output.visible,'Color','w','Units','pixels', ...
        'Position',[100 100 680 500]);
    ax = axes(fig);
    plot(ax,E,P,'Color',color,'LineWidth',1.75);
    format_spectrum_axes(ax,E,P,cfg);
end

function fig = plot_moap_contributions(spectrum, cfg, mechanism)
    E = spectrum.energy_meV(:);
    colors = [0 0.3 0.8; 0.85 0.1 0.1; 0.1 0.55 0.1; 0.65 0.1 0.65];
    fig = figure('Visible',cfg.output.visible,'Color','w','Units','pixels', ...
        'Position',[100 100 720 520]);
    ax = axes(fig); hold(ax,'on');
    max_curves = 2*numel(cfg.oap.photon_orders);
    handles = gobjects(max_curves,1);
    labels = cell(max_curves,1);
    all_curves = nan(numel(E),max_curves);
    count = 0;
    for j = 1:numel(cfg.oap.photon_orders)
        ell = cfg.oap.photon_orders(j);
        key = sprintf('order_%d',ell);
        color = colors(mod(j-1,size(colors,1))+1,:);
        if strcmpi(mechanism,'optical')
            curves = {spectrum.mechanism.optical.(key).emission_plot(:), ...
                spectrum.mechanism.optical.(key).absorption_plot(:)};
            styles = {'-','--'};
            names = {sprintf('$\\ell=%d$ emission',ell), ...
                sprintf('$\\ell=%d$ absorption',ell)};
        else
            curves = {spectrum.mechanism.piezoelectric.(key).total_plot(:)};
            styles = {'-'};
            names = {sprintf('$\\ell=%d$',ell)};
        end
        for k = 1:numel(curves)
            count = count+1;
            handles(count) = plot(ax,E,curves{k},styles{k}, ...
                'Color',color,'LineWidth',1.65);
            labels{count} = names{k};
            all_curves(:,count) = curves{k};
        end
    end
    legend(ax,handles(1:count),labels(1:count),'Location','best','Interpreter','latex');
    format_spectrum_axes(ax,E,all_curves(:,1:count),cfg);
end

function format_spectrum_axes(ax, E, P, cfg)
    xlabel(ax,'Photon energy $\hbar\Omega$ (meV)','Interpreter','latex');
    ylabel(ax,'MOAP (arb. units)','Interpreter','latex');
    valid = isfinite(E);
    if cfg.oap.plot_full_energy_range || ~any(isfinite(P) & P > 0)
        xlim(ax,[min(E(valid)) max(E(valid))]);
    else
        Pmax = max(P(isfinite(P)));
        active = E(isfinite(P) & P >= 1e-5*Pmax);
        xpad = 0.035*(max(active)-min(active));
        xlim(ax,[max(min(E(valid)),min(active)-xpad), ...
            min(max(E(valid)),max(active)+xpad)]);
    end
    finiteP = P(isfinite(P));
    if ~isempty(finiteP) && max(finiteP) > 0 && min(finiteP) >= -1e-10*max(finiteP)
        ylim(ax,[0 1.06*max(finiteP)]);
    end
    format_axes(ax);
end
