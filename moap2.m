function results = moap2(cfg)
    close all;
    if nargin == 0 || isempty(cfg)
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
        cfg.structure.auto_domain = false;
        cfg.structure.domain_factor = 3;
        cfg.structure.Nz = 1502;
        cfg.structure.n_states = 6;
        cfg.doping.Nd_sheet_cm2 = 1e12;
        cfg.doping.width_nm = 2;
        cfg.fields.B_T = 5;
        cfg.fields.E_kVcm = 25;
        cfg.fields.kx_inv_m = 0;
        cfg.temperature_K = 100;
        cfg.sp.max_outer = 100;
        cfg.sp.mix = .4;
        cfg.sp.tol_VH_relative = 1e-3;
        cfg.sp.tol_EF_meV = 1e-3;
        cfg.laser.E0_kVcm = 4.5;
        cfg.laser.a0_nm = 7.5;
        cfg.oap.transitions = [1 2];
        cfg.oap.photon_orders = [1 2];
        cfg.oap.photon_energy_meV = 10:.01:150;
        cfg.oap.qz_max_inv_nm = 2;
        cfg.oap.Nqz = 101;
        cfg.oap.screening_density_mode = 'from_sheet';
        cfg.oap.fixed_ne_cm3 = 1e18;
        cfg.oap.active_thickness_nm = 10;
        cfg.oap.plot_normalization = 'none';
        cfg.oap.plot_full_energy_range = true;
        cfg.oap.series.s_max = 2;
        cfg.oap.series.v_max = 3;
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
        cfg.constants.meV = cfg.constants.e*1e-3;
        cfg.constants.nm = 1e-9;
    end
    if ~exist(cfg.output.directory, 'dir')
        mkdir(cfg.output.directory);
    end
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
    E = spectrum.energy_meV;
    metrics.total = profile_fwhm(E,spectrum.total_raw);
    for mechanism = {'optical','piezoelectric'}
        mk = mechanism{1};
        metrics.(mk).total = profile_fwhm(E,spectrum.mechanism.(mk).total_raw);
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            metrics.(mk).(ok).emission = profile_fwhm(E,spectrum.mechanism.(mk).(ok).emission_raw);
            metrics.(mk).(ok).absorption = profile_fwhm(E,spectrum.mechanism.(mk).(ok).absorption_raw);
            metrics.(mk).(ok).total = profile_fwhm(E,spectrum.mechanism.(mk).(ok).total_raw);
        end
    end
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
        mechanism = {}; order_col = []; process = {}; peak_meV = []; fwhm_meV = []; hwhm_meV = [];
        for mech = {'optical','piezoelectric'}
            mk = mech{1};
            for order = cfg.oap.photon_orders
                ok = sprintf('order_%d',order);
                for proc = {'emission','absorption'}
                    item = metrics.(mk).(ok).(proc{1});
                    mechanism{end+1,1} = mk; order_col(end+1,1) = order; process{end+1,1} = proc{1}; %#ok<AGROW>
                    peak_meV(end+1,1) = item.peak_x; fwhm_meV(end+1,1) = item.fwhm; hwhm_meV(end+1,1) = item.hwhm; %#ok<AGROW>
                end
            end
        end
        writetable(table(mechanism,order_col,process,peak_meV,fwhm_meV,hwhm_meV,'VariableNames',{'Mechanism','PhotonOrder','Process','Peak_meV','FWHM_meV','HWHM_meV'}),fullfile(cfg.output.directory,'single_metrics.csv'));
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
    results.config = cfg;
    results.single = case_result;
    if cfg.output.save_mat
        save(fullfile(cfg.output.directory, 'project_results.mat'), 'results', 'cfg', '-v7.3');
    end
    fprintf('\n[MOAP] Completed successfully.\n[MOAP] Results folder: %s\n',cfg.output.directory);
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
        td(t).Q_half_inv_m = sum(abs(Iqz).^2 .* wqz);
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
    if strcmpi(cfg.oap.screening_density_mode,'fixed')
        qd = sqrt(max(cfg.oap.fixed_ne_cm3*1e6,0)*c.e^2/(cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
    else
        thickness = cfg.oap.active_thickness_nm*c.nm;
        if thickness <= 0, thickness = max(cfg.structure.Lz_nm*c.nm,sp.dz_m); end
        qd = sqrt(max(sp.Nd_sheet_m2/thickness,0)*c.e^2/(cfg.material.eps_static*c.eps0*c.kB*cfg.temperature_K));
    end
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
        P0_LO = E0^4*e^4*me^2*(hw0/hbar)*sqrt(kappa0)/(256*pi^3*hbar^6*eps0*kappa_LO*(Eph/hbar));
        P0_PE = E0^4*e^4*me^2*cfg.material.piezo_kappa2*kBT*sqrt(kappa0)/(256*pi^3*hbar^7*kappa0*eps0*(Eph/hbar));
        for it = 1:numel(td)
            eps_alpha = sp.E_J(td(it).initial);
            DeltaE = td(it).deltaE_J;
            transition_factor = abs(td(it).dipole_m)^2 * td(it).Q_half_inv_m;
            for ell = orders
                order_key = sprintf('order_%d', ell);
                order_factor = a0^(ell) / (2^(2*ell) * factorial(ell)^2);
                Pplus = P0_LO * transition_factor * order_factor * (N0 + 1.0) * compute_T_taylor_bessel_pdf(eps_alpha, sp.EF_J, DeltaE + hw0 - ell*Eph, ell, qd, cfg);
                Pminus = P0_LO * transition_factor * order_factor * N0 * compute_T_taylor_bessel_pdf(eps_alpha, sp.EF_J, DeltaE - hw0 - ell*Eph, ell, qd, cfg);
                spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;
                spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;
                spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;
                spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + P0_PE*transition_factor*order_factor*compute_T_taylor_bessel_pdf(eps_alpha,sp.EF_J,DeltaE-ell*Eph,ell,qd,cfg);
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
    scale = 1;
    if strcmpi(cfg.oap.plot_normalization,'global_max')
        scale = max(spectrum.total_raw);
        if ~isfinite(scale) || scale <= 0, scale = 1; end
    end
    spectrum.mechanism.optical.total_plot = spectrum.mechanism.optical.total_raw/scale;
    spectrum.mechanism.piezoelectric.total_plot = spectrum.mechanism.piezoelectric.total_raw/scale;
    for ell = orders
        ok = sprintf('order_%d',ell);
        spectrum.mechanism.optical.(ok).emission_plot = spectrum.mechanism.optical.(ok).emission_raw/scale;
        spectrum.mechanism.optical.(ok).absorption_plot = spectrum.mechanism.optical.(ok).absorption_raw/scale;
        spectrum.mechanism.optical.(ok).total_plot = spectrum.mechanism.optical.(ok).total_raw/scale;
        spectrum.mechanism.piezoelectric.(ok).total_plot = spectrum.mechanism.piezoelectric.(ok).total_raw/scale;
    end
end
function Tval = compute_T_taylor_bessel_pdf(eps_alpha_J, EF_J, D_J, ell, qd, cfg)
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
            LambdaKsum = besselk(abs(s-ell), v*absD/(2.0*kBT), 0) + besselk(abs(s-ell+1), v*absD/(2.0*kBT), 0);
            % Eq. (32) multiplied by Lambda
            LambdaJ = (hbar^2 * qd^2 / (2.0 * me)) * (hbar^2 * qd^2 / (2.0 * me * absD))^(s-ell) * LambdaKsum;
            Tval = Tval + (-1)^(s+v-1) * (s+1) * exp(-v*(eps_alpha_J+D_J/2.0-EF_J)/kBT) * qd_power * LambdaJ;
        end
    end
end
function m = profile_fwhm(x, y)
    x = x(:).'; y = y(:).';
    m = struct('peak',NaN,'peak_x',NaN,'fwhm',NaN,'hwhm',NaN, 'left_x',NaN,'right_x',NaN,'left_hwhm',NaN, 'right_hwhm',NaN,'half_level',NaN,'resolved',false);
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
    if y(il+1) == y(il), m.left_x = (x(il)+x(il+1))/2;
    else, m.left_x = x(il)+(h-y(il))*(x(il+1)-x(il))/(y(il+1)-y(il)); end
    ir = ir-1;
    if y(ir+1) == y(ir), m.right_x = (x(ir)+x(ir+1))/2;
    else, m.right_x = x(ir)+(h-y(ir))*(x(ir+1)-x(ir))/(y(ir+1)-y(ir)); end
    m.left_hwhm = m.peak_x-m.left_x;
    m.right_hwhm = m.right_x-m.peak_x;
    m.fwhm = m.right_x-m.left_x;
    m.hwhm = m.fwhm/2;
    m.resolved = true;
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
    if strcmp(mk,'optical'), curve_color = [.35 0 .55];
    else, curve_color = [.75 .1 .1]; end
    E = spectrum.energy_meV(:);
    P = spectrum.mechanism.(mk).total_plot(:);
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
