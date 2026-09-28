function out = solve_schrodinger_poisson(cfg)
%SOLVE_SCHRODINGER_POISSON Self-consistent Eqs. (1)-(6).

    c = cfg.constants;

    if cfg.structure.auto_domain
        domain_nm = max(cfg.structure.domain_nm, ...
            cfg.structure.domain_factor * cfg.structure.Lz_nm);
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

    % PDF audit: diagonalize exactly the number of subbands requested by
    % cfg.structure.n_states; the previous hard-coded value 12 made the
    % configured state count ineffective.
    n_states = cfg.structure.n_states;
    if isfield(cfg.sp, 'max_iter')
        max_iter = cfg.sp.max_iter;
    else
        % Backward-compatible use of the original default configuration.
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
        H = build_sp_hamiltonian(z_in, dz, Vconf(2:end-1), ...
                                 VH_old(2:end-1), cfg);
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

        % PDF loop: U_H^(m+1)=beta*U_H,new+(1-beta)*U_H^m, and convergence
        % is measured between consecutive mixed potentials.  A 1 meV floor
        % makes the first relative residual well-defined when U_H^0=0.
        VH_next = cfg.sp.mix .* VH_new + (1.0-cfg.sp.mix) .* VH_old;
        dVH_meV = max(abs(VH_next - VH_old)) / c.meV;
        VH_scale = max(max(abs(VH_old)), c.meV);
        dVH_relative = max(abs(VH_next - VH_old)) / VH_scale;

        if isfinite(EF_old) && isfinite(EF)
            dEF_meV = abs(EF - EF_old) / c.meV;
        end
    

        if dVH_relative <= tol_VH_relative && ...
                dEF_meV <= cfg.sp.tol_EF_meV
            converged = true;
            VH_old = VH_next;
            break;
        end
        VH_old = VH_next;
       
        EF_old = EF;
    end

    % One final diagonalization with the accepted potential.
    H = build_sp_hamiltonian(z_in, dz, Vconf(2:end-1), ...
                             VH_old(2:end-1), cfg);
    [Psi_in, E] = lowest_eigenpairs(H, n_states);
    Psi_in = normalize_wavefunctions(Psi_in, dz);
    Psi = zeros(numel(z), n_states);
    Psi(2:end-1, :) = Psi_in;
    if Nd_sheet > 0
        EF = solve_fermi_level(Nd_sheet, E, cfg);
        Ni = subband_sheet_populations(EF, E, cfg);
    else
        % PDF charge-neutral empty-well limit: no donor charge means no
        % electron sheet population and therefore no Hartree source.
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
