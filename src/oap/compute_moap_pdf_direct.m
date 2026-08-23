function spectrum = compute_moap_pdf_direct(sp, td, cfg)
%COMPUTE_MOAP_PDF_DIRECT
% direct numerical evaluation of Eq. (20)/(28).


    c = cfg.constants;

    hbar = c.hbar;
    e    = c.e;
    kBT  = c.kB * cfg.temperature_K;
    eps0 = c.eps0;
    cfg.oap.normalization_area_m2 = 1.0;
    me = cfg.material.mstar_rel * c.m0;

    kappa0   = cfg.material.eps_static;
    kappainf = cfg.material.eps_high;

    % Eq. (7):
   
    kappa_LO = 1.0 / (1.0/kappainf - 1.0/kappa0);

    hw0 = cfg.material.LO_phonon_meV * c.meV;
    omega0 = hw0 / hbar;

    % Bose population of LO phonon
    N0 = 1.0 / expm1(hw0 / kBT);

    % Strong EM field
    E0 = cfg.laser.E0_kVcm * 1.0e5;

    % Normalization area S
    S = cfg.oap.normalization_area_m2;

    if ~(isfinite(S) && S > 0)
        error(['cfg.oap.normalization_area_m2 must be explicitly ', ...
               'specified and positive.']);
    end

    orders = cfg.oap.photon_orders(:).';

    if ~isequal(orders, [1 2])
        error(['Tinh_toan_chi_tiet.pdf derives the final model ', ...
               'for photon orders l = [1 2] only.']);
    end

    Ephot = cfg.oap.photon_energy_meV(:).' * c.meV;
    nE = numel(Ephot);

    [qd, ne3d] = debye_wavevector(sp, cfg);

    spectrum = initialize_spectrum_struct(Ephot, cfg);

    spectrum.meta.model = 'direct_t';
    spectrum.meta.qd_inv_m = qd;
    spectrum.meta.qd_inv_nm = qd * c.nm;
    spectrum.meta.screening_ne_m3 = ne3d;
    spectrum.meta.normalization_area_m2 = S;

    % Keep these arrays for later audit.
    spectrum.meta.a0_m = zeros(1, nE);
    spectrum.meta.Omega_rad_s = zeros(1, nE);
    spectrum.meta.P0_optical = zeros(1, nE);
    spectrum.meta.P0_piezoelectric = zeros(1, nE);

    for iE = 1:nE

        Eph = Ephot(iE);

        % hbar Omega = photon energy
        Omega = Eph / hbar;

        a0 = e * E0 / (me * Omega^2);

        spectrum.meta.Omega_rad_s(iE) = Omega;
        spectrum.meta.a0_m(iE) = a0;

        % ================================================================
        % Optical prefactor -- Eq. (22), PDF-exact convention
        % ================================================================
    
        P0_LO = E0^4 * e^4 * me^2 * omega0 * S * sqrt(kappa0)/ (256.0 * pi^3 * hbar^6 * eps0 * kappa_LO * Omega^2);

        % ================================================================
        % Piezoelectric prefactor -- Eq. (27)
        % ================================================================
        
        P0_PE = E0^4 * e^4 * me^2 * cfg.material.piezo_kappa2 * kBT * S * sqrt(kappa0) / (256.0 * pi^3 * hbar^7 * kappa0 * eps0 * Omega^2);

        spectrum.meta.P0_optical(iE) = P0_LO;
        spectrum.meta.P0_piezoelectric(iE) = P0_PE;

        for it = 1:numel(td)

            alpha = td(it).initial;

            eps_alpha = sp.E_J(alpha);
            EF = sp.EF_J;

            DeltaE = td(it).deltaE_J;

            % |Q_ab|^2 K_ab
           
            if isfield(td(it), 'Q2_ab_m2')
                Q2 = td(it).Q2_ab_m2;
            else
                Q2 = abs(td(it).dipole_m)^2;
            end

            if isfield(td(it), 'K_ab_inv_m')
                Kab = td(it).K_ab_inv_m;
            else
                Kab = td(it).Q_half_inv_m;
            end

            transition_factor = Q2 * Kab;

            transition_total_LO = 0.0;
            transition_total_PE = 0.0;

            for ell = orders

                order_key = sprintf('order_%d', ell);

                % Exactly:
                %
                % ell = 1 -> a0^2 / 4
                % ell = 2 -> a0^4 / 64
                %
                order_factor = a0^(2*ell) / (2^(2*ell) * factorial(ell)^2);

                % ========================================================
                % 1. ELECTRON - LO PHONON
                % ========================================================
               
                Dplus = DeltaE + hw0 - ell * Eph;
                Dminus = DeltaE - hw0 - ell * Eph;

                Tplus = compute_T_direct_pdf(eps_alpha, EF, Dplus, ell, qd, cfg);

                Tminus = compute_T_direct_pdf(eps_alpha, EF, Dminus, ell, qd, cfg);

                Pplus = P0_LO * transition_factor * order_factor * (N0 + 1.0) * Tplus;

                Pminus = P0_LO * transition_factor * order_factor * N0 * Tminus;

                spectrum.mechanism.optical.(order_key).emission_raw(iE) = spectrum.mechanism.optical.(order_key).emission_raw(iE) + Pplus;

                spectrum.mechanism.optical.(order_key).absorption_raw(iE) = spectrum.mechanism.optical.(order_key).absorption_raw(iE) + Pminus;

                spectrum.mechanism.optical.(order_key).total_raw(iE) = spectrum.mechanism.optical.(order_key).total_raw(iE) + Pplus + Pminus;

                transition_total_LO = transition_total_LO + Pplus + Pminus;

                % ========================================================
                % 2. ELECTRON - PIEZOELECTRIC PHONON
                % ========================================================
                
                Dpiezo = DeltaE - ell * Eph;

                Tpiezo = compute_T_direct_pdf(eps_alpha, EF, Dpiezo, ell, qd, cfg);

                Ppiezo = P0_PE * transition_factor * order_factor * Tpiezo;

                spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) = spectrum.mechanism.piezoelectric.(order_key).total_raw(iE) + Ppiezo;

                transition_total_PE = transition_total_PE + Ppiezo;
            end

            tr_key = sprintf('transition_%d_%d', ...
                             td(it).initial, td(it).final);

            if ~isfield(spectrum.mechanism.optical.transitions, tr_key)
                spectrum.mechanism.optical.transitions.(tr_key) = zeros(1, nE);

                spectrum.mechanism.piezoelectric.transitions.(tr_key) = zeros(1, nE);
            end

            spectrum.mechanism.optical.transitions.(tr_key)(iE) = transition_total_LO;

            spectrum.mechanism.piezoelectric.transitions.(tr_key)(iE) = transition_total_PE;
        end
    end

    % ====================================================================
    % Sum mechanisms
    % ====================================================================

    for ell = orders
        order_key = sprintf('order_%d', ell);

        spectrum.mechanism.optical.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.optical.(order_key).total_raw;

        spectrum.mechanism.piezoelectric.total_raw = spectrum.mechanism.piezoelectric.total_raw + spectrum.mechanism.piezoelectric.(order_key).total_raw;
    end

    spectrum.total_raw = spectrum.mechanism.optical.total_raw + spectrum.mechanism.piezoelectric.total_raw;

    % plot_normalization must be 'none'.
    
    spectrum = normalize_spectrum_for_plot(spectrum, cfg);
end

function Tval = compute_T_direct_pdf(eps_alpha_J, EF_J, D_J, ell, qd, cfg)

    c = cfg.constants;

    hbar = c.hbar;
    kBT  = c.kB * cfg.temperature_K;
    me   = cfg.material.mstar_rel * c.m0;

    if ~(isscalar(ell) && any(ell == [1 2]))
        error('ell must be 1 or 2.');
    end

    if ~(isfinite(qd) && qd > 0)
        error('Debye wavevector qd must be finite and positive.');
    end

    integrand = @(q) local_integrand( q, eps_alpha_J, EF_J, D_J, ell, qd, me, hbar, kBT);

    Tval = integral( integrand, 0.0, Inf, 'RelTol', cfg.oap.direct.RelTol, 'AbsTol', cfg.oap.direct.AbsTol);
end


function y = local_integrand(q, eps_alpha, EF, D, ell, qd, me, hbar, kBT)

    y = zeros(size(q));

    mask = (q > 0);

    qq = q(mask);

    eps_star = eps_alpha + hbar^2 .* qq.^2 ./ (8.0 * me) + D ./ 2.0 + me .* D.^2 ./ (2.0 * hbar^2 .* qq.^2);

    x = (eps_star - EF) ./ kBT;

    F = stable_fermi(x);

    screening = qq.^(2*ell + 1) ./ (qq.^2 + qd^2).^2;

    Babs = hbar^2 .* qq.^2 ./ (2.0 * me) + abs(D);

    y(mask) = screening .* Babs .* F;
end


function F = stable_fermi(x)

    F = zeros(size(x));

    pos = (x >= 0);

    t = exp(-x(pos));
    F(pos) = t ./ (1.0 + t);

    t = exp(x(~pos));
    F(~pos) = 1.0 ./ (1.0 + t);
end