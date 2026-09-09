function J2 = landau_form_factor_squared(qperp, cfg)
%LANDAU_FORM_FACTOR_SQUARED Optional legacy in-plane factor.
%
% IMPORTANT FOR THE CURRENT QW MODEL:
% The physical magnetic field is in-plane and is already included in the
% z-dependent Hamiltonian through
%
%   (hbar*kx + e*B*z)^2/(2*m*)
%
% Therefore the production setting must normally be
%
%   cfg.oap.inplane_factor = 'none'
%
% The Landau form factor is kept only for legacy/comparison calculations.

    c = cfg.constants;

    mode = lower(strtrim(cfg.oap.inplane_factor));

    % ---------------------------------------------------------------------
    % Current physical model: no additional in-plane factor.
    % This check MUST occur before the B=0 fallback.
    % ---------------------------------------------------------------------
    if strcmp(mode, 'none')
        J2 = ones(size(qperp));
        return;
    end

    % ---------------------------------------------------------------------
    % Optional thermal recoil cutoff.
    % ---------------------------------------------------------------------
    if strcmp(mode, 'thermal')
        mstar = cfg.material.mstar_rel * c.m0;

        exponent = ...
            -cfg.oap.thermal_cutoff_factor .* ...
            c.hbar^2 .* qperp.^2 ./ ...
            (2*mstar*c.kB*cfg.temperature_K);

        J2 = exp(exponent);
        return;
    end

    % ---------------------------------------------------------------------
    % Legacy Landau model.
    % ---------------------------------------------------------------------
    if ~strcmp(mode, 'landau')
        error('Unknown cfg.oap.inplane_factor: %s', ...
              cfg.oap.inplane_factor);
    end

    B = cfg.fields.B_T;

    if B <= 1e-12
        % A Landau length is undefined at B=0.
        % Use the old thermal fallback only in explicit legacy mode.
        mstar = cfg.material.mstar_rel * c.m0;

        exponent = ...
            -cfg.oap.thermal_cutoff_factor .* ...
            c.hbar^2 .* qperp.^2 ./ ...
            (2*mstar*c.kB*cfg.temperature_K);

        J2 = exp(exponent);
        return;
    end

    Ni = cfg.oap.landau_initial;
    Nf = cfg.oap.landau_final;

    lB = sqrt(c.hbar/(c.e*B));
    u = 0.5*lB^2 .* qperp.^2;

    n1 = max(Ni, Nf);
    n2 = min(Ni, Nf);
    dN = n1 - n2;

    L = generalized_laguerre_numeric(n2, dN, u);

    log_pref = gammaln(n2+1) - gammaln(n1+1);

    J2 = exp(log_pref-u) .* ...
         (u.^dN) .* ...
         (L.^2);
end