function spectrum = compute_moap_spectrum(sp, td, cfg)
%COMPUTE_MOAP_SPECTRUM Select an OAP kernel.

    switch lower(cfg.oap.model)
        % Preserve the pre-rename public model names while routing them to
        % the current source files; the physics kernels are not conflated.
        case 'lorentzian'
            spectrum = compute_moap_lorentzian(sp, td, cfg);
        case 'direct_exact'
            spectrum = compute_moap_direct(sp, td, cfg);
        case 'direct_approx'
            spectrum = compute_moap_direct_approx(sp, td, cfg);
        case {'van','analytical_series','analytic_series'}
            spectrum = compute_moap_van(sp, td, cfg);
        case 'nhan'
            spectrum = compute_moap_nhan(sp, td, cfg);
        otherwise
            error('Unsupported OAP model: %s', cfg.oap.model);
    end
end
