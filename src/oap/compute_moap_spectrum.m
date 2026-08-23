function spectrum = compute_moap_spectrum(sp, td, cfg)
%COMPUTE_MOAP_SPECTRUM Select direct integral or analytical series.

    switch lower(cfg.oap.model)
        case 'direct_q_integral'
            spectrum = compute_moap_direct(sp, td, cfg);
        case 'analytical_series'
            spectrum = compute_moap_analytical_series(sp, td, cfg);
        case 'direct_t'
            spectrum = compute_moap_pdf_direct(sp, td, cfg);
        case 'taylor_bessel'
            spectrum = compute_moap_pdf_taylor_bessel(sp, td, cfg);            
        otherwise
            error('Unsupported OAP model: %s', cfg.oap.model);
    end
end
