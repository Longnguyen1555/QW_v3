function case_result = run_single_case(cfg, tag)


    if nargin < 2
        tag = 'case';
    end

    if cfg.run.verbose

        fprintf('\n');
        fprintf('============================================================\n');
        fprintf('[MOAP] Running case: %s\n', tag);

        fprintf( ...
            '[MOAP] B = %.4g T, E = %.4g kV/cm, T = %.4g K\n', ...
            cfg.fields.B_T, ...
            cfg.fields.E_kVcm, ...
            cfg.temperature_K);

        fprintf( ...
            '[MOAP] Lz = %.4g nm, alpha = %.4g, Nd = %.4g cm^-2\n', ...
            cfg.structure.Lz_nm, ...
            cfg.structure.alpha, ...
            cfg.doping.Nd_sheet_cm2);

    end

    sp = solve_schrodinger_poisson(cfg);

    if ~sp.converged
        warning( ...
            'Schrodinger-Poisson solver reached max_iter without convergence.');
    end

    transition_data = compute_transition_data(sp, cfg);

    spectrum = compute_moap_spectrum(sp, transition_data, cfg);

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


    figures = { ...
        fig_structure, ...
        fig_spectrum_optical, ...
        fig_spectrum_piezoelectric, ...
        fig_contributions_optical, ...
        fig_contributions_piezoelectric ...
    };


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