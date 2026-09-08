function cfg = set_sweep_parameter(cfg, parameter, value)
%SET_SWEEP_PARAMETER Apply one supported sweep parameter.

    switch lower(parameter)
        case 'b_t'
            cfg.fields.B_T = value;
        case 't_k'
            cfg.temperature_K = value;
        case 'lz_nm'
            cfg.structure.Lz_nm = value;
        case 'e_kvcm'
            cfg.fields.E_kVcm = value;
        case 'nd_sheet_cm2'
            cfg.doping.Nd_sheet_cm2 = value;
        case 'alpha'
            error('QW:UnsupportedAlphaSweep', [ ...
                'An alpha sweep is unsupported for the Manning/sech ', ...
                'potential and would not change its physics. Sweep ', ...
                'manning_prefactor explicitly instead.']);
        case 'manning_prefactor'
            cfg.structure.manning_prefactor = value;
        case 'u0_mev'
            cfg.structure.U0_meV = value;
        otherwise
            error('Unsupported sweep parameter: %s', parameter);
    end
end
