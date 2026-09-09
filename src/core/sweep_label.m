function label = sweep_label(parameter, value)
%SWEEP_LABEL Legend text for sweep plots.
    switch lower(parameter)
        case 'b_t'
            label = sprintf('B = %.3g T', value);
        case 't_k'
            label = sprintf('T = %.3g K', value);
        case 'lz_nm'
            label = sprintf('L_z = %.3g nm', value);
        case 'e_kvcm'
            label = sprintf('E = %.3g kV/cm', value);
        case 'nd_sheet_cm2'
            label = sprintf('N_D = %.3g cm^{-2}', value);
        case 'manning_prefactor'
            label = sprintf('\\nu = %.3g', value);
        case 'u0_mev'
            label = sprintf('U_0 = %.3g meV', value);
        otherwise
            label = sprintf('%s = %.3g', parameter, value);
    end
end
