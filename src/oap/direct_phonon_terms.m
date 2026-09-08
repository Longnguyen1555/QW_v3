function terms = direct_phonon_terms(mechanism, qperp, qd, cfg)
%DIRECT_PHONON_TERMS Coupling times Bose factors for direct_q_integral.
%
% q in the coupling is intentionally q_perp, not sqrt(q_perp^2+q_z^2).
% q_z is handled separately through the confinement form factor.

    qperp = qperp(:);
    [C2V, hw, Nph] = phonon_coupling_density( ...
        mechanism, qperp, qd, cfg);

    switch lower(mechanism)
        case 'optical'
            % LO dispersion: hw_q=hw0; emission N0+1, absorption N0.
            absorption_C2V = C2V .* Nph;
            emission_C2V = C2V .* (Nph + 1.0);

        case {'piezoelectric', 'piezo'}
            % Source-document correction: its written Eq. (18) reuses
            % optical N0/omega0 despite introducing sound speed c_s.  The
            % physical piezoelectric acoustic model instead uses
            % omega_q=c_s*q_perp and N_q=1/expm1(hbar*c_s*q_perp/kBT).
            %
            % At q_perp=0, N_q diverges.  We take the limit of the complete
            % screened product C2V*N_q: for qd>0 it scales as q_perp^2 and
            % is zero.  We do not replace the Bose population itself by 0.
            absorption_C2V = zeros(size(C2V));
            positive_q = (qperp > 0);
            absorption_C2V(positive_q) = ...
                C2V(positive_q) .* Nph(positive_q);
            emission_C2V = absorption_C2V + C2V;

        otherwise
            error('QW:DirectMechanism', ...
                'Unsupported phonon mechanism: %s', mechanism);
    end

    terms = struct();
    terms.C2V = C2V;
    terms.hw_J = hw;
    terms.Nph = Nph;
    terms.emission_C2V = emission_C2V;
    terms.absorption_C2V = absorption_C2V;
end
