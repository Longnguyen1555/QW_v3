function [C2V, hw, Nph] = phonon_coupling_density(mechanism, q, qd, cfg)
%PHONON_COUPLING_DENSITY Return V*|C_q|^2, phonon energy and Bose factor.
%
% Optical Eq. (2): q^2/(q^2+qd^2)^2 with constant LO energy.
% Piezo Eq. (18): q^3/(q^2+qd^2)^2 with acoustic hw=hbar*c_s*q.
% The caller supplies q=q_perp for the document approximation used by
% direct_q_integral.  The normalization volume cancels the q-state density.

    c = cfg.constants;
    epss = cfg.material.eps_static;

    denominator = (q.^2 + qd.^2).^2;
    denominator(denominator == 0) = Inf;

    switch lower(mechanism)
        case 'optical'
            hw0 = cfg.material.LO_phonon_meV * c.meV;
            dielectric = (1/cfg.material.eps_high - 1/epss);
            C2V = c.e^2*hw0/(2*c.eps0) .* dielectric .* ...
                  q.^2 ./ denominator;
            hw = hw0 .* ones(size(q));

        case {'piezoelectric','piezo'}
            s = cfg.material.sound_speed_mps;
            kappa2 = cfg.material.piezo_kappa2;
            C2V = kappa2*c.hbar*c.e^2*s/(2*epss*c.eps0) .* ...
                  q.^3 ./ denominator;
            hw = c.hbar*s.*q;

        otherwise
            error('Unsupported phonon mechanism: %s', mechanism);
    end

    if any(hw(:) < 0) || ~(isfinite(cfg.temperature_K) && cfg.temperature_K > 0)
        error('QW:PhononPopulation', ...
            'Phonon energies must be nonnegative and T must be positive.');
    end

    x = hw ./ (c.kB*cfg.temperature_K);
    Nph = zeros(size(hw));
    positive = (x > 0);
    Nph(positive) = 1 ./ expm1(x(positive));
    Nph(x == 0) = Inf;
end
