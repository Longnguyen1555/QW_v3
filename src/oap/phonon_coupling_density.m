function [C2V, hw, Nph] = phonon_coupling_density(mechanism, q, qd, cfg)
%PHONON_COUPLING_DENSITY Return V*|C_q|^2, phonon energy and Bose factor.
%
% Optical: target-paper Eq. (21)
% Piezoelectric: target-paper Eq. (13)
% The normalization volume cancels the q-state density in the continuum sum.

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

    Nph = zeros(size(hw));
    mask = hw > 0;
    x = hw(mask)/(c.kB*cfg.temperature_K);
    Nph(mask) = 1 ./ expm1(x);
    Nph(~isfinite(Nph)) = 0;
end
