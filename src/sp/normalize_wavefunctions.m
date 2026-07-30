function Psi = normalize_wavefunctions(Psi, dz)
%NORMALIZE_WAVEFUNCTIONS Enforce integral |psi_i|^2 dz = 1.
    norms = sqrt(sum(abs(Psi).^2, 1) .* dz);
    if any(norms <= 0 | ~isfinite(norms))
        error('Invalid wavefunction norm.');
    end
    Psi = bsxfun(@rdivide, Psi, norms);
end
