function [Psi, E] = lowest_eigenpairs(H, n_states)
%LOWEST_EIGENPAIRS Robust lowest eigenpairs for a symmetric sparse matrix.

    H = (H + H') / 2;
    opts.issym = true;
    opts.isreal = true;
    opts.tol = 1e-10;
    opts.maxit = 4000;

    try
        [Psi, D] = eigs(H, n_states, 'smallestreal', opts);
    catch
        try
            [Psi, D] = eigs(H, n_states, 'sa', opts);
        catch
            [Vfull, Dfull] = eig(full(H));
            evals = real(diag(Dfull));
            [evals, idx] = sort(evals, 'ascend');
            Psi = real(Vfull(:, idx(1:n_states)));
            D = diag(evals(1:n_states));
        end
    end

    E = real(diag(D));
    [E, idx] = sort(E, 'ascend');
    Psi = real(Psi(:, idx));

    if any(~isfinite(E))
        error('Non-finite eigenvalues returned by the eigensolver.');
    end
end
