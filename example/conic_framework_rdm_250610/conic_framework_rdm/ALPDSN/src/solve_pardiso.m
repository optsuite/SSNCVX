function X = solve_pardiso(A, B, pardiso_info)
    verbose = false;
    
    % Analyze the matrix and compute a symbolic factorization.
    info = pardisoinit(-2, 0);
    info = pardisoreorder(tril(A), info, verbose);
    % fprintf('The factors have %d nonzero entries.\n', info.iparm(18));

    % Compute the numeric factorization.
    info = pardisofactor(tril(A), info, verbose);
    % fprintf('The matrix has %d positive and %d negative eigenvalues.\n', ...
    %     info.iparm(22), info.iparm(23));

    % Compute the solutions X using the symbolic factorization.
    [X, info] = pardisosolve(tril(A), B, info, verbose);
    % fprintf('PARDISO performed %d iterative refinement steps.\n', info.iparm(7));

    % Free the Pardiso data structures.
    pardisofree(info);
    clear info;

end