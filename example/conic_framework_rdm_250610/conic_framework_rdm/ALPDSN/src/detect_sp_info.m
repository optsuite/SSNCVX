function sp_info = detect_sp_info(At, K, Lchol, system_opt, socp_formula, sys2_sparse_strategy, threshold, fid)
    % For system_opt = 2, we detect the dense and sparse rows and columns.
    % For system_opt = 0, we assume the matrix is sparse.

    if system_opt ~= 0 && system_opt ~= 2
        error('system_opt must be 0 or 2');
    end

    if nargin < 8
        fid = 1;
    end

    % Get the number of columns (m) from the first block
    m = size(At{1}, 2);

    % Initialize variables
    n = 0;
    nnz_per_column = zeros(1, m);

    % Calculate total number of rows (n) and non-zero elements per column
    for i = 1:length(At)
        n = n + size(At{i}, 1);
        if system_opt == 2
            nnz_per_column = nnz_per_column + sum(spones(At{i}), 1);
        else
            % do nothing
        end
    end

    % Calculate sparsity for each column
    sparsity = nnz_per_column / n;

    row_den = find(sparsity > 0.1);
    row_sp = find(sparsity <= 0.1);
    dense_propotion = length(row_den) / (length(row_den) + length(row_sp));
    if dense_propotion > 0.95
        row_den = 1:m;
        row_sp = zeros(1, 0);
    elseif dense_propotion < 0.01
        row_den = zeros(1, 0);
        row_sp = 1:m;
    end

    At_row_den = cellfun(@(x) full(x(:, row_den)), At, 'UniformOutput', false);
    At_row_sp = cellfun(@(x) sparse(x(:, row_sp)), At, 'UniformOutput', false);
    At_row_sp_merged = cell2mat(At_row_sp);

    % Detect dense & sparse columns of new At
    m_sp = length(row_sp);
    AHAt_mat22_At_new = sparse(0, m_sp);

    for p =1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'q')
            if socp_formula == 1
                nrows = size(At{p}, 1);
                ncols = length(cone.size);

                Prow_index = 1:nrows;
                Pcol_index = repelem(1:ncols, cone.size);
                Pvalue1 = rand(nrows, 1) + 1;
                Pvalue2 = -Pvalue1;
                
                index = 1;
                for j = 1:ncols
                    Pvalue1(j) = sqrt(2) / 2;
                    Pvalue2(j) = sqrt(2) / 2;
                    index = index + cone.size(j);
                end
                

                Q1 = sparse(Prow_index, Pcol_index, Pvalue1, nrows, ncols);
                Q2 = sparse(Prow_index, Pcol_index, Pvalue2, nrows, ncols);

                temp_sp1 = At_row_sp{p}' * Q1;
                temp_sp2 = At_row_sp{p}' * Q2;

                AHAt_mat22_At_new = [AHAt_mat22_At_new; temp_sp1'; temp_sp2'];

            elseif socp_formula == 2 || socp_formula == 3
                nrows = size(At{p}, 1);
                ncols = length(cone.size);

                Prow_index = 1:nrows;
                Pcol_index = repelem(1:ncols, cone.size);
                Pvalue1 = rand(nrows, 1) + 1;

                Q1 = sparse(Prow_index, Pcol_index, Pvalue1, nrows, ncols);

                temp_sp1 = At_row_sp{p}' * Q1;

                AHAt_mat22_At_new = [AHAt_mat22_At_new; temp_sp1'];
            end
        end

    end

    if isnumeric(Lchol)
        Rperm = speye(m);
    else
        Rperm = Lchol.Rperm;
    end

    AHAt_mat22_At = [At_row_sp_merged; AHAt_mat22_At_new; Rperm(:, row_sp)];

    if system_opt == 2
        nnz_per_row = sum(spones(AHAt_mat22_At), 2);
    else
        nnz_per_row = zeros(size(AHAt_mat22_At, 1), 1);
    end

    if (sys2_sparse_strategy == 0)
        sparsity = nnz_per_row / size(AHAt_mat22_At, 2);
        col_den = find(sparsity > threshold);
        col_sp = find(sparsity <= threshold);
    else
        col_den = find(nnz_per_row >= sys2_sparse_strategy);
        col_sp = find(nnz_per_row < sys2_sparse_strategy);
    end

    % augmented_lhs = [AHAt.mat22.At_sp' .* AHAt.mat22.diag_sp' * AHAt.mat22.At_sp, sparse(AHAt.mat22.At_den' .* C_vec'); sparse(C_vec .* AHAt.mat22.At_den), -spdiags(1./AHAt.mat22.diag_den .* C_vec .* C_vec, 0, n_den, n_den)];
    AHAt_mat22_At_den = AHAt_mat22_At(col_den, :);
    AHAt_mat22_At_sp = AHAt_mat22_At(col_sp, :);
    AHAt_mat22_diag_sp = rand(length(col_sp), 1) + 1;
    lhs_structure = [AHAt_mat22_At_sp' .* AHAt_mat22_diag_sp' * AHAt_mat22_At_sp, AHAt_mat22_At_den'; AHAt_mat22_At_den, speye(length(col_den))];
    lhs_structure = spones(lhs_structure);

    sp_info = struct('row_den', row_den, ...
                     'row_sp', row_sp, ...
                     'row_index', [row_den, row_sp], ...
                     'col_den', col_den, ...
                     'col_sp', col_sp, ...
                     'At_row_den', {At_row_den}, ...
                     'At_row_sp', {At_row_sp}, ...
                     'At_row_sp_merged', At_row_sp_merged, ...
                     'AHAt_mat22_At', AHAt_mat22_At, ...
                     'lhs_structure', lhs_structure);

    if system_opt == 2
        fprintf(fid, 'System_opt = 2\n');
        fprintf(fid, 'Number of dense rows: %d\n', length(row_den));
        fprintf(fid, 'Number of sparse rows: %d\n', length(row_sp));
        fprintf(fid, 'Number of dense columns (new At): %d\n', length(col_den));
        fprintf(fid, 'Number of sparse columns (new At): %d\n', length(col_sp));
        fprintf(fid, 'Number of non-zero elements in lhs_structure: %d\n', nnz(lhs_structure));
    elseif system_opt == 0
        fprintf(fid, 'System_opt = 0\n');
        fprintf(fid, 'Number of non-zero elements in lhs_structure: %d\n', nnz(lhs_structure));
    end
end