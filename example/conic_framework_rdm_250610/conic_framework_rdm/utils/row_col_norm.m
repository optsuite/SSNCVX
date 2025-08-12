
% compute the norm of each row and each column of a sparse matrix
% [row_norm, col_norm] = row_col_norm(A, rol_norm_opt, col_norm_opt)
% Input
%     A: sparse or dense matrix of size m x n
%     row_norm_opt: 1 for 1-norm, 2 for 2-norm, 3 for inf-norm
%     col_norm_opt: 1 for 1-norm, 2 for 2-norm, 3 for inf-norm

% Output
%     row_norm: norm of each row of A, of size m x 1
%     col_norm: norm of each column of A, of size n x 1



function [row_norm, col_norm] = row_col_norm(A, row_norm_opt, col_norm_opt)
    if issparse(A)
        [row_norm, col_norm] = mexrow_col_norm_sparse(A, row_norm_opt, col_norm_opt);
    else
        if row_norm_opt == 1
            row_norm = vecnorm(A, 1, 2);
        elseif row_norm_opt == 2
            row_norm = vecnorm(A, 2, 2);
        elseif row_norm_opt == 3
            row_norm = vecnorm(A, inf, 2);
        end
        row_norm = reshape(row_norm, [], 1);

        if col_norm_opt == 1
            col_norm = vecnorm(A, 1, 1);
        elseif col_norm_opt == 2
            col_norm = vecnorm(A, 2, 1);
        elseif col_norm_opt == 3
            col_norm = vecnorm(A, inf, 1);
        end
        col_norm = reshape(col_norm, [], 1);
    end

end