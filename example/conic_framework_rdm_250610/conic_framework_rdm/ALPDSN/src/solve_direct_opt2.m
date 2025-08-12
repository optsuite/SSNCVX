function [dy, relres, flag] = solve_direct_opt2(K, At, rhs1, iHWy, Lchol, trans)
    % lhs* dy = rhs1
    % lhs = iHWy.sig * Rt^{-1} * S * A * iHWy.Dsch * At * Sinv * R^{-1} + epsilon * I
    % Sinv * y is equivalent to y(perm) = y
    % S * y is equivalent to y = y(perm)
    % lhs * dy = rhs is equivalent to 
    % (iHWy.sig  * A * iHWy.Dsch * At  + epsilon * Sinv * Rt * R * S) * (Sinv * R^{-1} * dy) = Sinv * Rt * rhs
    % compute RtR = Sinv * Rt * R * S 

    sp_info = trans.sp_info;
    direct_solver = trans.direct_solver;

    % unzip sp_info
    row_den = sp_info.row_den;
    row_sp = sp_info.row_sp;
    row_index = sp_info.row_index;
    m_den = length(row_den);
    m_sp = length(row_sp);
    At_row_den = sp_info.At_row_den;
    At_row_sp = sp_info.At_row_sp;
    At_row_sp_merged = sp_info.At_row_sp_merged;

    col_den = sp_info.col_den;
    col_sp = sp_info.col_sp;
    n_den = length(col_den);
    n_sp = length(col_sp);

    % count number of columns of A, and number of cones
    n_col = 0;
    num_cone = 0;
    for p = 1:length(K)
        if strcmp(K{p}.type, 'q')
            num_cone = num_cone + length(K{p}.size);
        end
        n_col = n_col + size(At{p}, 1);
    end

    if trans.socp_formula == 1
        n_total = n_col + 2 * num_cone + m_den + m_sp;
    elseif trans.socp_formula == 2 || trans.socp_formula == 3
        n_total = n_col + num_cone + m_den + m_sp;
    end

    n_count_new = 0;
    n_count_diag = 0;

    % fprintf("total number of columns of A: %d, number of cones: %d\n", n_total, num_cone);

    % compute lhs = iHWy.sig * A * iHWy.Dsch * At  + epsilon * Sinv * Rt * R * S
    par = iHWy; 
    AHAt = struct();
    AHAt.mat11 = zeros(m_den, m_den);
    AHAt.mat21 = zeros(m_sp, m_den);
    AHAt.mat22 = struct();
    AHAt.mat22.At_new = sparse(0, m_sp);
    AHAt.mat22.diag = zeros(n_total, 1);

    for p =1: length(K)
        cone = K{p};
        n_col_block = size(At{p}, 1);
        if strcmp(cone.type, 's')
            % do nothing
        elseif strcmp(cone.type, 'q')
            if trans.socp_formula == 1
                num_sub_block = length(cone.size);

                Q1 = par.Q1{p};
                Q2 = par.Q2{p};
            
                % temp1_den = At_row_den{p}' * Q1;
                % temp1_sp = At_row_sp{p}' * Q1;
                % temp2_den = At_row_den{p}' * Q2;
                % temp2_sp = At_row_sp{p}' * Q2;

                temp1_den_t = Q1' * At_row_den{p};
                temp1_sp_t = Q1' * At_row_sp{p};
                temp2_den_t = Q2' * At_row_den{p};
                temp2_sp_t = Q2' * At_row_sp{p};

                mat11_1 = temp1_den_t' .* par.Dsch1{p}' * temp1_den_t;
                mat11_2 = temp2_den_t' .* par.Dsch2{p}' * temp2_den_t;
                mat11_3 = At_row_den{p}' .* (repelem(par.shift{p}, cone.size, 1))' * At_row_den{p};

                mat21_1 = temp1_sp_t' .* par.Dsch1{p}' * temp1_den_t;
                mat21_2 = temp2_sp_t' .* par.Dsch2{p}' * temp2_den_t;
                mat21_3 = At_row_sp{p}' .* (repelem(par.shift{p}, cone.size, 1))' * At_row_den{p};

                AHAt.mat11 = AHAt.mat11 + mat11_1 + mat11_2 + mat11_3;
                AHAt.mat21 = AHAt.mat21 + mat21_1 + mat21_2 + mat21_3;

                AHAt.mat22.At_new = [AHAt.mat22.At_new; temp1_sp_t; temp2_sp_t];

                AHAt.mat22.diag(n_count_diag + 1: n_count_diag + n_col_block) = repelem(par.shift{p}, cone.size, 1);
                AHAt.mat22.diag(n_count_new + n_col + 1: n_count_new + n_col + num_sub_block) = par.Dsch1{p};
                AHAt.mat22.diag(n_count_new + n_col + num_sub_block + 1: n_count_new + n_col + 2 * num_sub_block) = par.Dsch2{p};

                n_count_new = n_count_new + 2 * num_sub_block;
                n_count_diag = n_count_diag + n_col_block;

            elseif trans.socp_formula == 2 || trans.socp_formula == 3
                num_sub_block = length(cone.size);

                Q3 = par.Q3{p};

                temp1_den_t = Q3' * At_row_den{p};
                temp1_sp_t = Q3' * At_row_sp{p};

                mat11_1 = temp1_den_t' .* par.Dsch3{p}' * temp1_den_t;
                mat11_3 = At_row_den{p}' .* par.shift3{p}' * At_row_den{p};

                mat21_1 = temp1_sp_t' .* par.Dsch3{p}' * temp1_den_t;
                mat21_3 = At_row_sp{p}' .* par.shift3{p}' * At_row_den{p};

                AHAt.mat11 = AHAt.mat11 + mat11_1 + mat11_3;
                AHAt.mat21 = AHAt.mat21 + mat21_1 + mat21_3;

                AHAt.mat22.At_new = [AHAt.mat22.At_new; temp1_sp_t];

                AHAt.mat22.diag(n_count_diag + 1: n_count_diag + n_col_block) = par.shift3{p};
                AHAt.mat22.diag(n_count_new + n_col + 1: n_count_new + n_col + num_sub_block) = par.Dsch3{p};

                n_count_new = n_count_new + num_sub_block;
                n_count_diag = n_count_diag + n_col_block;
            end

        elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'b') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
            if (~isempty(par.Dsch12{p}))
                diag_den = par.Dsch12{p}';
                diag_sp = sparse(diag_den);

                mat11_1 = At_row_den{p}' .* diag_den * At_row_den{p};
                mat21_1 = At_row_sp{p}' .* diag_sp * At_row_den{p};

                AHAt.mat11 = AHAt.mat11 + mat11_1;
                AHAt.mat21 = AHAt.mat21 + mat21_1;

                AHAt.mat22.diag(n_count_diag + 1: n_count_diag + n_col_block) = diag_sp;

                n_count_diag = n_count_diag + n_col_block;
            end
        elseif strcmp(cone.type, 'u')
        
        end
    end

    % AHAt.mat11 = par.sig * AHAt.mat11;
    % AHAt.mat12 = par.sig * AHAt.mat12;
    % AHAt.mat22 = par.sig * AHAt.mat22;

    offset = par.epsilon * Lchol.RtR;
    % speye_m = speye(m_den + m_sp);
    % offset = par.epsilon * speye_m;
    AHAt.mat11 = AHAt.mat11 + offset(row_den, row_den);
    AHAt.mat21 = AHAt.mat21 + offset(row_sp, row_den);

    n_col_block = m_den + m_sp;
    AHAt.mat22.At = [At_row_sp_merged; AHAt.mat22.At_new; Lchol.Rperm(:, row_sp)];
    % AHAt.mat22.At = [At_row_sp_merged; AHAt.mat22.At_new; speye_m(:, row_sp)];
    AHAt.mat22.diag(n_total - n_col_block + 1: n_total) = par.epsilon * ones(n_col_block, 1);
    % small_eps = 1e-10;
    % AHAt.mat22.diag(abs(AHAt.mat22.diag) < small_eps) = rand(sum(abs(AHAt.mat22.diag) < small_eps), 1) * small_eps;

    % print n_den and n_sp
    % fprintf('Number of dense columns: %d, number of sparse columns: %d\n', n_den, n_sp);

    AHAt.mat22.At_den = AHAt.mat22.At(col_den, :);
    AHAt.mat22.At_sp = AHAt.mat22.At(col_sp, :);
    AHAt.mat22.diag_den = AHAt.mat22.diag(col_den);
    AHAt.mat22.diag_sp = AHAt.mat22.diag(col_sp);

    % AHAt.mat22.mat = AHAt.mat22.At_den' .* AHAt.mat22.diag_den' * AHAt.mat22.At_den + ...
    %     AHAt.mat22.At_sp' .* AHAt.mat22.diag_sp' * AHAt.mat22.At_sp;
    
    % compute rhs1_temp = Sinv * Rt * rhs1
    rhs1_temp = zeros(size(rhs1));
    rhs1_temp(Lchol.perm) = Lchol.Rt * rhs1;
    % rhs1_temp = rhs1;

    % solve lhs * ( Sinv * R^{-1} * dy) = rhs
    % dy = AHAt.mat \ rhs1_temp;
    rhs1_sp = rhs1_temp(row_sp, :);
    rhs1_den = rhs1_temp(row_den, :);
    mat21_tilde = [AHAt.mat21, rhs1_sp];

    % mat22_inv_mat21_tilde = AHAt.mat22.mat \ mat21_tilde;
    C_vec = sqrt(abs(AHAt.mat22.diag_den));
    augmented_lhs = [AHAt.mat22.At_sp' .* AHAt.mat22.diag_sp' * AHAt.mat22.At_sp, sparse(AHAt.mat22.At_den' .* C_vec'); sparse(C_vec .* AHAt.mat22.At_den), -spdiags(1./AHAt.mat22.diag_den .* C_vec .* C_vec, 0, n_den, n_den)];

    % %% print nnz of augmented_lhs
    % fprintf('nnz of augmented_lhs: %d, %d, %d. Sum: %d.\n', nnz(AHAt.mat22.At_sp' .* AHAt.mat22.diag_sp' * AHAt.mat22.At_sp), nnz(sparse(AHAt.mat22.At_den' .* C_vec')), nnz(-spdiags(1./AHAt.mat22.diag_den .* C_vec .* C_vec, 0, n_den, n_den)), nnz(augmented_lhs));

    % %% for testing
    % detect_zeroA_nonzeroB(AHAt);

    augmented_rhs = [mat21_tilde; zeros(n_den, m_den + 1)];

    if strcmp(direct_solver, 'backslash')
        augmented_x = augmented_lhs \ augmented_rhs;
    elseif strcmp(direct_solver, 'ldl')
        [L_ldl, D_ldl, P_ldl] = ldl(augmented_lhs, 'vector'); % 'vector' for permutation vector
        % d = diag(D_ldl);
        % tolerance = 1e-10;
        % isPSD = all(d >= -tolerance);
        % if isPSD
        %     fprintf('LHS is positive semi-definite.\n');
        % else
        %     fprintf('LHS is NOT positive semi-definite. The smallest pseudo eigenvalue is: %e\n', full(min(d)));
        % end
        bp_ldl = augmented_rhs(P_ldl, :);
        y_ldl = L_ldl \ bp_ldl;
        z_ldl = D_ldl \ y_ldl;
        x_permuted = L_ldl' \ z_ldl;
        augmented_x = zeros(size(augmented_rhs, 1), size(augmented_rhs, 2));
        augmented_x(P_ldl, :) = x_permuted;
    elseif strcmp(direct_solver, 'pardiso')
        augmented_x = solve_pardiso(augmented_lhs, augmented_rhs, trans.pardiso_info);
    end

    % fprintf('norm of residual: %e\n', norm(augmented_lhs * augmented_x - augmented_rhs, 'fro'));
    mat22_inv_mat21_tilde = augmented_x(1:m_sp, :);

    mat22_inv_mat21 = mat22_inv_mat21_tilde(:, 1:end-1);
    mat22_inv_rhs1_sp = mat22_inv_mat21_tilde(:, end);
    reduced_lhs = AHAt.mat11 - AHAt.mat21' * mat22_inv_mat21;
    reduced_rhs = rhs1_den - mat22_inv_mat21' * rhs1_sp;
    dy_den = reduced_lhs \ reduced_rhs;

    % dy_sp = AHAt.mat22.mat \ (rhs1_sp - AHAt.mat21 * dy_den);
    dy_sp = mat22_inv_rhs1_sp - mat22_inv_mat21 * dy_den;
    dy = [dy_den; dy_sp];
    dy(row_index) = dy;

    % recover dy 
    dy = Lchol.R * dy(Lchol.perm);

    % check residual
    relres = norm(rhs1 - matvec_y2mit(K, At, iHWy, dy, Lchol) ) / (1 + norm(rhs1));
    
    % fprintf("resdue of computing dy: %e\n", relres);
    flag = 1;
end


function detect_zeroA_nonzeroB(AHAt)
    % Assuming AHAt.mat22.At_sp is the sparse matrix and AHAt.mat22.diag_sp is the full vector
    At_sp = AHAt.mat22.At_sp;      % Sparse matrix of size m x n
    diag_sp = AHAt.mat22.diag_sp;  % Full vector of length n

    % Step 1: Compute A and B
    A = At_sp' .* diag_sp' * At_sp;    % Matrix A
    B = spones(At_sp)' * spones(At_sp);           % Sparsity pattern matrix B

    % Step 2: Find positions where A is zero but B is non-zero
    % Since A is sparse, (A == 0) checks both stored zeros and implicit zeros
    % To ensure we only consider positions where B is non-zero:
    % Extract non-zero positions in B
    [B_i, B_j] = find(B);

    % Extract corresponding values in A
    A_vals = A(sub2ind(size(A), B_i, B_j));

    % Identify positions where A is zero
    zeroA_nonzeroB_idx = find(A_vals == 0);

    % Get the row and column indices of these positions
    zeroA_nonzeroB_i = B_i(zeroA_nonzeroB_idx);
    zeroA_nonzeroB_j = B_j(zeroA_nonzeroB_idx);

    % Step 3: Iterate through each (i, j) to find and print contributing terms
    fprintf('Number of zeroA_nonzeroB_idx: %d.\n', length(zeroA_nonzeroB_idx));

    count = 0;

    for idx = 1:length(zeroA_nonzeroB_i)
        if count == 3
            break
        end

        i = zeroA_nonzeroB_i(idx);
        j = zeroA_nonzeroB_j(idx);
        
        % Find all k where At_sp(i, k) and At_sp(j, k) are non-zero
        [k_i, ~, ~] = find(At_sp(:, i));
        [k_j, ~, ~] = find(At_sp(:, j));
        k_common = intersect(k_i, k_j);
        
        % If no common k found, skip (though B(i,j) should indicate at least one)
        if isempty(k_common)
            fprintf('No common k found for (%d, %d).\n', i, j);
            continue;
        end
        
        % Extract c(k), v_i(k), and v_j(k)
        c_k = diag_sp(k_common);                % Coefficients c(k)
        v_i = At_sp(k_common, i);               % v_i(k)
        v_j = At_sp(k_common, j);               % v_j(k)
        
        % Compute the contributing terms
        terms = full(c_k) .* full(v_i) .* full(v_j);
        
        % Only consider non-zero terms to avoid unnecessary output
        nonzero_terms = terms ~= 0;
        terms = terms(nonzero_terms);
        k_common = k_common(nonzero_terms);
        v_i = v_i(nonzero_terms);
        v_j = v_j(nonzero_terms);

        if length(terms) == 0
            continue
        end

        count = count + 1;
        
        % Print the results
        fprintf('Position (%d, %d):\n', i, j);
        for t = 1:length(terms)
            fprintf('  Term %d: c(%d) * v(%d, %d) * v(%d, %d) = %.16f * %.16f * %.16f = %.16f\n', ...
                t, full(k_common(t)), i, full(k_common(t)), j, full(k_common(t)), full(c_k(t)), full(v_i(t)), full(v_j(t)), full(terms(t)));
        end
        fprintf('  Sum of terms: %.16f\n\n', sum(terms));
    end
end