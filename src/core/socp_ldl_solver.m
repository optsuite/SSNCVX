

function [dy,relres,flag] = socp_ldl_solver(opts,rhsz,optcg)

% function [dy, relres, flag] = solve_direct_opt2(K, At, rhs1, iHWy, Lchol, trans)
    K = opts.K;
    At = opts.At;
    rhs1 = rhsz;
    iHWy = optcg;
    m = size(At{1}, 2);
    Lchol.isidentity = 1;
    Lchol.Rperm = speye(m);
    Lchol.RtR = speye(m);
    Lchol.R = speye(m);
    Lchol.Rt = speye(m);
    Lchol.perm = 1:m;

    trans.system_opt = 2;
    trans.socp_formula = 1;
    trans.sys2_sparse_strategy = 0;
    trans.threshold = 0.1;
    trans.direct_solver = 'ldl';

    sp_info = detect_sp_info(At, K, Lchol, trans.system_opt, trans.socp_formula, trans.sys2_sparse_strategy, trans.threshold);
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
    rhs1_sp = rhs1_temp(row_sp,:);
    rhs1_den = rhs1_temp(row_den,:);
    mat21_tilde = [AHAt.mat21, rhs1_sp];

    % mat22_inv_mat21_tilde = AHAt.mat22.mat \ mat21_tilde;
    C_vec = sqrt(abs(AHAt.mat22.diag_den));
    lhs.mat11 = AHAt.mat22.At_sp' .* AHAt.mat22.diag_sp' * AHAt.mat22.At_sp;
    lhs.mat12 = sparse(AHAt.mat22.At_den' .* C_vec');
    lhs.mat22 = -spdiags(1./AHAt.mat22.diag_den .* C_vec .* C_vec, 0, n_den, n_den);
    augmented_lhs = [lhs.mat11, lhs.mat12; lhs.mat12', lhs.mat22];
    augmented_rhs = [mat21_tilde; zeros(n_den, m_den + 1)];

    if strcmp(direct_solver, 'backslash')
        augmented_x = augmented_lhs \ augmented_rhs;
    elseif strcmp(direct_solver, 'ldl')
        [L_ldl, D_ldl, P_ldl] = ldl(augmented_lhs, 'vector'); % 'vector' for permutation vector
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


function sp_info = detect_sp_info(At, K, Lchol, system_opt, socp_formula, sys2_sparse_strategy, threshold, fid)

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

    At_row_den = cellfun(@(x) full(x(:, row_den)), At, 'UniformOutput', false)';
    At_row_sp = cellfun(@(x) sparse(x(:, row_sp)), At, 'UniformOutput', false)';
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

    AHAt_mat22_At = [At_row_sp_merged; AHAt_mat22_At_new; Lchol.Rperm(:, row_sp)];

    if system_opt == 2
        nnz_per_row = sum(spones(AHAt_mat22_At), 2);
    else
        nnz_per_row = zeros(size(AHAt_mat22_At, 1), 1);
    end

    if sys2_sparse_strategy == 0
        sparsity = nnz_per_row / size(AHAt_mat22_At, 2);
        col_den = find(sparsity > threshold);
        col_sp = find(sparsity <= threshold);
    else
        col_den = find(nnz_per_row > sys2_sparse_strategy);
        col_sp = find(nnz_per_row <= sys2_sparse_strategy);
    end

    % AHAt_mat22_At_den = AHAt_mat22_At(col_den, :);
    % AHAt_mat22_At_sp = AHAt_mat22_At(col_sp, :);
    % AHAt_mat22_diag_sp = rand(length(col_sp), 1) + 1;
    % lhs_structure = [AHAt_mat22_At_sp' .* AHAt_mat22_diag_sp' * AHAt_mat22_At_sp, AHAt_mat22_At_den'; AHAt_mat22_At_den, speye(length(col_den))];
    % lhs_structure = spones(lhs_structure);

    sp_info = struct('row_den', row_den, ...
                     'row_sp', row_sp, ...
                     'row_index', [row_den, row_sp], ...
                     'col_den', col_den, ...
                     'col_sp', col_sp, ...
                     'At_row_den', {At_row_den}, ...
                     'At_row_sp', {At_row_sp}, ...
                     'At_row_sp_merged', At_row_sp_merged);
                     %'AHAt_mat22_At', AHAt_mat22_At, ...
                     %'lhs_structure', lhs_structure);

    % if system_opt == 2
    %     fprintf(fid, 'System_opt = 2\n');
    %     fprintf(fid, 'Number of dense rows: %d\n', length(row_den));
    %     fprintf(fid, 'Number of sparse rows: %d\n', length(row_sp));
    %     fprintf(fid, 'Number of dense columns (new At): %d\n', length(col_den));
    %     fprintf(fid, 'Number of sparse columns (new At): %d\n', length(col_sp));
    %     % fprintf(fid, 'Number of non-zero elements in lhs_structure: %d\n', nnz(lhs_structure));
    % elseif system_opt == 0
    %     fprintf(fid, 'System_opt = 0\n');
    %     % fprintf(fid, 'Number of non-zero elements in lhs_structure: %d\n', nnz(lhs_structure));
    % end
end


function By = matvec_y2mit(K,At,par,y,AL)

    if (nargin < 5); AL = []; end;
    if isempty(AL); existAL = 0; else; existAL = 1; end
    
    
    N = length(y);
    if (norm(y) == 0); By = zeros(N,1); return; end
    %%
    yorg = y;
    if ~AL.isidentity && (existAL)
        if strcmp(AL.matfct_options,'chol')
            y(AL.perm) = AL.R \ y;
        elseif strcmp(AL.matfct_options,'spcholmatlab')
            y(AL.perm) = mexbwsolve(AL.Rt,y);
        end
    end
    %%
    sumlen = 0;
    By = zeros(N,1);
    for p = 1:length(K)
        cone = K{p};
        % sumlen = sumlen + length(cone.size);
        if strcmp(cone.type,'s')
            cumsumlen = [0 cumsum(cone.size)];
            cumsumsquare = (cone.size.^2 + cone.size)/2;
            cumsumsquare = [0 cumsum(cumsumsquare )];
            %          Aty = Atyfun_sdpnal(pblk,At{p},y);
            Aty = Atyfun(cone,At{p},y);
            for  jj = 1:length(cone.size)
                n = cone.size(jj);
                rr = size(par.P1{sumlen+jj},2);
                if (rr > 0 && rr < n)
                    if (rr <= n/2)
                        tmp0 = par.P1{sumlen+jj}'*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                        tmp1 = (tmp0*par.P1{sumlen+jj})*par.P1{sumlen+jj}';
                        tmp2 = par.Dsch12{sumlen+jj}.*(tmp0*par.P2{sumlen+jj});
                        tmp2 = tmp2*par.P2t{sumlen+jj};
                        tmp3 = par.P1{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
                        tmpcone = cone;
                        tmpcone.size = cone.size(jj);
                        %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},tmp3+tmp3');
                        By = By + par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),tmp3+tmp3');
                    else
                        tmp0 = par.P2t{sumlen+jj}*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                        tmp1 = (tmp0*par.P2{sumlen+jj})*par.P2t{sumlen+jj};
                        tmp2 = (par.Dsch1{sumlen+jj}-par.Dsch12{sumlen+jj}').*(tmp0*par.P1{sumlen+jj});
                        tmp2 = tmp2*par.P1t{sumlen+jj};
                        tmp3 = par.P2{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
                        tmpcone = cone;
                        tmpcone.size = cone.size(jj);
                        %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},par.Dsch1{p}*Aty-tmp3-tmp3');
                        By = By + par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),par.Dsch1{sumlen+jj}*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))-tmp3-tmp3');
                    end
                elseif (rr == n)
                    tmpcone = cone;
                    tmpcone.size = cone.size(jj);
                    %             By = By + par.sig*AXfun_sdpnal(pblk,At{p},Aty);
                    By = By + par.Dsch1{sumlen+jj}*par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1)));
                end
            end
        elseif strcmp(cone.type,'q')
            if ~isfield(par, 'socp_formula') || par.socp_formula == 1 || par.socp_formula == 2
                Aty = At{p}*y;
                tmp1 = par.Q1{p}'*Aty.*par.Dsch1{sumlen+1};
                tmp1 = repelem(tmp1,cone.size,1).*par.P1{sumlen+1};
                tmp2 = par.Q2{p}'*Aty.*par.Dsch2{sumlen+1};
                tmp2 = repelem(tmp2,cone.size,1).*par.P2{sumlen+1};
                tmp3 = repelem(par.shift{p},cone.size,1).*Aty;
                tmp = tmp1+tmp2+tmp3;
                By = By + par.sig*(tmp'*At{p})';
            elseif par.socp_formula == 3
                Aty = At{p}*y;
                tmp1 = par.Q3{p}'*Aty.*par.Dsch3{p};
                tmp1 = repelem(tmp1,cone.size,1).*par.P3{p};
                tmp3 = par.shift3{p}.*Aty;
                tmp = tmp1+tmp3;
                By = By + par.sig*(tmp'*At{p})';
            end
        elseif strcmp(cone.type,'l') || strcmp(cone.type,'b') || strcmp(cone.type,'u') || strcmp(cone.type,'b2l')
            if (~isempty(par.Dsch12{sumlen+1}))
                tmp = par.Dsch12{sumlen+1}.*(At{p}*y);
                By = By + par.sig*(tmp'*At{p})';
            end
        elseif strcmp(cone.type,'u')
            tmp = At{p}*y;
            By = By + par.sig*(tmp'*At{p})';
        end
        sumlen = sumlen + strcmp(cone.type, 's') * length(cone.size) + (1 - strcmp(cone.type, 's')) * 1;
    end
    if ~AL.isidentity && (existAL)
        if strcmp(AL.matfct_options,'chol')
            By = AL.Rt \ By(AL.perm);
        elseif strcmp(AL.matfct_options,'spcholmatlab')
            By = mexfwsolve(AL.R,By(AL.perm));
        end
    end
    %
    %    if (par.use_proximal);
    %       By = By + par.H2.*(y/par.sighat);
    %    else
    %       sighat = max([1e4,10*par.sig]);
    %       By = By + y/sighat;
    %    end
    
    %    if isfield(par,'epsilon')
    By = By + par.epsilon*yorg;
    
    %    end
    %%************************************************************************
end