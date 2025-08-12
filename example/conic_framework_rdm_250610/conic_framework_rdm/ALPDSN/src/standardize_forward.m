function [std_model, transform] = standardize_forward(model)
    %% standardize a general socp 

    % input problem is
    % min c' * x + c0
    % s.t.  lc <= A * x <= uc
    %     F * x + g in D
    %     lx <= x <= ux
    %

    % output problem is
    % min c' * x + c0
    % s.t. A * x = b
    %       x \in K
    % here K also contains box boundaries

    % assert(check_gen_model(model)) ;

    transform = struct;

    assert(all(model.lc <= model.uc));
    % remove rows with lc = -inf, uc = inf (i.e., unconstrained rows)
    unconstrained_rows = find(model.lc == -inf & model.uc == inf);
    model.A(unconstrained_rows, :) = [];
    model.lc(unconstrained_rows) = [];
    model.uc(unconstrained_rows) = [];
    clear unconstrained_rows;

    % % move rows of A with only one non-zero entry to lx and ux
    %% warning: this code is useful for beam30, beam7 (optsuite presolved)
    % boundary_rows = find(sum(model.A ~= 0, 2) == 1);
    % if numel(boundary_rows) > 10
    % A_temp = model.A(boundary_rows, :);
    % lc_temp = model.lc(boundary_rows);
    % uc_temp = model.uc(boundary_rows);

    % % Find the columns corresponding to the non-zero entries of boundary_rows
    % [rows, cols, nonzeros_entries] = find(A_temp);

    % % Logical index for positive entries in model.A
    % positive_idx = find(nonzeros_entries > 0);
    % non_positive_idx = find(nonzeros_entries < 0);

    % % Update lx ux for positive entries
    % new_lx = lc_temp(rows(positive_idx)) ./ nonzeros_entries(positive_idx);
    % [uniquecols, ~, ic] = unique(cols(positive_idx));
    % if ~isempty(ic)
    %     new_lx = accumarray(ic, new_lx, [], @max);
    %     model.lx(uniquecols) = max(model.lx(uniquecols), new_lx);
    %     new_ux = uc_temp(rows(positive_idx)) ./ nonzeros_entries(positive_idx);
    %     new_ux = accumarray(ic, new_ux, [], @min);
    %     model.ux(uniquecols) = min(model.ux(uniquecols), new_ux);
    % end
    % % Update lx ux for non-positive entries
    % new_lx = uc_temp(rows(non_positive_idx)) ./ nonzeros_entries(non_positive_idx);
    % [uniquecols, ~, ic] = unique(cols(non_positive_idx));
    % if ~ isempty(ic)
    %     new_lx = accumarray(ic, new_lx, [], @max);
    %     model.lx(uniquecols) = max(model.lx(uniquecols), new_lx);
    %     new_ux = lc_temp(rows(non_positive_idx)) ./ nonzeros_entries(non_positive_idx);
    %     new_ux = accumarray(ic, new_ux, [], @min);
    %     model.ux(uniquecols) = min(model.ux(uniquecols), new_ux);
    % end
    % model.A(boundary_rows, :) = [];
    % model.lc(boundary_rows) = [];
    % model.uc(boundary_rows) = [];
    % clear boundary_rows A_temp lc_temp uc_temp rows cols nonzeros_entries positive_idx non_positive_idx new_lx new_ux uniquecols ic;
    % end
    %% classify constraints and variables
    con_idx.e = find(model.lc == model.uc);                   % equation constraints
    con_idx.l = find( (model.lc ~= -inf & model.uc == inf) ); % lower bound constraints
    con_idx.u = find( (model.lc == -inf & model.uc ~= inf) ); % upper bound constraints
    con_idx.b = find( (model.lc ~= -inf & model.uc ~= inf & model.lc ~= model.uc) ); % boxed constraints

    %% we observe that in real-world data, most constraints are either equation or one-sided
    if ~isempty(con_idx.b)
        error('Boxed constraints are not supported yet');
    end
    b = zeros(size(model.A, 1), 1);
    b(con_idx.e) = model.lc(con_idx.e);
    % convert lower bound constraints to upper bound constraints
    model.A(con_idx.l, :) = - model.A(con_idx.l, :);
    b(con_idx.l) = - model.lc(con_idx.l);
    b(con_idx.u) = model.uc(con_idx.u);
    con_idx.u = union(con_idx.u, con_idx.l);
    con_idx.l = [];

    % convert inequality constraints to equality constraints by add slack variables
    % A * x <= b ==>  A * x + s = b, s >= 0
    A_slack = sparse(con_idx.u, 1:numel(con_idx.u) , ones(numel(con_idx.u), 1), size(model.A, 1), numel(con_idx.u));
    model.A = [model.A, A_slack];
    model.c = [model.c; zeros(numel(con_idx.u), 1)];
    model.F = [model.F, sparse(size(model.F, 1), numel(con_idx.u))];
    model.lx = [model.lx; zeros(numel(con_idx.u), 1)];
    model.ux = [model.ux; inf(numel(con_idx.u), 1)];
    
    var_idx.l = find( (model.lx ~= -inf & model.ux == inf) ); % lower bound variables
    var_idx.u = find( (model.lx == -inf & model.ux ~= inf) ); % upper bound variables
    var_idx.b = find( (model.lx ~= -inf & model.ux ~= inf & model.lx ~= model.ux) ); % boxed variables
    var_idx.f = find( (model.lx == -inf & model.ux == inf) ); % free variables
    var_idx.fix = find( (model.lx == model.ux) );             % fixed variables

    % convert upper bound variables to lower bound variables
    model.A(:, var_idx.u) = - model.A(:, var_idx.u);
    model.c(var_idx.u) = - model.c(var_idx.u);
    model.F(:, var_idx.u) = - model.F(:, var_idx.u);
    model.lx(var_idx.u) = - model.ux(var_idx.u);
    model.ux(var_idx.u) = inf;
    


    %% we find in most cases F is identy of identity concatenated with a zero matrixand together with g = 0, i.e., 
    % F = I or F = [I, 0] or F = [0, I]. We detect this and remove F and g
    std_model = struct;
    nonzero_cols = find(any(model.F ~= 0, 1));
    zero_cols = setdiff(1:size(model.F, 2), nonzero_cols);
    if is_identity(model.F(:, nonzero_cols)) && norm(model.g) < 1e-16
        if isempty(zero_cols)
            % fprintf('F is identity\n');
            std_model.K = model.D;
            std_model.At = MatCell({model.A'});
            std_model.c = MatCell({model.c});
        else
            % fprintf('F is identity concatenated with a zero matrix\n');
            std_model.K = [model.D;
            Cone({BasicCone('u', numel(zero_cols))})];
            std_model.At = MatCell({model.A(:, nonzero_cols)';
                        model.A(:, zero_cols)' });
            std_model.c = MatCell({model.c(nonzero_cols);
                    model.c(zero_cols)});
        end

        std_model.b = b;
        model.lx = [model.lx(nonzero_cols); model.lx(zero_cols)];
        model.ux = [model.ux(nonzero_cols); model.ux(zero_cols)];
    else
        % for general F or g, we need to add a new variable y = F * x + g
        std_model.K = [Cone({BasicCone('u', size(model.F, 2))}); % x free
                       model.D];                               % y \in D
        std_model.At = MatCell({[model.A', model.F'] ;
                        [sparse(size(model.F, 1), size(model.A, 1)), - speye(size(model.F, 1))] });
        std_model.c = MatCell({model.c; 
                       zeros(size(model.F, 1), 1)});
        std_model.b = [b; - model.g];
        model.lx = [model.lx; -inf(size(model.F, 2), 1)];
        model.ux = [model.ux; inf(size(model.F, 2), 1)];
        
    end

    %% now we try to remove variables bounds or put them into K
    i = 1;
    c0 = model.c0;
    p = 1;
    while p <= length(std_model.K)
        cone = std_model.K{p};
        len = sum(cone.size);
        if strcmp(cone.type, 'u')
            if all(model.lx(i: i + len - 1) == -inf) && all(model.ux(i: i + len - 1) == inf)
                % free variables
            else % bounded variables
                lx = model.lx(i: i + len - 1);
                ux = model.ux(i: i + len - 1);
                var_idx = struct; % note that upper bound variables has been converted to lower bound variables
                var_idx.l = find( (lx ~= -inf & ux == inf) ); % lower bound variables
                var_idx.b = find( (lx ~= -inf & ux ~= inf) & (lx ~= ux )); % boxed variables
                var_idx.f = find( (lx == -inf & ux == inf) ); % free variables
                var_idx.fix = find( lx == ux ); % fixed variables
                if ~isempty(var_idx.fix)
                    % fix the fixed variables
                    c0 = c0 + std_model.c{p}(var_idx.fix)'*lx(var_idx.fix);
                    std_model.b = std_model.b - std_model.At{p}(var_idx.fix, :)'*lx(var_idx.fix);
                end
                
                % split the cone into several cones
                n_cones = ( ~ isempty(var_idx.l) ) + ( ~ isempty(var_idx.b) ) + ( ~ isempty(var_idx.f) );
                new_cones = cell(0, 1);;
                new_At = cell(0, 1);
                new_c = cell(0, 1);
                if ~ isempty(var_idx.l)
                    c0 = c0 + std_model.c{p}(var_idx.l)'*lx(var_idx.l);
                    std_model.b = std_model.b - std_model.At{p}(var_idx.l, :)'*lx(var_idx.l);
                    new_cones = [new_cones; Cone({BasicCone('l', numel(var_idx.l))})];
                    new_At = [new_At; {std_model.At{p}(var_idx.l, :)}];
                    new_c = [new_c; {std_model.c{p}(var_idx.l)}];
                end
                if ~ isempty(var_idx.f)
                    new_cones = [new_cones; Cone({BasicCone('u', numel(var_idx.f))})];
                    new_At = [new_At; {std_model.At{p}(var_idx.f, :)}];
                    new_c = [new_c; {std_model.c{p}(var_idx.f)}];
                end
                if ~ isempty(var_idx.b)
                    new_cones = [new_cones; Cone({BasicCone('b', numel(var_idx.b), [lx(var_idx.b), ux(var_idx.b)])})];
                    new_At = [new_At; {std_model.At{p}(var_idx.b, :)}];
                    new_c = [new_c; {std_model.c{p}(var_idx.b)}];
                end

                std_model.K = [std_model.K(1: p - 1); new_cones; std_model.K(p + 1: end)];
                std_model.At = [std_model.At(1: p - 1); new_At; std_model.At(p + 1: end)];
                std_model.c = [std_model.c(1: p - 1); new_c; std_model.c(p + 1: end)];
                p = p + n_cones - 1;

            end
        else % other cones in D
            if all(model.lx(i: i + len - 1) == -inf) && all(model.ux(i: i + len - 1) == inf)
                % no other bounds
            elseif strcmp(std_model.K{p}.type, 'q') && all(model.lx(i: i + len - 1) == 0 | model.lx(i: i + len - 1) == -inf) && all(model.ux(i: i + len - 1) == inf) && model.lx(i) == 0 && all(model.lx(i + std_model.K{p}.size(1: end-1) ) == 0)
                % quadratic cone with only the head index have lower bound 0 in each sub-cone
            else
                std_model.K{p}.type = [std_model.K{p}.type,'+', 'b'];
                std_model.K{p}.params = [model.lx(i: i + len - 1), model.ux(i: i + len - 1)];
            end
        end
        p = p + 1;
        i = i + len;
    end

    % remove empty rows in A
    m = size(std_model.b, 1);
    row_normA = zeros(m, 1);
    for p = 1: length(std_model.K)
        row_normA = row_normA + (vecnorm(std_model.At{p}, 2, 1)') .^ 2;
    end
    row_normA = sqrt(row_normA);

    empty_row = find(row_normA == 0);
    if ~isempty(empty_row)
        % fprintf("remove %d empty rows in A.\n", length(empty_row) );
        for p = 1: length(std_model.K)
            std_model.At{p}(:, empty_row) = [];
        end
        std_model.b(empty_row) = [];
    end
    clear row_normA empty_row;

  
end


function flag = is_identity(F)
    if size(F, 1) ~= size(F, 2)
        flag = false;
        return;
    end
    flag = (norm(F - speye(size(F, 1)), 'fro') < 1e-16);
end