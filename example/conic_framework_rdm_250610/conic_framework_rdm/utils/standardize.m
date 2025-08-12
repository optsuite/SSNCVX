function model = standardize(model)
    %% standardize a general socp 

    % input problem is
    % min c' * x + c0
    % s.t.  lc <= A * x <= uc
    %     F * x + g in D
    %     lx <= x <= ux

    %% classify constraints and variables

    t0 = tic;
    model.standardize = struct;


    
    assert(all(model.lc <= model.uc));
    % remove rows with lc = -inf, uc = inf (i.e., unconstrained rows)
    unconstrained_rows = find(model.lc == -inf & model.uc == inf);
    model.A(unconstrained_rows, :) = [];
    model.lc(unconstrained_rows) = [];
    model.uc(unconstrained_rows) = [];

    % remove cols with lx == ux (i.e., fixed variables)
    fixed_cols = find(model.lx == model.ux);
    model.g = model.g + model.F(:, fixed_cols) * model.lx(fixed_cols);
    model.F(:, fixed_cols) = [];
    model.c0 = model.c0 + sum(model.lx(fixed_cols) .* model.c(fixed_cols));
    model.c(fixed_cols) = [];
    model.A(:, fixed_cols) = [];
    model.lx(fixed_cols) = [];
    model.ux(fixed_cols) = [];
    
    %% classify constraints and variables
    model.con_idx_e = find(model.lc == model.uc);                   % equation constraints
    model.con_idx_ie = find(model.lc ~= model.uc);                  % inequality constraints
    model.con_idx_l = find( (model.lc ~= -inf & model.uc == inf) ); % lower bound constraints
    model.con_idx_u = find( (model.lc == -inf & model.uc ~= inf) ); % upper bound constraints
    model.con_idx_b = find( (model.lc ~= -inf & model.uc ~= inf & model.lc ~= model.uc) ); % boxed constraints

    model.n_con = size(model.lc, 1);
    model.n_con_e = length(model.con_idx_e);
    model.n_con_ie = length(model.con_idx_ie);
    model.n_con_l = length(model.con_idx_l);
    model.n_con_u = length(model.con_idx_u);
    model.n_con_b = length(model.con_idx_b);

    model.var_idx_l = find( (model.lx ~= -inf & model.ux == inf) ); % lower bound variables
    model.var_idx_u = find( (model.lx == -inf & model.ux ~= inf) ); % upper bound variables
    model.var_idx_b = find( (model.lx ~= -inf & model.ux ~= inf) ); % boxed variables
    model.var_idx_f = find( (model.lx == -inf & model.ux == inf) ); % free variables
    model.var_idx_nf = find( ~ (model.lx == -inf & model.ux == inf) ); % not free variables


    model.n_var = size(model.lx, 1);
    model.n_var_l = length(model.var_idx_l);
    model.n_var_u = length(model.var_idx_u);
    model.n_var_b = length(model.var_idx_b);
    model.n_var_f = length(model.var_idx_f);
    model.n_var_nf = length(model.var_idx_nf);
    model.n_var_soc = size(model.F, 1);

    %% Introduce slack variables
    %compute transformation matrix for constraints and variables
    model.standardize.Lambda_w = ones(model.n_con, 1);
    model.standardize.Lambda_w(model.con_idx_u) = -1;
    model.standardize.Lambda_w_ie = model.standardize.Lambda_w(model.con_idx_ie);
    model.standardize.Lambda_x = ones(model.n_var, 1);
    model.standardize.Lambda_x(model.var_idx_u) = -1;
    model.standardize.h_w = model.lc;
    model.standardize.h_w(model.con_idx_u) = model.uc(model.con_idx_u);
    model.standardize.h_x = model.lx;
    model.standardize.h_x(model.var_idx_u) = model.ux(model.var_idx_u);
    model.standardize.h_x(model.var_idx_f) = 0;

    %standardize constraints and variables
    model.c = model.standardize.Lambda_x .* model.c;
    model.c0 = model.c0 + model.standardize.h_x' * model.c;
    model.b = model.standardize.Lambda_w .* (model.standardize.h_w - model.A * model.standardize.h_x);
    model.A = spdiag(model.standardize.Lambda_w) * model.A * spdiag(model.standardize.Lambda_x);
    model.g = model.g + model.F * model.standardize.h_x;
    model.F = model.F * spdiag(model.standardize.Lambda_x); 


    % Here we obtain the following standardized form SOCP
    % min c' * x + c0
    % s.t. A * x - w = b
    %      F * x - y = g
    %      x(var_idx_box) + x_bar(var_idx_box) = ux(var_idx_box) - lx(var_idx_box)
    %      w(con_idx_box) + w_bar(con_idx_box) = uc(con_idx_box) - lc(con_idx_box)
    %      w >= 0, w_bar(con_idx_box) >= 0
    %      x(var_idx_nf) >= 0, x_bar(var_idx_box) >= 0, x(var_idx_f) free
    %      y in D


    %% Utilize the structure of F
    % find digonal subblock in F
    F_pattern = (model.F ~= 0);
    nnz_per_row = full(sum(F_pattern, 2));
    cadidate_row = find(nnz_per_row == 1);
    nnz_per_col = sum(F_pattern, 1);
    cadidate_col = intersect(find(nnz_per_col == 1), model.var_idx_f);

    Omega = reshape(cadidate_col, length(cadidate_col), 1);
    Phi = zeros(length(cadidate_col), 1);
    for i = 1: length(cadidate_col)
        Phi(i) = find(F_pattern(:, cadidate_col(i)));
    end
    clear candiate_row candiate_col nnz_per_row nnz_per_col F_pattern


    Phi_bar = setdiff([1: size(model.F, 1)]', Phi);
    Omega_bar = setdiff(model.var_idx_f, Omega);

    %% eliminate
    model.F_Phi_Omega = diag(model.F(Phi, Omega));
    model.A(:, Omega) = model.A(:, Omega) * spdiag(1./ model.F_Phi_Omega);
    model.F = model.F(Phi_bar, :);
    model.F(:, Omega) = 0;
    %% q part
    A_q = sparse(model.n_con, model.n_var_soc);
    A_q(:, Phi) = model.A(:, Omega);
    F_q = sparse(size(model.F, 1), model.n_var_soc);
    F_q(:, Phi_bar) = - speye(length(Phi_bar));

    %% Omega is the free variables
    model.var_idx_f = Omega_bar;
    model.var_idx_q = Omega;
    clear Omega_bar Omega 

    %% rearrange varaibles
    model.A = [model.A(:, model.var_idx_nf), model.A(:, model.var_idx_f), model.A(:, model.var_idx_q)];
    model.F = [model.F(:, model.var_idx_nf), model.F(:, model.var_idx_f), model.F(:, model.var_idx_q)];
    model.lx = [model.lx(model.var_idx_nf); model.lx(model.var_idx_f); model.lx(model.var_idx_q)];
    model.ux = [model.ux(model.var_idx_nf); model.ux(model.var_idx_f); model.ux(model.var_idx_q)];
    model.c = [model.c(model.var_idx_nf); model.c(model.var_idx_f); model.c(model.var_idx_q)];

    model.var_idx_nf = [1: length(model.var_idx_nf)]';
    model.var_idx_f = [length(model.var_idx_nf) + 1: length(model.var_idx_nf) + length(model.var_idx_f)]';
    model.var_idx_q = [length(model.var_idx_nf) + length(model.var_idx_f) + 1: length(model.var_idx_nf) + length(model.var_idx_f) + length(model.var_idx_q)]';




    %% concatenate A and F
    model.Abar = [model.A; 
                   model.F];


    % for csc format , use At is faster for computing A * D * At in newton system
    model.Atbar = model.Abar';

    % fprintf("[after stadandize] Abar norm: %e\n", norm(model.Abar, 'fro'));


    %% compute variable index in not free part
    model.var_nf_idx_b = find( (model.lx(model.var_idx_nf) ~= -inf & model.ux(model.var_idx_nf) ~= inf & model.lx(model.var_idx_nf) ~= model.ux(model.var_idx_nf)) ); % boxed variables
    model.Atbar_b = model.Abar(:, model.var_nf_idx_b)';



    %% compute constraint index in inequility part
    model.lc_ie = model.lc(model.con_idx_ie);
    model.uc_ie = model.uc(model.con_idx_ie);
    model.con_ie_idx_l = find( (model.lc_ie ~= -inf & model.uc_ie == inf) ); % lower bound constraints
    model.con_ie_idx_u = find( (model.lc_ie == -inf & model.uc_ie ~= inf) ); % upper bound constraints
    model.con_ie_idx_b = find( (model.lc_ie ~= -inf & model.uc_ie ~= inf & model.lc_ie ~= model.uc_ie) ); % boxed constraints

    %% build the internal representation of the problem

    temp1 = sparse(model.n_con, model.n_con_ie);
    temp1(model.con_idx_ie, :) = speye(model.n_con_ie);
    temp2 = sparse(model.n_var_b, model.n_var_nf);
    temp2(:, model.var_nf_idx_b) = speye(model.n_var_b);
    temp3 = sparse(model.n_con_b, model.n_con);
    temp3(:, model.con_idx_b) = speye(model.n_con_b);
    temp3 = temp3(:, model.con_idx_ie);

    % "int" stand for "internal" here
    model.At_int = [model.Atbar,                      temp2', sparse(model.n_con_b, model.n_var_nf + length(model.var_idx_q) + model.n_var_soc);
                    -temp1', sparse(size(model.F, 1), model.n_con_ie + model.n_var_b + model.n_con_b), sparse(model.n_var_b, model.n_var_f), sparse(model.n_con_b,  model.n_var_b);
                    sparse(model.n_con, model.n_var_b + model.n_con_b), sparse(model.n_var_b, model.n_con_ie + model.n_var_soc), speye(model.n_var_b), temp3';
                    sparse(model.n_var_b + model.n_con_b, model.n_var_nf + length(model.var_idx_q) + model.n_var_soc), speye(model.n_con_b), sparse(model.n_var_b, model.n_var_f), sparse(model.n_con_b, model.n_con_ie + model.n_var_soc)];    model.b_int = [model.b + model.A(:, model.var_idx_q) * model.g(Phi, :);
                        - model.g(Phi_bar); 
                        model.ux(model.var_idx_b) - model.lx(model.var_idx_b); 
                        model.uc(model.con_idx_b) - model.lc(model.con_idx_b)];
    c_y = zeros(model.n_var_soc, 1);
    c_y(Phi) = model.c(model.var_idx_q) ./ model.F_Phi_Omega;
    model.c0_int = model.c0 - sum(c_y(Phi) .* model.g(Phi));
    model.c_int = [model.c(model.var_idx_nf);
                   model.c(model.var_idx_f);
                    c_y;
                    zeros(model.n_con_ie, 1);
                    zeros(model.n_var_b, 1);
                    zeros(model.n_con_b, 1)];

    model.blk_int = [{'l', model.n_var_nf;    % x_nf
                  'u', length(model.var_idx_f)};     % x_f
                   model.blk_D;           % y
                  {'l', model.n_con_ie;   % w
                  'l', model.n_var_b;    % x_bar(var_idx_box)
                  'l', model.n_con_b}];   % w_bar(con_idx_box)

    
    % remove empty blocks
    model.blk_int = model.blk_int(~cellfun(@(cone_size) sum(cone_size) == 0, model.blk_int(:, 2)), :);
    
    
    model.int_idx_f = [model.n_var_nf + 1 : model.n_var_nf + length(Phi_bar)]'; % index of free variables in internal model, i.e, Omega_bar
    model.int_idx_q = [model.n_var_nf + length(Phi_bar) + 1 : model.n_var_nf + length(Phi_bar) + model.n_var_soc]'; % index of quadratic variables in internal model
    model.int_idx_l =   [[1: model.n_var_nf], ...
                [model.n_var_nf + length(Phi_bar) + model.n_var_soc + 1 : model.n_var_nf + length(Phi_bar) + model.n_var_soc + model.n_con_ie]]'; % index of linear variables in internal model

    model.standardize.time = toc(t0);
end