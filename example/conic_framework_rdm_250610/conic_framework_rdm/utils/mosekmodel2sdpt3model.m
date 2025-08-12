function sdpt3model = mosekmodel2sdpt3model(mosekmodel)
    %% transform a mosek model to sdpt3 model

    % input problem is
    % min c' * x + cfix
    % s.t.  blc <= a * x <= buc
    %     f * x + g in accs
    %     blx <= x <= bux


    % tranform accs to blk
    assert(mod(length(mosekmodel.accs), 2) == 0, "length(mosekmodel.accs) must be even.");
    accs_reshape = reshape(mosekmodel.accs, 2, [])';
    cone_type = accs_reshape(:, 1);
    for i = 1:length(cone_type)
        assert(cone_type(i) >= 0 && cone_type(i) <= 5, " unsupported cone type." + cone_type(i));
    end
    % cone_type = 0: free cone
    % cone_type = 1: zero cone
    % cone_type = 2: nonnegative cone
    % cone_type = 3: nonpositive cone
    % cone_type = 4: quadratic cone
    % cone_type = 5: rotated quadratic cone     
    
    row_type = repelem(cone_type, accs_reshape(:, 2));
    row_idx0 = find(row_type == 0);
    row_idx1 = find(row_type == 1);
    row_idx2 = find(row_type == 2);
    row_idx3 = find(row_type == 3);
    row_idx4 = find(row_type == 4);
    row_idx5 = find(row_type == 5);
    cone_idx0 = find(cone_type == 0);
    cone_idx1 = find(cone_type == 1);
    cone_idx2 = find(cone_type == 2);
    cone_idx3 = find(cone_type == 3);
    cone_idx4 = find(cone_type == 4);
    cone_idx5 = find(cone_type == 5);

    % create a new my model
    model = struct();

    % move linear blocks (including zero cone) to A matrix
    model.A = [mosekmodel.a; mosekmodel.f(row_idx1, :); mosekmodel.f(row_idx2, :); mosekmodel.f(row_idx3, :)];
    model.lc = [mosekmodel.blc; - mosekmodel.g(row_idx1); - mosekmodel.g(row_idx2); -inf(length(row_idx3), 1)];
    model.uc = [mosekmodel.buc; - mosekmodel.g(row_idx1); inf(length(row_idx2), 1); - mosekmodel.g(row_idx3)];

    % concatenate the second order cones 
    model.F = mosekmodel.f([row_idx4; row_idx5], :);
    model.g = mosekmodel.g([row_idx4; row_idx5]);
    model.blk_D = {'q', reshape(accs_reshape(cone_idx4, 2), 1, []);
                'r', reshape(accs_reshape(cone_idx5, 2), 1, [])};

    % remove empty blocks
    model.blk_D = model.blk_D(~cellfun(@(cone_size) sum(cone_size) == 0, model.blk_D(:, 2)), :);

    % copy the rest of the data
    model.c0 = mosekmodel.cfix;
    model.c = mosekmodel.c;
    model.lx = mosekmodel.blx;
    model.ux = mosekmodel.bux;

    % standardize my model to be standard form
    std_model = standardize(model);


    sdpt3model = struct();
    sdpt3model.At_mat = std_model.At_int;
    sdpt3model.C_mat = std_model.c_int;
    sdpt3model.b = std_model.b_int;
    sdpt3model.blk = std_model.blk_int;

    % transform At, C to cell
    sdpt3model.At = mat2cell(sdpt3model.At_mat, cellfun(@(cone_size) sum(cone_size), sdpt3model.blk(:, 2)), size(sdpt3model.At_mat, 2));
    sdpt3model.C = mat2cell(sdpt3model.C_mat, cellfun(@(cone_size) sum(cone_size), sdpt3model.blk(:, 2)), 1);
 

    
   
end