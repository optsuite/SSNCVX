function [model, model_original] =  data2model(dataset,probname,dir_data,options)
model_original = struct();
if dataset == "theta"
    file = fullfile(dir_data, dataset, probname + "_Alt.mat");
    [blk, At, C, b] = thetaread(file);
%     [b,At,cnz] = data_process(blk,At,b);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
elseif dataset == "RDM"
    [~, basename1, basename2] = fileparts(probname);
    basename = [basename1, basename2];
    file = [dir_data, '/RDM/', probname, '.mat'];
    load(file);
    [b,At,cnz] = data_process(blk,At,b);
    model.At = MatCell(At);
%     model.At = At;
    model.b = b;
     model.C =C;
     model.blk = blk;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
    model.basename = basename;
elseif dataset == "R1TA"
    [~, basename1, basename2] = fileparts(probname);
    basename = [basename1, basename2];
    file = [dir_data, '/R1TA/content/', probname, '.mat'];
%     [b,At,cnz] = data_process(blk,At,b);
    load(file);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
elseif dataset == "thetaplus"
    dataset2 = 'theta';
    file = fullfile(dir_data, dataset2, probname + "_Alt.mat");
    [blk, At, C, b] = thetaread(file);
    
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
    model.L = {0};
    model.U = {inf};
elseif dataset == "qap"
    file = fullfile(dir_data, dataset, probname + ".dat");
    [A, B] = qapread(file);
    [blk,At,C,b,Ascale,Bscale] = qapAW_sdpnal(A, B);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
elseif dataset == "rcp"
    K = options.K;
     file = [dir_data,'rcp/' ,probname, '.data'];
    [blk,At,C,b] = rcpread_sdpnal(file, K, options.n0);
     model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
    model.L = 0;
    model.U = inf;
elseif dataset == "biq"
    file = [dir_data, '/biq_all/', probname, '.sparse'];
%      file = [dir_data,'\rcp\' ,probname, '.data'];
        Q = biqread_sdpnal(file);
    [blk,At,A,C,b,Bt,B,d] = biq_ineq_sdpnal(Q);
     model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
    model.L = {0};
    model.U = {inf};
elseif dataset == "fap"
    file = [dir_data, '/fap/', probname, '.dat'];
    [blk,At,C,b,L,U] = fapread_lu_sdpnal(file);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = Cone.fromblk(blk);
    model.L = L;
    model.U = U;
elseif dataset == "sedumi"
    try
    model2 = load([dir_data,'/sedumi/',probname,'.mat']);
    if isfield(model2,'At')
    [blk,At,C,b] = read_sedumi(model2.At',model2.b,model2.c(:),model2.K);
    [b,At,cnz] = data_process(blk,At,b);
    elseif isfield(model2,'A')
    [blk,At,C,b] = read_sedumi(model2.A,model2.b,model2.c(:),model2.K);
%     [b,At,cnz] = data_process(blk,At,b);
    elseif isfield(model2,'AT')
    [blk,At,C,b] = read_sedumi(model2.AT',model2.b,model2.c(:),model2.K);
    end
    catch
    [blk,At,C,b] = read_sdpa([dir_data,'/sedumi/',probname,'.dat-s']);  
    [blk_new,At_new,C_new, c, hists] = rdm_detect_block(blk,At,C,b);
    At = At_new;
    blk = blk_new;
    C = C_new;
    b = c;
    [b,At,cnz] = data_process(blk,At,b);
    end

    model.K = Cone.fromblk(blk);
    model.C = MatCell(C);
    
    model.At = MatCell(At);
    model.b = b;
elseif dataset == "DIMACS"
    model2 = load(fullfile(dir_data, dataset, probname + ".mat"));
    if isfield(model2,'At')
        [blk,At,C,b] = read_sedumi(model2.At',model2.b,model2.c(:),model2.K);
    elseif isfield(model2,'A')
        [blk,At,C,b] = read_sedumi(model2.A,model2.b,model2.c(:),model2.K);
    elseif isfield(model2,'AT')
        [blk,At,C,b] = read_sedumi(model2.AT',model2.b,model2.c(:),model2.K);
    end
    model.K = Cone.fromblk(blk);
    model.C = MatCell(C);
    
    model.At = MatCell(At);
    model.b = b;
    model_original = model;
elseif dataset == "CBLIB"
    if isfile(fullfile(dir_data, dataset, probname + ".mat"))
        load(fullfile(dir_data, dataset, probname + ".mat"));
        model = standardize(model);

        mystdmodel = struct();
        mystdmodel.At_mat = model.At_int;
        mystdmodel.b = model.b_int;
        mystdmodel.K = Cone.fromblk(model.blk_int);
        mystdmodel.At = MatCell.vert_split(mystdmodel.At_mat, size(mystdmodel.K));
        mystdmodel.C = MatCell.vert_split(model.C, size(mystdmodel.K));


        model = mystdmodel;
    elseif isfile(fullfile(dir_data, dataset, probname + ".cbf.gz"))
        filename = fullfile(dir_data, dataset, probname + ".cbf.gz");
        [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
        prob = res.prob;
        model = mosekmodel2sdpt3model(prob);
        model.At = MatCell(model.At);
        model.C = MatCell(model.C);
        model.b = model.b;
        model.K = Cone.fromblk(model.K);
        model.name = probname;
    end
elseif dataset == "sdplib"
    if isfile(fullfile(dir_data, dataset, probname + ".dat-s"))
        [blk, At, C, b] = read_sdpa(fullfile(dir_data, dataset, probname + ".dat-s"));

        model = struct();
        model.b = b;
        model.K = Cone.fromblk(blk);
        model.At = MatCell(At);
        model.C = MatCell(C);

    else
        error("No such file: %s", fullfile(dir_data, dataset, probname + ".dat-s"));
    end
elseif dataset == "gen"
    if isfile(fullfile(dir_data, dataset, probname + ".mat"))
        model_orig = load(fullfile(dir_data, dataset, probname + ".mat"));
        [model, transform] = standardize_forward(model_orig.model);
        model.C = model.c;
        model = rmfield(model, 'c');
        model_original = model;
        model = convert_b_cone(model);
    else
        error("No such file: %s", fullfile(dir_data, dataset, probname + ".mat"));
    end
elseif dataset == "pass"
    model_orig = load(fullfile(dir_data, dataset, probname + ".mat"));
    % model = model_orig.model;
    model = model_orig;
    model.C = model.c;
    model = rmfield(model, 'c');
    model_original = model;
    model = convert_b_cone(model);
end
if isempty(model_original)
    model_original = model;
end
end


% Function: convert_b_cone
function model = convert_b_cone(model)
    % Initialize new model components
    new_K = {};
    new_C = {};
    new_At = {};
    new_b = model.b;
    exist_b_cone = false;
    b_cone_dim = 0;
    
    % Iterate over each block in the original model
    for i = 1:length(model.K)
        if ~strcmp(model.K{i}.type, 'b')
            new_K{end+1} = model.K{i};
            new_C{end+1} = model.C{i};
            if exist_b_cone
                new_At{end+1} = [model.At{i}, sparse(size(model.At{i},1), b_cone_dim)];
            else
                new_At{end+1} = model.At{i};
            end
        else
            % Handle 'b' type blocks by splitting into two 'l' blocks
            m = model.K{i}.size;              % Number of box-constrained variables
            l = model.K{i}.params(:,1);       % Lower bounds (m x 1)
            u = model.K{i}.params(:,2);       % Upper bounds (m x 1)
            
            exist_b_cone = true;
            b_cone_dim = m;

            % Append new constraints x1 + x2 = u - l
            new_b = [new_b - model.At{i}' * l; u - l];
            
            % Append m zero rows to all existing At blocks (excluding the new ones)
            for j = 1:(i-1)
                new_At{j} = [new_At{j}, sparse(size(new_At{j},1), m)];
            end
            
            % Split the 'b' block into x1 and x2
            new_K{end+1} = struct('type', 'b2l', 'size', m * 2);
            new_C{end+1} = [model.C{i}; zeros(m, 1)];
            A_x1 = [model.At{i}, speye(m); sparse(m, size(model.At{i},2)), speye(m)];
            new_At{end+1} = A_x1;
        end
    end
    
    % Update the model with the transformed components
    model.K = new_K';
    model.C = new_C';
    model.At = new_At';
    model.b = new_b;
end