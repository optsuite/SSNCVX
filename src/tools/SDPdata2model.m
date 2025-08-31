function model =  data2model(dataset,probname,dir_data,options)
if strcmp(dataset,'theta')
    file = fullfile(dir_data, 'SDP', probname + "_Alt.mat");
    [blk, At, C, b] = thetaread(file);
%     [b,At,cnz] = data_process(blk,At,b);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
elseif strcmp(dataset,'fap')   
    file = [dir_data, '/SDP/', probname, '.dat'];
    [blk,At,C,b,L,U] = fapread_lu_sdpnal(file);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
    model.L = L;
    model.U = U;
    model.K = fromblk(blk);
    model.C = MatCell(C);
    model.At = MatCell(At);
    model.b = b;
elseif dataset == "DIMACS"
    model_orig = load(fullfile(dir_data, dataset, probname + ".mat"));
    model_orig.C = reshape(model_orig.c, [], 1);
    model_orig.b = reshape(model_orig.b, [], 1);
    if isfield(model_orig, 'A')
        model_orig.At = model_orig.A';
    else
        model_orig.A = model_orig.At';
    end

    K{1} = BasicCone('l', model_orig.K.l);
    K{2} = BasicCone('q', model_orig.K.q);
    C{1} = model_orig.C(1 : model_orig.K.l);
    C{2} = model_orig.C(model_orig.K.l + 1 : end);
    At{1} = model_orig.At(1 : model_orig.K.l, :);
    At{2} = model_orig.At(model_orig.K.l + 1 : end, :);
    b = model_orig.b;
    model.At = MatCell(At);
    model.C = MatCell(C);
    model.b = b;
    model.K = Cone(K);
    model.name = probname;
    clear At K b C
elseif dataset == "CBLIB"
    if isfile(fullfile(dir_data, dataset, probname + ".mat"))
        load(fullfile(dir_data, dataset, probname + ".mat"));
        model = standardize(model);

        mystdmodel = struct();
        mystdmodel.At_mat = model.At_int;
        mystdmodel.b = model.b_int;
        mystdmodel.K = fromblk(model.blk_int);
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
        model.K = fromblk(model.K);
        model.name = probname;
    end
elseif dataset == "SPCA"
    C = load(fullfile(dir_data, dataset, probname + ".mat"));
    model.C = C.Sigma;
end
end