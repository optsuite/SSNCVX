function model =  data2model(dataset,probname,dir_data,options)
if strcmp(dataset,'theta')
    file = fullfile(dir_data, 'SDP', probname + "_Alt.mat");
    [blk, At, C, b] = thetaread(file);
%     [b,At,cnz] = data_process(blk,At,b);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
elseif strcmp(dataset,'RDM')
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
    model.K = fromblk(blk);
    model.basename = basename;
elseif strcmp(dataset,'R1TA')
    [~, basename1, basename2] = fileparts(probname);
    basename = [basename1, basename2];
    file = [dir_data, '/R1TA/content/', probname, '.mat'];
%     [b,At,cnz] = data_process(blk,At,b);
    load(file);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
elseif strcmp(dataset,'thetaplus') 
    dataset2 = 'theta';
    file = fullfile(dir_data, dataset2, probname + "_Alt.mat");
    [blk, At, C, b] = thetaread(file);
    
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
    model.L = {0};
    model.U = {inf};
elseif strcmp(dataset,'qap')  
    file = fullfile(dir_data, dataset, probname + ".dat");
    [A, B] = qapread(file);
    [blk,At,C,b,Ascale,Bscale] = qapAW_sdpnal(A, B);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
elseif strcmp(dataset,'rcp')  
    K = options.K;
     file = [dir_data,'rcp/' ,probname, '.data'];
    [blk,At,C,b,W0] = rcpread(file, K, options.n0);
     model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
    model.L = 0;
    model.U = inf;
    model.W0 = W0;
elseif strcmp(dataset,'biq')   
    file = [dir_data, '/biq_all/', probname, '.sparse'];
%      file = [dir_data,'\rcp\' ,probname, '.data'];
        Q = biqread_sdpnal(file);
    [blk,At,A,C,b,Bt,B,d] = biq_ineq_sdpnal(Q);
     model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
    model.L = {0};
    model.U = {inf};
elseif strcmp(dataset,'fap')   
    file = [dir_data, '/SDP/', probname, '.dat'];
    [blk,At,C,b,L,U] = fapread_lu_sdpnal(file);
    model.At = MatCell(At);
    model.b = b;
    model.C = MatCell(C);
    model.K = fromblk(blk);
    model.L = L;
    model.U = U;
elseif strcmp(dataset,'sedumi')    
    try
    model2 = load([dir_data,'/sedumi/',probname,'.mat']);
    if isfield(model2,'At')
    [blk,At,C,b] = read_sedumi(model2.At',model2.b,model2.c(:),model2.K);
%     [b,At,cnz] = data_process(blk,At,b);
    elseif isfield(model2,'A')
    [blk,At,C,b] = read_sedumi(model2.A,model2.b,model2.c(:),model2.K);
%     [b,At,cnz] = data_process(blk,At,b);
    elseif isfield(model2,'AT')
    [blk,At,C,b] = read_sedumi(model2.AT',model2.b,model2.c(:),model2.K);
    end
    catch
    [blk,At,C,b] = read_sdpa([dir_data,'/sedumi/',probname,'.dat-s']);
%      [b,At,cnz] = data_process(blk,At,b);
%     ttmp =zeros( size(At{1},1),1);
%     ttmp = 1.414*svec_sdpnal(blk,C);
%     At{1}= [At{1} sparse(ttmp{1})];
%     ttmpb = 2726.286738;
%     b = [b; ttmpb];
%     [blk_new,At_new,C_new, c, hists] = rdm_detect_block(blk,At,C,b);
%     At = At_new;
%     blk = blk_new;
%     C = C_new;
%     b = c;
    [b,At,cnz] = data_process(blk,At,b);
    end

%     for i = 1:length(model2.K.s)
%     blk{i,1} = 's';
%     blk{i,2} = model2.K.s;
%     C{i} = reshape(model2.c,model2.K.s(i),model2.K.s(i));
% %     for j = 1:size(model2.b,1)
% %         Acell{j} = reshape(model2.A(j,:),model2.K.s(i),model2.K.s(i));
% %     end
%     index = [];
%     for j = 1:model2.K.s(i)
%         index = [index (j-1)*model2.K.s(i)+j:(j-1)*model2.K.s(i)+model2.K.s(i)];
%     end
%     Attmp = model2.A(:,index)';
%     At{i} = Attmp;
%     end
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
elseif dataset == "sdplib"
    if isfile(fullfile(dir_data, dataset, probname + ".dat-s"))
        [blk, At, C, b] = read_sdpa(fullfile(dir_data, dataset, probname + ".dat-s"));

        model = struct();
        model.b = b;
        model.K = fromblk(blk);
        model.At = MatCell(At);
        model.C = MatCell(C);

    else
        error("No such file: %s", fullfile(dir_data, dataset, probname + ".dat-s"));
    end
elseif dataset == "SPCA"
    C = load(fullfile(dir_data, dataset, probname + ".mat"));
    model.C = C.Sigma;
end
end