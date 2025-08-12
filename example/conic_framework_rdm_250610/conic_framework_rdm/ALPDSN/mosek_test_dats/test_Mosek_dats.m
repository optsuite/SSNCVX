addpath('/bicmr/home/optsuite/opt/mosek/10.1/toolbox/r2017a');
addpath('/bicmr/home/optsuite/wry/sdp/script/theta/mexfiles_from_SDPNALplus');
setenv('PATH', [getenv('PATH') ';/bicmr/home/optsuite/opt/mosek/10.1/tools/platform/linux64x86/bin']);
Datapath = '/bicmr/home/optsuite/wry/sdp/dats/';
addpath(genpath('./'))
% Datapath = '/public/shared/sdp_data/mittelmann/';
Datapath = ['E:\seafile\Seafile\ALDSDP\code\ssnsdp-beta-2nd\sdp_data\sedumi\'];
data_dir = addpath_data();
%addpath('/opt/mosek/10.1/toolbox/r2017a');
%addpath('/mnt/d/git/gitlab/conic_programming_framework-main/ALPDSN/mexfiles_from_SDPNALplus');
%setenv('PATH', [getenv('PATH') ';/opt/mosek/10.1/tools/platform/linux64x86/bin']);
%Datapath = '/mnt/d/dats/';
%
probnames = {'CH2_1A1_STO-6GN8r14g1T2','AlH_1-Sigma+_STO-6GN14r20g1T2_5','Bex2_1_5','BH2_2A1_STO-6GN7r14g1T2','Bst_jcbpaf2','fap09',...
    'H3O+_1-A1_STO-6GN10r16g1T2_5','NH2-.1A1.STO6G.r14','NH3_1-A1_STO-6GN10r16g1T2_5','NH4+.1A1.STO6GN.r18','G60_mb','G60mc'};
% probnames = {'G60_mb','G60mc'};
fid = fopen('test_Mosek_sedumi.txt', 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'prob', 'm', 'n', 'X_rank', 'S_rank','pinf', 'dinf', 'gap', 'iter', 'time', 'status');
fclose(fid);
tol = 1e-7;
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = tol;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = tol;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = tol;
for k = [3:5 7 9 11 12]
    fname = probnames{k};
    fprintf('*************************************************************************************\n');
    fprintf('Solving problem %s ...\n\n', fname);
    [blk,At,C,b] = read_sdpa(fullfile(Datapath, [fname '.dat-s']));
    prob = convert_RDM_to_mosek(At,blk,C,b);
    [r,res] = mosekopt('minimize info param',prob,param);
    try
        model2 = load([data_dir,'/sedumi/',fname,'.mat']);
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
        [blk,At,C,b] = read_sdpa([data_dir,'/sedumi/',fname,'.dat-s']);

%         [b,At,cnz] = data_process(blk,At,b);
    end
    model.K = Cone.fromblk(blk);
    model.C = MatCell(C);

    model.At = MatCell(At);
    model.b = b;
    kk = length(model.K);
    lenk = 0;
    modeltmp = model;
    %
    model.K = {};
    para = 0;
    for ii = 1:kk
        tmpC2 = modeltmp.C{ii};
        Attmp = modeltmp.At{ii};
        cone2 = modeltmp.K{ii};
        cumsum_size = [0 cumsum(cone2.size)];
        sqcumsum_size = [0 cumsum((cone2.size.^2 + cone2.size)/2  )];
        for jj = 1 : length(modeltmp.K{ii}.size)
            cone2.size = modeltmp.K{ii}.size(jj);
            if strcmp(modeltmp.K{ii}.type,'s')
                model.C{lenk+jj,1} = tmpC2(cumsum_size(jj)+1:cumsum_size(jj+1),cumsum_size(jj)+1:cumsum_size(jj+1));
                model.At{lenk+jj,1} = Attmp(sqcumsum_size(jj)+1:sqcumsum_size(jj+1),:) ;
                model.K{lenk+jj,1} = cone2;
            else
                model.C{lenk+jj,1} = tmpC2;
                model.At{lenk+jj,1} = Attmp ;
                model.K{lenk+jj,1} = cone2;
            end
        end
        lenk = lenk + length(modeltmp.K{ii}.size);
    end



    if strcmp(modeltmp.K{1}.type,'l')
        kk = 1;

        X_bar = res.sol.itr.barx;
        X = cell(1, length(prob.bardim)+kk);
        S_bar = res.sol.itr.bars;
        S = cell(1, length(prob.bardim));
        startIndex = 1;
        %fprintf('X_bar.length=%i.\n', length(X_bar));
        X_bar = res.sol.itr.barx;
        X = cell(1, length(prob.bardim)+kk);
        S_bar = res.sol.itr.bars;
        S = cell(1, length(prob.bardim)+kk);
        X{kk} = res.sol.itr.xx;
        S{kk} = res.sol.itr.slx;
        startIndex = 1;
    else
        kk = 0;
        X_bar = res.sol.itr.barx;
        X = cell(1, length(prob.bardim)+kk);
        S_bar = res.sol.itr.bars;
        S = cell(1, length(prob.bardim));
        startIndex = 1;
        %fprintf('X_bar.length=%i.\n', length(X_bar));
        X_bar = res.sol.itr.barx;
        X = cell(1, length(prob.bardim)+kk);
        S_bar = res.sol.itr.bars;
        S = cell(1, length(prob.bardim)+kk);
        startIndex = 1;
    end
    Xvec = [res.sol.itr.xx];
    Svec = [res.sol.itr.slx];
    dot = 0;
    %fprintf('X_bar.length=%i.\n', length(X_bar));
    for i = kk+1:kk+length(prob.bardim)
        dim = prob.bardim(i-kk);
        %fprintf('dim=%i.\n', dim);
        tmp_size = dim * (dim+1)/2;
        endIndex = startIndex + tmp_size - 1;
        %fprintf('endIndex=%i.\n', endIndex);
        % 生成对称矩阵
        X_full = zeros(dim);
        X_full(tril(true(dim))) = X_bar(startIndex:endIndex);
        X_full = X_full + tril(X_full,-1)';
        S_full = zeros(dim);
        S_full(tril(true(dim))) = S_bar(startIndex:endIndex);
        S_full = S_full + S_full' - diag(diag(S_full));
        dot = dot + sum(sum(X_full*S_full));
        % 将对称矩阵存储到 X
        %         spy(X_full)
        X{i} = X_full;
        S{i} = S_full;
        %         blk{i,1} = 's';
        %         blk{i,2} = size(X_full,1);
        Xvec = [Xvec ; X_full(:)];
        Svec = [Svec ; S_full(:)];

        startIndex = endIndex + 1;
    end
    %         Xbb = svec(blk,X);
    y = res.sol.itr.y;
    %     model.A = [];
    %     for i = 1:length(At)
    %         model.A = [model.A; At{1}];
    %     end
    pinf = norm(AXfun(model.K,model.At,X)  -model.b)/(1+norm(model.b))
    y = res.sol.itr.y;
    dinf = norm(Atyfun(model.K,model.At,y) + S' - model.C)/(1+norm(model.C))
    etaC1 = dot/(1+norm(Xvec) + norm(Svec))
    etaC2 = abs(Xvec'*Svec)/(1+norm(Xvec) + norm(Svec))
    new_blk = {};
    idx = 1;



    X_rank =0;
    S_rank =0;
    eigxsum = 0;
    eigssum = 0;
    % 	field = fieldnames(model.K);
    x_index = 1;
    for j = 1:length(model.K)
        % 		field = fields{j};
        value = model.K{j}.size;
        if strcmp(model.K{j}.type,'s')
            for m_index = 1:length(value)
                X_eig =eig(X{j});
                eigxsum = eigxsum + norm(X_eig(X_eig<1e-16))^2;
                S_eig =eig(S{j});
                eigssum = eigssum + norm(S_eig(S_eig<1e-16))^2;
                X_rank = X_rank + sum(X_eig-max(X_eig)*1e-6>0);
                S_rank = S_rank + sum(S_eig-max(S_eig)*1e-6>0);
                % 				x_index =  x_index + 1;
            end
        elseif strcmp(model.K{j}.type,'l')
            X_eig = res.sol.itr.xx;
            eigxsum = eigxsum + norm(X_eig(X_eig<1e-16))^2;
            S_eig = res.sol.itr.snx;
            eigssum = eigssum + norm(S_eig(S_eig<1e-16))^2;
            X_rank = X_rank + sum(X_eig-max(X_eig)*1e-6>0);
            S_rank = S_rank + sum(S_eig-max(S_eig)*1e-6>0);
        end
    end
    etaK1 =  sqrt(eigxsum);
    etaK2 =  sqrt(eigssum);
    iter = res.info.MSK_IINF_INTPNT_ITER;
    %     time = res.info.MSK_DINF_OPTIMIZER_TIME;
    %     pinf = res.info.MSK_DINF_INTPNT_PRIMAL_FEAS;
    %     dinf = res.info.MSK_DINF_INTPNT_DUAL_FEAS;
    rel_gap = abs(res.info.MSK_DINF_INTPNT_PRIMAL_OBJ - res.info.MSK_DINF_INTPNT_DUAL_OBJ) / (1 + abs(res.info.MSK_DINF_INTPNT_PRIMAL_OBJ) + abs(res.info.MSK_DINF_INTPNT_DUAL_OBJ));
    status = 'f';
    if strcmp(res.sol.itr.solsta, 'OPTIMAL')
        status = 's';
    end
    out = struct;
    out_file = [fname '_res.mat'];
    load(strcat('./res/', out_file));
    %     out.time = time;
    out.m = length(b);
    %     out.n = n;
    out.X_rank = X_rank;
    out.S_rank = S_rank;
    out.pinf = pinf;
    out.dinf = dinf;
    out.gap = rel_gap;
    out.etaK1 = etaK1;
    out.etaK2 = etaK2;
    out.etaC1 = etaC1;
    out.etaC2 = etaC2;
    out.iter = iter;
    save(strcat('./res2/', out_file),'out');
    iter = res.info.MSK_IINF_INTPNT_ITER;
    time = out.time;
    %     pinf = res.info.MSK_DINF_INTPNT_PRIMAL_FEAS;
    %     dinf = res.info.MSK_DINF_INTPNT_DUAL_FEAS;
    %     rel_gap = abs(res.info.MSK_DINF_INTPNT_PRIMAL_OBJ - res.info.MSK_DINF_INTPNT_DUAL_OBJ) / (1 + abs(res.info.MSK_DINF_INTPNT_PRIMAL_OBJ) + abs(res.info.MSK_DINF_INTPNT_DUAL_OBJ));
    status = 'f';
    if strcmp(res.sol.itr.solsta, 'OPTIMAL')
        status = 's';
    end
    %     out.time = time;
    out_file = [fname '_res.mat'];
    save(out_file,'out');
    fid = fopen('test_Mosek_sedumi.txt', 'a');
    fprintf(fid, '%s\t%d\t%d\t%d\t%d\t%.4e\t%.4e\t%.4e\t%d\t%.2f\t%s\n', fname, prob.m, prob.n, X_rank, S_rank, pinf, dinf, rel_gap, iter, time, status);
    fclose(fid);
end
fprintf('%i problem(s) solved by MOSEK.\n', numel(probnames));
%     blktmp = blk;
%     blk={};
%     lenk = 0;
%     for ii = 1:length(blktmp)
%         tmpC2 = C{ii};
%         Attmp = At{ii};
%         [tpye cone2_size] =blktmp{ii,:};
%         cumsum_size = [0 cumsum(cone2_size)];
%         sqcumsum_size = [0 cumsum((cone2_size.^2 + cone2_size)/2  )];
%         for jj = 1 : length(cone2_size)
%             cone2_sizejj = cone2_size(jj);
%             if strcmp(blktmp{ii,1},'s')
%                 model.C{lenk+jj,1} = tmpC2(cumsum_size(jj)+1:cumsum_size(jj+1),cumsum_size(jj)+1:cumsum_size(jj+1));
%                 model.At{lenk+jj,1} = Attmp(sqcumsum_size(jj)+1:sqcumsum_size(jj+1),:) ;
%                 model.blk{lenk+jj,1} = {blktmp{ii,1}};
%                 model.blk{lenk+jj,2} = {cone2_sizejj};
%             else
%                 model.C{lenk+jj,1} = tmpC2;
%                 model.At{lenk+jj,1} = Attmp ;
%                 model.blk{lenk+jj,1} = {blktmp{ii,1}};
%                 model.blk{lenk+jj,2} = {cone2_sizejj};
%             end
%         end
%         lenk = lenk + length(cone2_size);
%     end