addpath('/bicmr/home/optsuite/opt/mosek/10.1/toolbox/r2017a');
setenv('PATH', [getenv('PATH') ';/bicmr/home/optsuite/opt/mosek/10.1/tools/platform/linux64x86/bin']);
Datapath = '/bicmr/home/optsuite/wry/sdp/mittelmann/';
Datapath = ['E:\seafile\Seafile\ALDSDP\code\ssnsdp-beta-2nd\sdp_data\sedumi\'];
addpath('/opt/mosek/10.1/toolbox/r2017a');
%setenv('PATH', [getenv('PATH') ';/opt/mosek/10.1/tools/platform/linux64x86/bin']);
%Datapath = '/mnt/d/mat/sdp/';
probnames = {'1dc.1024','1et.2048','1tc.2048','1zc.1024','G40_mb','G48_mb','G48mc','G55mc','G59mc','biggs', ...
            'broyden25','buck5','cancer_100','checker_1.5','chs_5000','cnhil10',  ...
			'cphil12','diamond_patch','e_moment_quadknap_17_100_0.5_2_2','e_moment_stable_17_0.5_2_2','foot',  ...
			'hand','ice_2.0','inc_1200','mater-6','neosfbr25','neosfbr30e8',  ...
			'neu1','neu1g','neu2','neu2c','neu2g','neu3','neu3g','p_auss2_3.0','rabmo','reimer5',  ...
			'ros_2000','rose15','sensor_1000','shmup4','shmup5','spar060-020-1_LS',  ...
			'swissroll','taha1a','taha1b','taha1c','theta12','theta102','theta123','tiger_texture','trto4','trto5',  ...
			'vibra4','vibra5','yalsdp','rendl1_2000_1e-6','biomedP', ...
	'prob_2_4_0','prob_2_4_1','torusg3-15','hamming_8_3_4','hamming_9_5_6'};
%probnames = {'nonc_500'};
%fid = fopen('test_Mosek_sedumi.txt', 'w');
%fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'prob', 'm', 'n', 'X_rank', 'S_rank','pinf', 'dinf', 'gap', 'iter', 'time', 'status');
%fclose(fid);
tol = 1e-7;
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = tol;
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = tol;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = tol;
for k =  59
%     if sum(k == [47 50 58 59 60]) == 1 
%         continue
%     end
    fname = probnames{k};
	fprintf('*************************************************************************************\n');
	fprintf('Solving problem %s ...\n\n', fname);
    model = load(fullfile(Datapath, [fname '.mat']));
    if ~isfield(model, 'A')
	if ~isfield(model, 'At')
            At = model.AT;
        else
	    At = model.At;
        end
            model.A = At';
    end
	fieldsToDelete = fieldnames(model.K);
	for i = 1:numel(fieldsToDelete)
		if model.K.(fieldsToDelete{i}) == 0
			model.K = rmfield(model.K, fieldsToDelete{i});
		end
	end
    prob = convert_sedumi2mosek(model.A, model.b, model.c, model.K);
	[m,n] = size(model.A);
	n=sum(prob.bardim);
    [r, res] = mosekopt('minimize info param', prob, param);
	
    X_bar = res.sol.itr.barx;
    X = cell(1, length(prob.bardim));
    S_bar = res.sol.itr.bars;
    S = cell(1, length(prob.bardim));
    startIndex = 1;
    Xvec = [res.sol.itr.xx];
    Svec = [res.sol.itr.slx];
    dot = 0;
    %fprintf('X_bar.length=%i.\n', length(X_bar));
    for i = 1:length(prob.bardim)
        dim = prob.bardim(i);
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
        dot = dot + sum(sum(X_full.*S_full));
        % 将对称矩阵存储到 X
%         spy(X_full)
        X{i} = X_full;
        S{i} = S_full;
        blk{i,1} = 's';
        blk{i,2} = size(X_full,1);
        Xvec = [Xvec ; X_full(:)];
        Svec = [Svec ; S_full(:)];
        
        startIndex = endIndex + 1;
    end
%         Xbb = svec(blk,X);
    pinf = norm(model.A*Xvec  -model.b)/(1+norm(model.b))
    y = res.sol.itr.y;
    dinf = norm(model.A'*y + Svec - model.c(:))/(1+norm(model.c))
    etaC1 = dot/(1+norm(Xvec) + norm(Svec))
    etaC2 = abs(Xvec'*Svec)/(1+abs(Xvec'*Svec))
	
	X_rank =0;
	S_rank =0;
    eigxsum = 0;
    eigssum = 0;
	fields = fieldnames(model.K);
	x_index = 1;
	for j = 1:numel(fields)
		field = fields{j};
		value = model.K.(field);
		if strcmp(field,'s')
			for m_index = 1:length(value)
				X_eig =eig(X{x_index});
                eigxsum = eigxsum + norm(X_eig(X_eig<1e-16))^2;
				S_eig =eig(S{x_index});
                eigssum = eigssum + norm(S_eig(S_eig<1e-16))^2;
				X_rank = X_rank + sum(X_eig-max(X_eig)*1e-6>0);
				S_rank = S_rank + sum(S_eig-max(S_eig)*1e-6>0);
				x_index =  x_index + 1;
			end
		elseif strcmp(field,'l')
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
    time = res.info.MSK_DINF_OPTIMIZER_TIME;
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
	out.m = m;
    out.n = n;
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
    fid = fopen('test_Mosek_sedumi.txt', 'a');
    fprintf(fid, '%s\t%d\t%d\t%d\t%d\t%.4e\t%.4e\t%.4e\t%d\t%.2f\t%s\n', fname, m, n, X_rank, S_rank, pinf, dinf, rel_gap, iter, time, status);
    fclose(fid);
end
fprintf('%i problem(s) solved by MOSEK.\n', numel(probnames));
