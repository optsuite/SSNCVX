%%% deprecated

%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-13 16:28:19
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%
function [X, out, y, S, Zb] = SSNCPplusallm(model, opts)

% SSNSDP: a semi-smooth Newton method for solving SDP.
%
%   (P) min <C,X> s.t. A(X)=b, X \in K, L <= X <= U
%   (D) min -<b,y> -<Z, L> + <Zb, U> s.t. A' * y + S + Z - Zb = C, S \in K^*, Z, Zb >= 0,
%   where K is a positive semidefinite cone and K^* is its dual cone
%
%%**********************************************************************
%

%% check input model
if ~isfield(model, 'At') || ~isfield(model, 'K') || ~isfield(model, 'C') || ~isfield(model, 'b') || ~isfield(model, 'L') || ~isfield(model, 'U')
    error("Input model is not valid");
end

%% make a copy
trans.isL = 0;
if ( ~isempty(model.L) && ~isempty(model.U) )
    trans.isL = 1;
    trans.cgtol = opts.cgtol;
    trans.cgratio = opts.cgratio;
    trans.resratio = opts.resratio;
    % else  isempty(model.L) && isempty(model.U)
    % trans.isL = 1;
end
At = model.At;
K = model.K;
b = model.b;
C = model.C;
L = model.L;
U = model.U;
%% Set the parameter
if nargin < 2; opts = []; end
opts = default_opts(opts);

if opts.log_path == ""
    fid = 1;
else
    fid = fopen(opts.log_path, 'w');
end



record = opts.record;

if record >=1
    print_std_model_info(fid, model);
end
% trans.sigma = opts.pmu;% sigma

maxits = opts.maxits;
tol = opts.tol;
tolscale = opts.tolscale;

muopts = opts.muopts;
cgopts = opts.cgopts;
NEWTopts = opts.NEWTopts;
trans.cgmin = opts.cgmin;
trans.cgmed = opts.cgmed;
trans.cgmax = opts.cgmax;
trans.sigxl = opts.sigxl;
trans.sigxm = opts.sigxm;
trans.sigxu = opts.sigxu;
trans.sigyl = opts.sigyl;
trans.sigym = opts.sigym;
trans.sigyu = opts.sigyu;
trans.sigzl = opts.sigzl ;
trans.save_path = strcat(opts.save_path,'process.txt');
str = [opts.basename, newline];
trans.sigzm =    opts.sigzm ;
trans.sigzu =    opts.sigzu ;
trans.sigql =    opts.sigql ;
trans.sigqm =   opts.sigqm;
trans.sigqu =    opts.sigqu ;
trans.cgres = 1;
trans.gfactor = opts.gfactor;

%% record original variable and preprocessing
trans.borg = b;
trans.Corg = C;
trans.Atorg = At;
trans.K = K;
trans.xLorg = L;
trans.xUorg = U;
trans.state = 1;


trans.precond = 0;
m = length(b);

trans.normborg = max(1, norm(trans.borg));
trans.normCorg = max(1, norm(trans.Corg));

[model, trans] = preprocess_SDP(model, opts, trans);
At = model.At;
K = model.K;
b = model.b;
C = model.C;
trans.Amap = @(X) AXmap(X, K, At, trans.Lchol);
trans.ATmap = @(y) Atymap(y, K, At, trans.Lchol);
trans.projectP = @(X) sdp_project_P(K, X, trans); 

trans.diagAAt = zeros(m, 1);

for p = 1:length(K)
    trans.diagAAt = trans.diagAAt + sum(At{p} .* At{p})';
end

trans.diagAAt = max(1e-4, full(trans.diagAAt));
trans.sigma = 1;

%% initial variable and operator
X.var = MatCell.zeros_like(C);
X.Avar = trans.Amap(X.var);
Zb.var = MatCell.zeros_like(C);
Zb.Avar = trans.Amap(Zb.var);
Z.var = MatCell.zeros_like(C);
Z.Avar = trans.Amap(Z.var);
W = Z.var;
y.var = zeros(trans.mdim, 1);
y.Avar = X.var;
S.var = MatCell.zeros_like(C);
S.Avar = trans.Amap(S.var);
q.var = MatCell.zeros_like(C);
q.Avar = trans.Amap(q.var);
trans.cgall = 0;
trans.bmAX = b - X.Avar;
trans.AC = trans.Amap(C);
trans.normb = max(1, norm(b));
trans.normC = max(1, norm(C));

trans = initial_NEWT_trans(trans, NEWTopts);
trans = initial_hists_trans(trans);

%%------------------------------------------------------------------
%% set up print format

out = struct;
trans.rescale = 1;
trans.ischange = 1;
trans.cgiter = 0;
doNEWT = 1;
trans.NEWT.totaliter = 0;
trans.ADMM.totaliter = 0;
retol = 1;
cstop = 0;
trans.tstart = clock;
trans.tic = tic;
for i = 1:length(K)
    trans.blk{i,1} = K{i}.type;
    trans.blk{i,2} = K{i}.size;
end
%% Compute projection
if trans.isL
    [F, S,trans] = compute_gradplus(K, X, y, Zb,q, b, C, trans);
else
    [F, S, X,trans] = compute_projection(K, X, y, b, C, trans);
end

trans = set_param(trans);
count = 0;

%% main iteration
if trans.isL
    for iter = 1:maxits


        trans.NEWT.iter = trans.NEWT.iter + 1;
        trans.NEWT.totaliter = trans.NEWT.totaliter + 1;

        if cgopts.CG_adapt 
            cgopts.CG_tol = max(min(0.1 * F.res, 0.1)*trans.cgratio, trans.cgtol)  ; % F.res
        end

        trans.NEWT.tau1 = 5*trans.NEWT.lambda*(F.FYres^trans.NEWT.sigPowy);
        trans.NEWT.tau2 = 5*trans.NEWT.lambda*(F.FZres^trans.NEWT.sigPowzb);
        trans.NEWT.tau3 = 5*trans.NEWT.lambda*(F.FXres^trans.NEWT.sigPowx);
        trans.NEWT.tau4 = 5*trans.NEWT.lambda*(F.Fqres)^trans.NEWT.sigPowq + 1e-7;


        [dy,dzb,dx,dq,trans] = gendirectionxz(F,K,At,trans,cgopts);
        trans.cgall = trans.cgall + trans.cgiter;

        ynew.var = y.var + dy;
        ynew.Avar = trans.ATmap(ynew.var);
        Zbnew.var = Zb.var + dzb;
        Xnew.var = X.var + dx;
        qnew.var = q.var + dq;

        [Fnew,Snew,trans] = compute_gradplus(K,Xnew,ynew,Zbnew,qnew,b,C,trans);

        [AX,AXorg] = trans.Amap(Xnew.var);
        Xnew.Avar = AX;
        Xnew.Avarorg = AXorg;

        %     tmpnrm = norm(ynew.Avar + Snew.var + Zbnew.var - C);

        %----------------------------------------------------------------------


        trans.NEWT.nrmd = norm(dx) + norm(dy) + norm(dzb) + norm(dq);
%         trans.NEWT.nrmd = sqrt(norm(dx)^2 + norm(dy)^2 + norm(dzb)^2 + norm(dq)^2);
        trans.NEWT.Fd = -sum(Fnew.FY .* dy);

        for k = 1:trans.nblock
            trans.NEWT.Fd = trans.NEWT.Fd - sum(sum(Fnew.FX{k} .* dx{k}));
            trans.NEWT.Fd = trans.NEWT.Fd - sum(sum(Fnew.FZ{k}.* dzb{1}));
            trans.NEWT.Fd = trans.NEWT.Fd - sum(sum(Fnew.Fq{k}.*dq{k}));
        end


        if F.res < trans.resratio
            trans.NEWT.rhs = Fnew.res * trans.NEWT.nrmd ^ 2;
        else
            trans.NEWT.rhs = trans.NEWT.nrmd ^ 2;
        end

        trans.NEWT.ratio = trans.NEWT.Fd / trans.NEWT.rhs;

        % condition on performing a newton step 解的不好的时候


        if Fnew.res < 5e1 * F.res || count <= 5 %5e1*F.res
            count = 0;
            y = ynew;
            F = Fnew;
            X = Xnew;
            Zb = Zbnew;
            S = Snew;
            q = qnew;
            if iter > 2
                trans = sigPow_CGmaxit_updateplus(F.res,trans,cgopts);
                cgopts.CG_maxit = trans.NEWT.CG_maxit;
            end
        elseif count > 5
                trans.NEWT.lambda = iter;
                count = 0;
        else
            %             trans.NEWT.tau1 = trans.NEWT.lambda*(F.FYres^trans.NEWT.sigPowy);
            %             trans.NEWT.tau2 = trans.NEWT.lambda*(F.FXres^trans.NEWT.sigPowx);
            trans.NEWT.lambda = trans.gfactor * trans.NEWT.lambda;
            count = count +1;

            %% print...   to be implemented

            continue;
        end
        %------------------------------------------------------------------
        % check optimality
        trans = record_optimal_NEWTplus(C,b,trans,X,y,S,Zb,iter);


        if (max(trans.hists.pinforg(end), trans.hists.dinforg(end) * tolscale) < tol * retol)

            if opts.tolscale < 100
                [Xnew,ynew,Snew,Zbnew,rec] =recover_var_chk2ndplus(Zb, X, y, S, trans);
                if max([rec.pinf,rec.dinf*tolscale,rec.K1dual,rec.C1]) < tol*retol
                    cstop = 1;
                    X = Xnew; y = ynew; S = Snew; Zb= Zbnew;
                    %                 norm(Zb{1}-Zbtemp{1},"fro")
                    %                 norm(Zb{1},"fro")
                    %                 norm(Zbtemp{1},"fro")
                    %                 norm(X{1}-Xtemp{1},"fro")
                    %                 norm(X{1},"fro")
                    %                 norm(Xtemp{1},"fro")
                    trans.rec = rec;
                else
                    cstop = 0;
                end
            else
                cstop = 1;
                [Xnew,ynew,Snew,Zbnew,rec] = recover_var_chk2ndplus(Zb,X,y,S,trans);
                X = Xnew; y = ynew; S = Snew;
                trans.rec = rec;
            end

        end
        log_info = {'%5s', 'iter', num2str(iter);
            '%7s', 'mu', num2str(trans.sigma, '%2.1e');
            '%8s', 'pobj', num2str(trans.hists.pobj(end), '%2.1e');
            '%8s', 'dobj', num2str(trans.hists.dobj(end), '%2.1e');
            '%7s', 'gap', num2str(trans.hists.gap(end), '%2.1e');
            '%7s', 'gaporg', num2str(trans.hists.gaporg(end), '%2.1e');
            '%7s', 'pinforg', num2str(trans.hists.pinforg(end), '%2.1e');
            '%7s', 'dinforg', num2str(trans.hists.dinforg(end), '%2.1e');
            '%7s', 'F', num2str(F.res, '%2.1e');
            '%7s', 'tau1', num2str(trans.NEWT.tau1, '%2.1e');
            '%7s', 'tau2', num2str(trans.NEWT.tau2, '%2.1e');
            '%7s', 'tau3', num2str(trans.NEWT.tau3, '%2.1e');
            '%7s', 'tau4', num2str(trans.NEWT.tau4, '%2.1e');
            '%7s', 'CGres', num2str(trans.cgres(end), '%2.1e');
                    '%7s', 'nt_res', num2str(trans.newton_res, '%2.1e');
            '%5s', 'CGit', num2str(trans.cgiter, '%5d');
            %         '%7s', 'ratio', num2str(trans.NEWT.ratio, '%2.1e')
            };

        %% define header for printing
        if  ~exist('str_head2', 'var')
            str_head2 = "";

            for i = 1:size(log_info, 1)
                str_head2 = str_head2 + sprintf(log_info{i, 1}, log_info{i, 2}) + " ";
            end
            str = [str, char(str_head2)];
            str = [str, newline];
            if record >=1
                fprintf(fid, '%s\n', str_head2);
            end
        end

        if  (cstop || ...
                iter == 1 || iter == maxits || mod(iter, NEWTopts.print_itr) == 0)

            if   mod(iter, 20 * NEWTopts.print_itr) == 0 && iter ~= maxits && ~cstop
                str = [str, char(str_head2)];
                str = [str, newline];
                if record >= 1
                    fprintf(fid, '%s\n', str_head2);
                end
            end

            str_info = "";

            for i = 1:size(log_info, 1)
                str_info = str_info + sprintf(log_info{i, 1}, log_info{i, 3}) + " ";
            end
            str = [str, char(str_info)];
            str = [str, newline];
            if record >= 1
                fprintf(fid, '%s\n', str_info);
            end

        end

        %% 更新kappa 这里用的lambdba表示
        if trans.NEWT.ratio >= NEWTopts.eta2
            trans.NEWT.lambda = max(NEWTopts.gamma1 * trans.NEWT.lambda, 1e-16); 
        elseif trans.NEWT.ratio >= NEWTopts.eta1
            trans.NEWT.lambda = NEWTopts.gamma2 * trans.NEWT.lambda;
        else
            trans.NEWT.lambda = NEWTopts.gamma3 * trans.NEWT.lambda;
        end

        if cstop
            out.status = sprintf('max(pinf,dinf) < %3.2e', tol);
            break;
        end

        if (iter == maxits) || toc(trans.tic) > 10000
            out.status = sprintf('reach the maximum iteration');
            [Xnew,ynew,Snew,Zbnew,rec] = recover_var_chk2ndplus(Zb,X,y,S,trans);
            X = Xnew; y = ynew; S = Snew; Zb= Zbnew;
            trans.rec = rec;
            break;
        end

        %     if ~doNEWT
        %         [muopts,trans] = set_param_mu(muopts,trans,opts,iter);
        %     end



        if muopts.adp_mu == 1
            if doNEWT

                pmup = trans.sigma;
                trans = mu_update_NEWTplus(trans, iter, muopts, opts);

                if trans.sigma ~= pmup
                    if record >= 1; fprintf(fid, '\n  -- mu updated: %f\n', trans.sigma); end
                    %                 [F, X, Z] = comp_res_fp_updmu(X, Z, K, b, C, trans.sigma,pmup,trans);
                end
            else
                %             trans = mu_update_ADMM(trans,iter,muopts);
            end
        end

    end
    trans.str = str;
    fid = fopen(trans.save_path,'a+');
    fprintf(fid,'%s',trans.str);
    out = generate_outs(fid, out, X, y, S, trans, opts, iter);
    out.cgall = trans.cgall;
else
    for iter = 1:maxits


        F.FXres = 0;
        for p = 1:length(F.FX); F.FXres = F.FXres + sum(sum(F.FX{p}.*F.FX{p})); end;
        F.FXres = sqrt(full(F.FXres));
        %      F.FYres =       ops_sdpnal(F.FY,'norm');
        %     F.FYres = norm(F.FY);

        F.FYres = full(sqrt(sum(sum(F.FY.*F.FY))));
        F.res = F.FXres + F.FYres;


        trans.NEWT.iter = trans.NEWT.iter + 1;
        trans.NEWT.totaliter = trans.NEWT.totaliter + 1;
        % parameter for the revised Newton system
        trans.NEWT.tau1 = trans.NEWT.lambda*(F.FYres^trans.NEWT.sigPowy);
        trans.NEWT.tau2 = trans.NEWT.lambda*(F.FXres^trans.NEWT.sigPowx);
        %         trans.NEWT.sig = 1 ;
        if cgopts.CG_adapt %更改CG的tol
            cgopts.CG_tol = max(min(0.1*F.res, 0.1)*0.1,cgopts.cgtolmin);% FZ.res
        end

        %         [d,trans] = gendirection(FZ,blk,At,trans,cgopts,NEWTopts,scale_data);
        %     [dx,dy,trans] = gendirectionxyori(F,trans.blk,At,trans,cgopts);
        [dx,dy,trans] = gendirectionxy(F, K, At, trans, cgopts, opts.method);
        trans.cgall = trans.cgall + trans.cgiter;

        %----------------------------------------------------------------------
        ynew = struct();
        ynew.var = y.var + dy;
        Xnew = struct();
        ynew.Avar = trans.ATmap(ynew.var);
        Xnew.var = X.var + dx ;


        [Fnew,S,Xnew,trans] = compute_projection(K,Xnew,ynew,b,C,trans);
        %     [Fnew,S,Xnew,trans] = compute_projection2(trans.blk,Xnew,ynew,b,C,trans);

        [AX,AXorg] = trans.Amap(Xnew.var);
        Xnew.Avar = AX;
        Xnew.Avarorg = AXorg;

        trans.NEWT.nrmd = norm(dx)+norm(dy);
        trans.NEWT.FZd = -sum(Fnew.FY .* dy);

        for k = 1:trans.nblock
            trans.NEWT.FZd = trans.NEWT.FZd - sum(sum(Fnew.FX{k} .* dx{k}));
        end
        trans.NEWT.FZdorg = -sum(F.FY.*dy);
        for k = 1:trans.nblock
            trans.NEWT.FZdorg = trans.NEWT.FZdorg - sum(sum(F.FX{k}.*dx{k}));
        end
        thetaFZdorg = trans.NEWT.FZdorg/F.res/trans.NEWT.nrmd;
        trans.NEWT.FZFZnew = sum(F.FY.*Fnew.FY);
        for k = 1:trans.nblock
            trans.NEWT.FZFZnew = trans.NEWT.FZFZnew + sum(sum(F.FX{k}.*Fnew.FX{k}));
        end
        thetaFZFZnew = trans.NEWT.FZFZnew/Fnew.res/F.res;

        if F.res < 1e-4
            trans.NEWT.rhs = Fnew.res*trans.NEWT.nrmd^2;
        else
            trans.NEWT.rhs = trans.NEWT.nrmd^2;
        end
        trans.NEWT.ratio = trans.NEWT.FZd/trans.NEWT.rhs;



        if Fnew.res < 5e1*F.res || count <= 5 %5e1*FZ.res
            count = 0;
            y = ynew;
            F = Fnew;
            X = Xnew;
            trans = sigPow_CGmaxit_update(F.res,trans,cgopts);
            cgopts.CG_maxit = trans.NEWT.CG_maxit;
        else
            %             trans.NEWT.tau1 = trans.NEWT.lambda*(FZ.FYres^trans.NEWT.sigPowy);
            %             trans.NEWT.tau2 = trans.NEWT.lambda*(FZ.FXres^trans.NEWT.sigPowx);
            trans.NEWT.lambda = trans.gfactor* trans.NEWT.lambda;
            count = count +1;
            continue;
        end

        %------------------------------------------------------------------
        % check optimality
        trans = record_optimal_NEWT(C,b,trans,X,y,S,iter);

        if (max(trans.hists.pinforg(end),trans.hists.dinforg(end)*tolscale) < tol*retol)
            if opts.tolscale < 100
                [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,F,X,y,S,W,trans);
                if max([rec.pinf,rec.dinf*tolscale,rec.Knorm]) < tol*retol
                    cstop = 1;
                    X = Xnew; y = ynew; S = Snew;
                    trans.rec = rec;
                else
                    cstop = 0;
                end
            else
                cstop = 1;
                [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,F,X,y,S,W,trans);
                X = Xnew; y = ynew; S = Snew;
                trans.rec = rec;
            end
        end


        log_info = {'%5s', 'iter', num2str(iter);
            '%7s', 'mu', num2str(trans.sigma, '%2.1e');
            '%8s', 'pobj', num2str(trans.hists.pobj(end), '%2.1e');
            '%8s', 'dobj', num2str(trans.hists.dobj(end), '%2.1e');
            '%7s', 'gap', num2str(trans.hists.gap(end), '%2.1e');
            '%7s', 'gaporg', num2str(trans.hists.gaporg(end), '%2.1e');
            '%7s', 'pinforg', num2str(trans.hists.pinforg(end), '%2.1e');
            '%7s', 'dinforg', num2str(trans.hists.dinforg(end), '%2.1e');
            '%7s', 'F', num2str(F.res, '%2.1e');
            '%7s', 'tau1', num2str(trans.NEWT.tau1, '%2.1e');
            '%7s', 'tau2', num2str(trans.NEWT.tau2, '%2.1e');
            '%7s', 'CGres', num2str(trans.cgres(end), '%2.1e');
            %         '%7s', 'nt_res', num2str(trans.newton_res, '%2.1e');
            '%5s', 'CGit', num2str(trans.cgiter, '%5d');
            % '%7s', 'ratio', num2str(trans.NEWT.ratio, '%2.1e')
            };

        %% define header for printing
        if  ~exist('str_head2', 'var')
            str_head2 = "";

            for i = 1:size(log_info, 1)
                str_head2 = str_head2 + sprintf(log_info{i, 1}, log_info{i, 2}) + " ";
            end
            str = [str, char(str_head2)];
            str = [str, newline];
            if record >=1
                fprintf(fid, '%s\n', str_head2);
            end
        end

        %% print iteration information
        if  (cstop || ...
                iter == 1 || iter == maxits || mod(iter, NEWTopts.print_itr) == 0)

            if   mod(iter, 20 * NEWTopts.print_itr) == 0 && iter ~= maxits && ~cstop
                str = [str, char(str_head2)];
                str = [str, newline];
                if record >= 1
                    fprintf(fid, '%s\n', str_head2);
                end
            end

            str_info = "";

            for i = 1:size(log_info, 1)
                str_info = str_info + sprintf(log_info{i, 1}, log_info{i, 3}) + " ";
            end
            str = [str, char(str_info)];
            str = [str, newline];
            if record >= 1
                fprintf(fid, '%s\n', str_info);
            end

        end
        %% 更新kappa 这里用的lambdba表示
        if trans.NEWT.ratio >= NEWTopts.eta2
            trans.NEWT.lambda = max(NEWTopts.gamma1*trans.NEWT.lambda, 1e-16); % ratio越大说明越好，大的话要减小kappa
        elseif trans.NEWT.ratio >= NEWTopts.eta1
            trans.NEWT.lambda = NEWTopts.gamma2*trans.NEWT.lambda;
        else
            trans.NEWT.lambda = NEWTopts.gamma3*trans.NEWT.lambda;
        end

        if cstop
            out.status = sprintf('max(pinf,dinf) < %3.2e',tol);
            break;
        end
        if (iter == maxits) || toc(trans.tic) > 20000
            out.status = sprintf('reach the maximum iteration');
            [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,F,X,y,S,W,trans);
            X = Xnew; y = ynew; S = Snew;
            trans.rec = rec;
            break;
        end

        %     if ~doNEWT
        %         [muopts,trans] = set_param_mu(muopts,trans,opts,iter);
        %     end


        if opts.tolscale >= 100
            if trans.NEWT.totaliter == 500
                retol = 3;
                tolscale = tolscale/2;
                if record ==1;fprintf('change the tol\n'); end
            elseif trans.NEWT.totaliter == 750
                retol = 6;
                tolscale = tolscale/2;
                if record ==1;fprintf('change the tol\n'); end
            end
        end

        if muopts.adp_mu == 1
            pmup = trans.sigma;
            trans = mu_update_NEWT(trans,iter,muopts,opts);

            if trans.sigma ~= pmup
                if record >=1; fprintf('  -- mu updated: %f\n',trans.sigma); end
                %                 [FZ, X, Z] = comp_res_fp_updmu(X, Z, blk, b, C, trans.sigma,pmup,trans);
            end
            %             trans = mu_update_ADMM(trans,iter,muopts);
        end

    end
    trans.str = str;
    fid = fopen(trans.save_path,'a+');
    fprintf(fid,'%s',trans.str);
    out = generate_outs(fid, out, X, y, S, trans, opts, iter);
    out.cgall = trans.cgall;
end


end



%% End main routine

function opts = default_opts(opts)
%% ------------------------------------------------------------------
if ~isfield(opts, 'log_path'); opts.log_path = ""; end
if ~isfield(opts, 'record'); opts.record = 1; end
if ~isfield(opts, 'scale_data'); opts.scale_data = 1; end
if ~isfield(opts, 'pmu'); opts.pmu = 1; end
if ~isfield(opts, 'sstruct'); opts.sstruct = 1; end
if ~isfield(opts, 'maxits'); opts.maxits = 1000; end
if ~isfield(opts, 'tol'); opts.tol = 1e-6; end
if ~isfield(opts, 'tolscale'); opts.tolscale = 1; end

%% -------------------------------------------------------------------------
% options for fixed-point algorithms
if ~isfield(opts, 'ADMMopts'); opts.ADMMopts = struct; end
ADMMopts = opts.ADMMopts;
if ~isfield(ADMMopts, 'print_itr'); ADMMopts.print_itr = 20; end
if ~isfield(ADMMopts, 'rho'); ADMMopts.rho = 1.618; end
if ~isfield(ADMMopts, 'imaxit'); ADMMopts.imaxit = 20000; end
opts.ADMMopts = ADMMopts;

%% ------------------------------------------------------------------------
% options for the semismooth Newton algorithm
if ~isfield(opts, 'NEWTopts'); opts.NEWTopts = struct; end
NEWTopts = opts.NEWTopts;
if ~isfield(NEWTopts,'sigPowx');   NEWTopts.sigPowx = 1;       end
if ~isfield(NEWTopts,'sigPowy');   NEWTopts.sigPowy = 1;       end
if ~isfield(NEWTopts,'sigPowz');   NEWTopts.sigPowz = 1;       end
if ~isfield(NEWTopts,'sigPowq');   NEWTopts.sigPowq = 1;       end
if ~isfield(NEWTopts, 'resFac'); NEWTopts.resFac = 0.98; end
if ~isfield(NEWTopts, 'eta1'); NEWTopts.eta1 = 1e-4; end
if ~isfield(NEWTopts, 'eta2'); NEWTopts.eta2 = 0.9; end
if ~isfield(NEWTopts, 'gamma1'); NEWTopts.gamma1 = opts.gamma1; end %可调
if ~isfield(NEWTopts, 'gamma2'); NEWTopts.gamma2 = opts.gamma2; end %可调
if ~isfield(NEWTopts, 'gamma3'); NEWTopts.gamma3 = opts.gamma3; end %可调
if ~isfield(NEWTopts, 'lambda'); NEWTopts.lambda = 1; end %adjust mu
if ~isfield(NEWTopts, 'print_itr'); NEWTopts.print_itr = 1; end
opts.NEWTopts = NEWTopts;

%% ------------------------------------------------------------------------
% parameter for newton system
if ~isfield(opts, 'method'); opts.method = 'iterative'; end


%% ------------------------------------------------------------------------
% parameter for CG
if ~isfield(opts, 'cgopts'); opts.cgopts = struct; end
cgopts = opts.cgopts;
if ~isfield(cgopts, 'CG_maxit'); cgopts.CG_maxit = 500; end %可调
if ~isfield(cgopts, 'CG_tol'); cgopts.CG_tol = 1e-2; end
if ~isfield(cgopts, 'cgtolmin'); cgopts.cgtolmin = opts.cgtol; end
if ~isfield(cgopts, 'CG_adapt'); cgopts.CG_adapt = 1; end
opts.cgopts = cgopts;

%% ------------------------------------------------------------------
% parameters for adjusting pmu
if ~isfield(opts, 'muopts'); opts.muopts = struct; end
muopts = opts.muopts;
if ~isfield(muopts, 'adp_mu'); muopts.adp_mu = 1; end %1 or 0
if ~isfield(muopts, 'NEWT'); muopts.NEWT = struct; end
if ~isfield(muopts.NEWT, 'mu_min'); muopts.NEWT.mu_min = 1e-6; end
if ~isfield(muopts.NEWT, 'mu_max'); muopts.NEWT.mu_max = 1e6; end
if ~isfield(muopts.NEWT, 'mu_update_itr'); muopts.NEWT.mu_update_itr = opts.muopts.mu_update_itr; end %10 可调
if ~isfield(muopts.NEWT, 'mu_delta'); muopts.NEWT.mu_delta = 5e-1; end % 可调
if ~isfield(muopts.NEWT, 'mu_fact'); muopts.NEWT.mu_fact = opts.muopts.mu_fact; end %5/3 可调


opts.muopts = muopts;
end


function [model_new,trans] = preprocess_SDP(model,opts,trans)

K = model.K;
At = model.At;
b = model.b;
C = model.C;


%% scale the data and cholesky decomposition
trans.mdim = length(b);
trans.nblock = length(K);
DA = speye(trans.mdim);
bscale = 1;
Cscale = 1;

if opts.scale_data
    % normA: norm of each row of A
    normA = zeros(trans.mdim,1);
    for k = 1:(trans.nblock)
        normA = normA + sum(At{k}.*At{k})';
    end
    normA = max(1,sqrt(normA));
    DA = spdiag(1./normA);
end
for k = 1:trans.nblock
    At{k} = At{k}*DA;
end
b = DA*b;
scale.DA = DA;

AAt = AAtfun(At);
m = size(AAt,1);
Lchol = struct;
if (nnz(AAt) < 0.2*m*m); use_spchol=1; else,  use_spchol=0; end
%use_spchol=0;
if (use_spchol)
    [Lchol.R, Lchol.p, Lchol.perm] = chol(AAt,'vector');
    Lchol.Rt = Lchol.R';
    Lchol.matfct_options = 'spcholmatlab';
else
    if issparse(AAt); AAt = full(AAt); end;
    Lchol.matfct_options = 'chol';
    Lchol.perm = 1:m;
    [Lchol.R, indef] = chol(AAt);
    Lchol.Rt = Lchol.R';
end
if trans.isL

end
% if isequal(Lchol.R ,eye(size(Lchol.R)),'fro') < 1e-10
if isdiag(Lchol.R) && (norm(diag(Lchol.R) - 1) < 1e-7)
    Lchol.isidentity = true;
else
    Lchol.isidentity = false;
end
b = fwsolve(Lchol,b);

if (opts.scale_data==1)
    bscale = max(1,norm(b));
    Cscale = max(1,norm(C));
end
trans.xL = trans.xLorg/bscale;
trans.xU = trans.xUorg/bscale;
b = b / bscale; %单位化
C = 1 / Cscale * C;

model_new.K = K;
model_new.At = At;
model_new.b = b;
model_new.C = C;

objscale = bscale*Cscale;
scale.bscale = bscale;
scale.Cscale = Cscale;
scale.objscale = objscale;
scale.scale_data = opts.scale_data;
trans.scale = scale;
trans.Lchol = Lchol;
end


function trans = initial_NEWT_trans(trans, NEWTopts)
trans.NEWT.lambda = NEWTopts.lambda;
trans.NEWT.sigPowx = trans.sigxu;
trans.NEWT.sigPowy = trans.sigyu;
trans.NEWT.sigPowzb = trans.sigzu;
trans.NEWT.sigPowq = trans.sigqu;
trans.NEWT.iter = 0;
trans.NEWT.swt = 0;
trans.NEWT.CG_maxit = trans.cgmin;
trans.NEWT.subiter = 0;
trans.NEWT.maxiter = 20;
trans.NEWT.lastres = inf;
end

function trans = initial_hists_trans(trans)
trans.hists.maxinf = inf;
end

function [AX, AXorg] = AXmap(X, K, At, Lchol)
AXorg = AXfun(K,At,X);
AX = fwsolve(Lchol, AXorg);
end

function Aty = Atymap(y, K, At, Lchol)

Aty = Atyfun(K, At, bwsolve(Lchol, y));
end



function [F,S,trans] = compute_gradplus(K,X,y,Zb,q,b,C,trans)
W = trans.ATmap(y.var) +  X.var / trans.sigma  - C + Zb.var;
[S0,F.par] = project(K,W);
F.par.precond = trans.precond; %加入预条件，为pre做准备

tmp2 = q.var   - Zb.var * trans.sigma;
F.Dh = MatCell.zeros_like(C);

for k = 1:trans.nblock
F.Dh{k} = double(tmp2{k}> trans.xL{k} & tmp2{k} < trans.xU{k} );
end
tmp2 = trans.projectP(tmp2)/trans.sigma;
% Zb.var = Zb.var - q.var/trans.sigma + tmp2; 
F.par.Dsch12 = F.par.Dsch2;
S.var = S0 - W;
S.Avar = trans.Amap(S.var);
F.FY = -b + trans.sigma*trans.Amap(S0);
F.FZ = trans.sigma * S0 - trans.sigma * tmp2;
F.FX = X.var / trans.sigma   - S0;
F.Fq = q.var / trans.sigma   - tmp2;
F.FYres = norm(F.FY);
F.FZres = norm(F.FZ);
F.FXres = norm(F.FX);
F.Fqres = norm(F.Fq);
F.res = F.FYres + F.FXres + F.FZres + F.Fqres;

% F.res = sqrt(F.FYres^2+F.FZres^2+F.FXres^2+F.Fqres^2);

end

function [F, S, X, trans] = compute_projection(K, X, y, b, C, trans)
if trans.nblock >= 2
    X.var = project(K, X.var);
    % project2_sdpnal()
end
tmp = C - X.var / trans.sigma ;
W = trans.ATmap(y.var) - tmp;

[S0, F.par] = project(K, W);


F.par.Dsch12 = F.par.Dsch2;
% F.par.shift = {0};
F.par.precond = trans.precond; %加入预条件，为pre做准备
% F.par.P1t{1} = F.par.P1{1}';
S.var = S0 - W;
S.Avar = trans.Amap(S.var);
F.FY = -b + trans.sigma * trans.Amap(S0);
F.FX = X.var / trans.sigma  - S0;
F.res = norm(F.FY) + norm(F.FX);
end


function trans = record_optimal_NEWTplus(C, b, trans, X, y, S, Zb, iter)
% check optimality
% bTy  = full(b'*(trans.AC + (Z.Avar-X.Avar)/trans.sigma));
bTy = full(b' * y.var);
trCX = full(dot(C, X.var));
dobj = trans.scale.objscale * bTy;
pobj = trans.scale.objscale * trCX;
gap = abs(trCX - bTy) / max(1, abs(trCX));
gaporg = abs(pobj - dobj) / max(1, abs(pobj));

bmAX = b - X.Avar;
bmAXnrm = norm(X.Avarorg ./ diag(trans.scale.DA) * trans.scale.bscale - trans.borg);

pinforg = bmAXnrm / (1 + norm(trans.borg));
pinf = norm(bmAX) / trans.normb;


y = y.var;

y = trans.scale.Cscale * (trans.scale.DA * bwsolve(trans.Lchol, y));
trans.AtymCSnrm = norm(Atyfun(trans.K, trans.Atorg,y) + trans.scale.Cscale * S.var - trans.Corg + trans.scale.Cscale * Zb.var);
dinf = trans.AtymCSnrm/(1+trans.normC);
dinforg = trans.AtymCSnrm/(1+norm(trans.Corg));
PX = norm(trans.projectP(X.var) - X.var);
PZb = norm(trans.projectP(- Zb.var));
etaP1 = PX/trans.normC;
etaP1dual = PZb/trans.normC;
etaP1org = PX/(1+norm(trans.Corg));
etaP1dualorg = PZb/(1+norm(trans.Corg));

trans.hists.pobj(iter) = pobj;
trans.hists.dobj(iter) = dobj;
trans.hists.gaporg(iter) = gaporg;
trans.hists.gap(iter) = gap;
trans.hists.pinf(iter) = pinf;
trans.hists.pinforg(iter) = pinforg;
trans.hists.dinf(iter) = dinf;
trans.hists.dinforg(iter) = dinforg;
trans.hists.pvd(iter) = pinf / dinf;
trans.hists.dvp(iter) = dinf / pinf;
trans.hists.pvdorg(iter) = pinforg / dinf;
trans.hists.dvporg(iter) = dinf / pinforg;
trans.hists.cgiter(iter) = trans.cgiter;
trans.hists.P1 = etaP1;
trans.hists.P1org = etaP1org;
trans.hists.P1dual = etaP1dual;
trans.hists.P1dualorg = etaP1dualorg;
maxinf = max(pinf, dinf);

trans.hists.isNEWT(iter) = 1;
trans.hists.maxinf = maxinf;

trans.bmAX = bmAX;
end

function trans = record_optimal_NEWT(C, b, trans, X, y, S, iter)
% check optimality
% bTy  = full(b'*(trans.AC + (Z.Avar-X.Avar)/trans.sigma));
bTy = full(b' * y.var);
trCX = full(dot(C, X.var));
dobj = trans.scale.objscale * bTy;
pobj = trans.scale.objscale * trCX;
gap = abs(trCX - bTy) / max(1, abs(trCX));
gaporg = abs(pobj - dobj) / max(1, abs(pobj));

bmAX = b - X.Avar;
bmAXnrm = norm(X.Avarorg ./ diag(trans.scale.DA) * trans.scale.bscale - trans.borg);
% pinforg = bmAXnrm/trans.normborg;
pinforg = bmAXnrm / (1 + norm(trans.borg));
pinf = norm(bmAX) / trans.normb;



y = y.var;
% S = S.var;
% scale = trans.scale;
y = trans.scale.Cscale * (trans.scale.DA * bwsolve(trans.Lchol, y));
% trans.AtymCSnrm = F.res/trans.sigma;
trans.AtymCSnrm = norm((Atyfun(trans.K, trans.Atorg, y) + trans.scale.Cscale * S.var) - trans.Corg);
dinf = trans.AtymCSnrm / trans.normC;
dinforg = trans.AtymCSnrm / (1 + norm(trans.Corg));

trans.hists.pobj(iter) = pobj;
trans.hists.dobj(iter) = dobj;
trans.hists.gaporg(iter) = gaporg;
trans.hists.gap(iter) = gap;
trans.hists.pinf(iter) = pinf;
trans.hists.pinforg(iter) = pinforg;
trans.hists.dinf(iter) = dinf;
trans.hists.dinforg(iter) = dinforg;
trans.hists.pvd(iter) = pinf / dinf;
trans.hists.dvp(iter) = dinf / pinf;
trans.hists.pvdorg(iter) = pinforg / dinf;
trans.hists.dvporg(iter) = dinf / pinforg;
trans.hists.cgiter(iter) = trans.cgiter;

maxinf = max(pinf, dinf);

trans.hists.isNEWT(iter) = 1;
trans.hists.maxinf = maxinf;

trans.bmAX = bmAX;
end

function [X, y, S, Z,rec] = recover_var_chk2ndplus(Z, X, y, S, trans)


y = y.var;


X = X.var;
S = S.var;
Z = Z.var;
% y = y.var;
%     end

scale = trans.scale;
X = scale.bscale * X;
y = scale.Cscale * (scale.DA * bwsolve(trans.Lchol, y));
S = scale.Cscale * S;
Z = scale.Cscale * Z;
XZnorm = norm(X)+norm(Z);
%% compute optimal index
AX = AXfun(trans.K, trans.Atorg, X);
bmAXnrm = norm(AX - trans.borg);
rec.pinf = bmAXnrm / (1 + norm(trans.borg));
AtymCSnrm = norm(Atyfun(trans.K, trans.Atorg, y) + S - trans.Corg + Z);

rec.dinf = AtymCSnrm / (1 + norm(trans.Corg));
pobj = full(dot(trans.Corg, X));


if norm(trans.xLorg) == 0
dobj = full(trans.borg'*y);
PX = norm(trans.projectP(X) - X) ;
PZb = abs(dot(X, Z));
else
mZ = -Z;
mZ= mZ{1};
mZ(abs(mZ)<1e-5) = 0;
mZP = max(mZ,0);
mZN = min(mZ,0);
trans.xUorg = min(trans.xUorg{1},10000);
trans.xLorg = max(trans.xLorg{1},-10000);

dobj = full(trans.borg'*y) - sum(mZP.*trans.xUorg,'all') - sum(mZN.*trans.xLorg,'all') ;

mZ = Z;
mZ= mZ{1};
mZ(abs(mZ)<1e-5) = 0;
mZP = max(mZ,0);
mZN = min(mZ,0);

zact = double(mZP>0).*trans.xLorg + double(mZN<0).*trans.xUorg;
XmPX = X{1} - zact;
% PX = ops_sdpnal(ops_sdpnal(trans.projectP(X),'-',X),'norm');
PX = ops_sdpnal(ops_sdpnal(min(max(X{1},trans.xLorg),trans.xUorg),'-',X),'norm');
% PZb = ops_sdpnal(trans.projectP(ops_sdpnal(-1,'*',Zb)),'norm');

% PZb = abs(ops_sdpnal(X,'inprod',Zb));
PZb = abs(sum(Z{1}.*XmPX,"all"));
end

etaP1org = PX/(1+XZnorm );
etaP1dualorg = PZb/(1+XZnorm);

rec.gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
rec.pobj = pobj;
rec.dobj = dobj;
rec.totaltime = etime(clock, trans.tstart);
trXS = dot(X, S);
trXZb = dot(X, Z);

normX = norm(X);
normS = norm(S);
normZ = norm(Z);
[~, Xpar] = project(trans.K, X);
[~, Spar] = project(trans.K, S);
eigXnorm = 0;
eigSnorm = 0;

for p = 1:length(trans.K)
    cone = trans.K{p};
    n = sum(cone.size);

    if strcmp(cone.type, 's')
        rx = length(Xpar.posidx{p});
        rs = length(Spar.posidx{p});
        xnegidx = (rx + 1):n;
        snegidx = (rs + 1):n;
    elseif strcmp(cone.type, 'q')
        xnegidx = setdiff([1:2 * length(cone.size)], Xpar.posidx{p});
        snegidx = setdiff([1:2 * length(cone.size)], Spar.posidx{p});
    elseif strcmp(cone.type, 'l')
        xnegidx = setdiff([1:n], Xpar.posidx{p});
        snegidx = setdiff([1:n], Spar.posidx{p});
    elseif strcmp(cone.type, 'u')
        % pass
    end

    eigXnorm = eigXnorm + norm(Xpar.dd{p}(xnegidx), 'fro') ^ 2;
    eigSnorm = eigSnorm + norm(Spar.dd{p}(snegidx), 'fro') ^ 2;
end

eigXnorm = sqrt(eigXnorm);
eigSnorm = sqrt(eigSnorm);
rec.K1 = eigXnorm / (1 + normX);
rec.K1dual = eigSnorm / (1 + normS);
rec.C1 = abs(trXS) / (1 + normX + normS);
rec.C2 = abs(trXZb)/(1+normX+normZ);
rec.P1 = etaP1org;
rec.P1dual = etaP1dualorg;
end

function [X, y, S, rec] = recover_var_chk2nd(Z, F, X, y, S, W, trans)


y = y.var;

X = X.var;
S = S.var;
%     end

scale = trans.scale;
X = scale.bscale * X;
y = scale.Cscale * (scale.DA * bwsolve(trans.Lchol, y));
S = scale.Cscale * S;

%% compute optimal index
AX = AXfun(trans.K, trans.Atorg, X);
bmAXnrm = norm(AX - trans.borg);
rec.pinf = bmAXnrm / (1 + norm(trans.borg));
AtymCSnrm = norm(Atyfun(trans.K, trans.Atorg, y) + S - trans.Corg);

rec.dinf = AtymCSnrm / (1 + norm(trans.Corg));
pobj = full(dot(trans.Corg, X));
dobj = full(trans.borg' * y);
rec.gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
rec.pobj = pobj;
rec.dobj = dobj;
rec.totaltime = etime(clock, trans.tstart);
trXS = dot(X, S);

normX = norm(X);
normS = norm(S);
[~, Xpar] = project(trans.K, X);
[~, Spar] = project(trans.K, S);
[Xtmp,~] = project(trans.K,  X - S);
Knorm = norm(Xtmp - X);
eigXnorm = 0;
eigSnorm = 0;

for p = 1:length(trans.K)
    cone = trans.K{p};
    n = sum(cone.size);

    if strcmp(cone.type, 's')
        rx = length(Xpar.posidx{p});
        rs = length(Spar.posidx{p});
        xnegidx = (rx + 1):n;
        snegidx = (rs + 1):n;
    elseif strcmp(cone.type, 'q')
        xnegidx = setdiff([1:2 * length(cone.size)], Xpar.posidx{p});
        snegidx = setdiff([1:2 * length(cone.size)], Spar.posidx{p});
    elseif strcmp(cone.type, 'l')
        xnegidx = setdiff([1:n], Xpar.posidx{p});
        snegidx = setdiff([1:n], Spar.posidx{p});
    elseif strcmp(cone.type, 'u')
        % pass
    end

    eigXnorm = eigXnorm + norm(Xpar.dd{p}(xnegidx), 'fro') ^ 2;
    eigSnorm = eigSnorm + norm(Spar.dd{p}(snegidx), 'fro') ^ 2;
end

eigXnorm = sqrt(eigXnorm);
eigSnorm = sqrt(eigSnorm);
rec.K1 = eigXnorm / (1 + normX);
rec.K1dual = eigSnorm / (1 + normS);
rec.C1 = abs(trXS) / (1 + normX + normS);
rec.Knorm = Knorm /(1 + normX + normS);

end


function out = generate_outs(fid, out, X, y, S, trans, opts, iter)

rec = trans.rec;

if opts.record
    fprintf(fid, '\n-------------------------------------------------------------------------------------------------\n') ;
    fprintf(fid, '%14s %14s %10s %10s %10s %10s %10s %12s\n', 'pboj', 'dobj', 'gap', 'pinf', 'dinf', 'C1', 'K1', 'K1dual');
    fprintf(fid, '%14.8e %14.8e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n', ...
        rec.pobj, rec.dobj, rec.gap, rec.pinf, rec.dinf, rec.C1, rec.K1, rec.K1dual) ;
    fprintf(fid, '\n-------------------------------------------------------------------------------------------------\n') ;
end

trans.hists.pobj(iter + 2) = rec.pobj;
trans.hists.dobj(iter + 2) = rec.dobj;
trans.hists.pinf(iter + 2) = rec.pinf;
trans.hists.pinforg(iter + 2) = rec.pinf;
trans.hists.dinf(iter + 2) = rec.dinf;
trans.hists.dinforg(iter + 2) = rec.dinf;
trans.hists.gap(iter + 2) = rec.gap;
trans.hists.gaporg(iter + 2) = rec.gap;
trans.hists.pvd(iter + 2) = rec.pinf / rec.dinf;
trans.hists.dvp(iter + 2) = rec.dinf / rec.pinf;
trans.hists.pvdorg(iter + 2) = rec.pinf / rec.dinf;
trans.hists.dvporg(iter + 2) = rec.dinf / rec.pinf;
%trans.hists.res(iter+2) = res;

out.pobj = rec.pobj;
out.dobj = rec.dobj;
out.gap  = rec.gap;
out.iter  = iter;
out.C1 = rec.C1;

out.K1 = rec.K1;
out.K1dual = rec.K1dual;

if trans.isL
    out.C2 = rec.C2;
    out.P1 = rec.P1;
    out.P1dual = rec.P1dual;
else
    out.Knorm = rec.Knorm;
end
out.pinf = rec.pinf;
out.dinf = rec.dinf;
out.rec = rec;
out.hists = trans.hists;
out.totaltime = rec.totaltime;

end

%判断是做ADMM还是牛顿步

function trans = set_param(trans)
trans.NEWT.iter = 0;
trans.NEWT.swt = 0;

trans.NEWT.subiter = trans.NEWT.subiter + 1;
trans.NEWT.maxiter = 10;
trans.NEWT.maxiter = 400;

end

function trans = sigPow_CGmaxit_update(res, trans,cgopts)

if res < 1e-4
    sigPowy = trans.sigyl;
    sigPowx = trans.sigxl;
else
    sigPowy = trans.sigyu;
    sigPowx = trans.sigxu;
end

% CG_maxit = trans.cgmax;
% if trans.state == 1  && min(trans.cgres(end),trans.cgres(end-1))>cgopts.CG_tol && trans.cgiter > trans.cgmin
%         trans.NEWT.CG_maxit = trans.cgmed;
%         trans.state = 2;
% elseif trans.state == 2  && trans.cgres(end)>max(cgopts.CG_tol,1e-6) && trans.cgiter > trans.cgmed
%         trans.NEWT.CG_maxit = trans.cgmax;
%         trans.state = 3;
% elseif trans.state == 3  && trans.cgres(end)>max(cgopts.CG_tol,1e-6) && trans.cgiter > trans.cgmax
%     trans.NEWT.CG_maxit = trans.cgmin;
%     trans.state = 1;
% end

% if trans.iter > 300
%     sigPowy      = 0.8;
%     sigPowx      = 0.3;
%     trans.NEWT.CG_maxit = 500;
% end

if trans.NEWT.CG_maxit == trans.cgmin  && min(trans.cgres(end),trans.cgres(end-1))>cgopts.CG_tol && trans.cgiter > trans.cgmin
    trans.NEWT.CG_maxit = trans.cgmax;
elseif trans.NEWT.CG_maxit == trans.cgmax  && trans.cgres(end)>max(cgopts.CG_tol,1e-6) && trans.cgiter > trans.cgmax
    trans.NEWT.CG_maxit = trans.cgmin;
end

trans.NEWT.sigPowx = sigPowx;
trans.NEWT.sigPowy = sigPowy;
end

function trans = sigPow_CGmaxit_updateplus(res,trans,cgopts)
if res < 1e-4
    sigPowy      = trans.sigyl;
    sigPowx      =trans.sigxl;
    sigPowzb     = trans.sigzl;
    sigPowq     = trans.sigql;
elseif res < 1e-2
    sigPowy      = trans.sigym;
    sigPowx      = trans.sigxm;
    sigPowzb     = trans.sigzm;
    sigPowq     = trans.sigqm;
else
    sigPowy      = trans.sigyu;
    sigPowx      = trans.sigxu;
    sigPowzb     = trans.sigzu;
    sigPowq     = trans.sigqu;
end



% CG_maxit = trans.cgmax;
if trans.NEWT.CG_maxit == trans.cgmin  && min(trans.cgres(end),trans.cgres(end-1))>cgopts.CG_tol && trans.cgiter > trans.cgmin
    trans.NEWT.CG_maxit = trans.cgmed;
elseif trans.NEWT.CG_maxit == trans.cgmed  && trans.cgres(end)>max(cgopts.CG_tol,1e-6) && trans.cgiter > trans.cgmed
    trans.NEWT.CG_maxit = trans.cgmax;
elseif trans.NEWT.CG_maxit == trans.cgmax  && trans.cgres(end)>max(cgopts.CG_tol,1e-6) && trans.cgiter > trans.cgmax
    trans.NEWT.CG_maxit = trans.cgmin;
end

trans.NEWT.sigPowx = sigPowx;
trans.NEWT.sigPowy = sigPowy;
trans.NEWT.sigPowzb = sigPowzb;
trans.NEWT.sigPowq = sigPowq;
% trans.NEWT.CG_maxit = CG_maxit;
% trans.NEWT.CG_maxit = CG_maxit;
end


function trans = mu_update_NEWT(trans, iter, muopts, opts)
smean = @geo_mean;
muNEWT = muopts.NEWT;

if mod(iter, muNEWT.mu_update_itr) == 0
    sitr = iter - muNEWT.mu_update_itr + 1;
    avg_pvd = smean(trans.hists.pvd(sitr:iter));
    avg_dvp = smean(trans.hists.dvp(sitr:iter));
    if iter < 300
        avg_dvp = smean(trans.hists.dvp(sitr:iter));
    else
        avg_dvp = smean(trans.hists.dvporg(sitr:iter));
        muNEWT.mu_fact = 1.2;
        muNEWT.mu_delta = 5;
    end
    if iter == 300
        trans.sigma = 1;
    end
    if avg_dvp > muNEWT.mu_delta
        trans.sigma = trans.sigma * (muNEWT.mu_fact);
    else
        trans.sigma = trans.sigma / (muNEWT.mu_fact);
    end

    trans.sigma = min(muNEWT.mu_max, max(muNEWT.mu_min, trans.sigma));
end

end


function trans = mu_update_NEWTplus(trans, iter, muopts, opts)
smean = @geo_mean;
muNEWT = muopts.NEWT;

if mod(iter, muNEWT.mu_update_itr) == 0
    sitr = iter - muNEWT.mu_update_itr + 1;
    avg_pvd = smean(trans.hists.pvd(sitr:iter));
    avg_dvp = smean(trans.hists.dvp(sitr:iter));
%     if iter < 300
%         avg_dvp = smean(trans.hists.dvp(sitr:iter));
%     else
%         avg_dvp = smean(trans.hists.dvporg(sitr:iter));
%         muNEWT.mu_fact = 1.2;
%         muNEWT.mu_delta = 5;
%     end
%     if iter == 300
%         trans.sigma = 1;
%     end
    if avg_dvp > muNEWT.mu_delta
        trans.sigma = trans.sigma * (muNEWT.mu_fact);
    else
        trans.sigma = trans.sigma / (muNEWT.mu_fact);
    end

    trans.sigma = min(muNEWT.mu_max, max(muNEWT.mu_min, trans.sigma));
end

end


% function [dx,dy,trans]= gendirectionxyori(Ftz,blk,At,trans,cgopts)
% nblock = trans.nblock;
% Amap  = trans.Amap;
% ATmap = trans.ATmap;
% sig = trans.sigma; %注意这里与SSNSDP的不同之处
% tau1 = trans.NEWT.tau1;
% tau2 = trans.NEWT.tau2;
% CG_maxit = cgopts.CG_maxit;
% pcgTol = cgopts.CG_tol;
% Lchol = trans.Lchol;
% Fx = Ftz.FX;
% Fy = Ftz.FY;
%
% iHW = Ftz.par;
% for k = 1:nblock
%     iHW.Dsch12{k} = (sig*iHW.Dsch12{k})./(1+tau2*sig-iHW.Dsch12{k});
% end
% ecoe1 = 1/tau2;
% N24Fx = DPhi_sdpnal2(blk,iHW,Fx,ecoe1);
% N24Fx = Amap(N24Fx);
% rhs1 = -Fy +N24Fx;
%
% iHWy = Ftz.par;
% iHWy.sig = 1;
% for k = 1:nblock
%     iHWy.Dsch12{k} = (iHWy.Dsch12{k}.*(1+sig*tau2)*sig)./(1+sig*tau2-iHWy.Dsch12{k});
% end
% ecoe2 = (1+sig*tau2)/tau2;
% iHWy.epsilon = tau1;
% L = struct();
% % L.invdiagM = trans.invdiagM;
% [dy,relres,flag] = psqmr_sdpnal2(@matvec_sdpnaly,blk,At,rhs1,ecoe2,iHWy,L,pcgTol,CG_maxit,Lchol);
% trans.cgres = relres;
% trans.flag = flag;
% trans.cgiter = length(relres);
%
% rhs2 = ATmap(dy);
% ecoe3 = 1;
% rhs2 = DPhi_sdpnal2(blk,Ftz.par,rhs2,ecoe3);
% for k = 1:nblock
%     rhs2{k} = rhs2{k}-Fx{k};
%     rhs2{k} = (rhs2{k})./(1+sig*tau2);
% end
% iHWx = Ftz.par;
% for k = 1:nblock
%     iHWx.Dsch12{k} = (iHWx.Dsch12{k}*sig)./(1+sig*tau2-iHWx.Dsch12{k});
% end
% %     dx= DPhi_sdpnalN4T(blk,iHWx,rhs2,sig,tau);
% dx= DPhi_sdpnal2(blk,iHWx,rhs2,ecoe1);
% for k = 1:nblock
%     dx{k} = dx{k} + sig * rhs2{k};
% end
% end
%
%
% function [FZ,S,X,trans] = compute_projection2(blk,X,y,b,C,trans)
% % if trans.Fzres < 1e-3
% [X.var,~] = project2_sdpnal(blk,X.var);
% % end
% tmp = ops_sdpnal(C,'-',ops_sdpnal(X.var,'/',trans.sigma));
% W = ops_sdpnal(trans.ATmap(y.var),'-',tmp);
% [S0,FZ.par] = project2_sdpnal(blk,W);
% FZ.par.precond = trans.precond; %加入预条件，为pre做准备
%
% S.var = ops_sdpnal(S0,'-',W);
% S.Avar = trans.Amap(S.var);
% FZ.FY = -b + trans.sigma*trans.Amap(S0);
% FZ.FX = ops_sdpnal(ops_sdpnal(X.var,'/',trans.sigma),'-',S0);
% FZ.res = ops_sdpnal(FZ.FY,'norm')+ops_sdpnal(FZ.FX,'norm');
% end