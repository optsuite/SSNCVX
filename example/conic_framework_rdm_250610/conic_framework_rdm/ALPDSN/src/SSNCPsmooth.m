function [X,out,y,S] = SSNCP(blk,At,C,b,opts)


% SSNSDP: a semi-smooth Newton method for solving SDP.
%
%   (P) max <C,X> s.t. A(X)=b, X \in K,
%   (D) min b^Ty s.t.  C = At(y) + S,
%   where K is a positive semidefinite cone.
%
%
% Acknowledgements: Some source codes are based on the implementation of
% SDPNAL, SDPNALplus and  ADMM+.  We thank  Defeng Sun, Kim-Chuan Toh,
% Xinyuan Zhao, Liuqin Yang for the tremendous efforts on developing these
% codes and kindly let us to use some parts of them.

%%**********************************************************************
%% SDPNAL+:
%% Copyright (c) 2014 by
%% Defeng Sun, Kim-Chuan Toh, Yancheng Yuan, Xinyuan Zhao
%% Corresponding author: Kim-Chuan Toh
%% Last modified 2016/06/01
%%**********************************************************************
%

addpath(genpath("./utils"));
%% Set the parameter
if nargin < 5; opts = []; end
opts = default_opts(opts);
record  = opts.record;
scale_data = opts.scale_data;
% trans.pmu = opts.pmu;% sigma
maxits = opts.maxits;
tol = opts.tol;
tolscale = opts.tolscale;

smooth_par_init = 1e-2;
smooth_par_decay = 0.1;
smooth_par_min = 1e-8;
zeta0 = 0.1;

muopts = opts.muopts;
cgopts = opts.cgopts;
NEWTopts = opts.NEWTopts;
trans.cgmin = opts.cgmin;
trans.cgmax = opts.cgmax;
trans.sigxl = opts.sigxl;
trans.sigxu = opts.sigxu;
trans.sigyl = opts.sigyl;
trans.sigyu = opts.sigyu;
trans.gfactor = opts.gfactor;
trans.gfactor2 = opts.gfactor2;  

%% record original variable and preprocessing
trans.borg = b;
trans.Corg = C;
trans.Atorg = At;
trans.blkorg = blk;
use_proximal  = 0; %% better choice: use_proximal  = 1;
% modify_sig = 0;
trans.precond = 0;
m = length(b);



trans.normborg = max(1,norm(trans.borg));
trans.normCorg = max(1,ops_sdpnal(trans.Corg,'norm'));
[blk,At,C,b,trans] = preprocess_SDP(blk,At,C,b,opts,trans);
trans.Amap  = @(X)AXmap(X,blk,At,trans.Lchol);
trans.ATmap = @(y)Atymap(y,blk,At,trans.Lchol);

trans.diagAAt = zeros(m,1);
for p = 1:size(blk,1); trans.diagAAt = trans.diagAAt + sum(At{p}.*At{p})'; end
trans.diagAAt = max(1e-4,full(trans.diagAAt));
trans.pmu = 1;

%% choose smoothing function

% huber function
smoothfun = @(lambda, mu) (lambda >= mu) .* (lambda - mu/2) + ...
    (lambda > 0 & lambda < mu) .* (lambda.^2 / (2*mu)) ;
d1smoothfun = @(lambda, mu) (lambda >= mu) .* 1 + ...
    (lambda < mu & lambda > 0) .* (lambda ./ mu) ;  % derivative with respect to lambda
d2smoothfun = @(lambda, mu) (lambda >= mu) .* (-0.5) + ...
    (lambda < mu & lambda > 0) .* (- 0.5 * lambda .^2 / mu^2) ;  % derivative with respect to mu



%% initial variable and operator
smooth_par = smooth_par_init;
trans.smooth_par = smooth_par;


X.var  = ops_sdpnal(C,'zeros');
X.Avar = trans.Amap(X.var);
Z.var  = ops_sdpnal(C,'zeros');
Z.Avar = trans.Amap(Z.var);
W = Z.var;
y.var  = zeros(trans.mdim,1);
y.Avar = X.var;
S.var = ops_sdpnal(C,'zeros');
S.Avar = trans.Amap(S.var);
% CS = cell(trans.nblock,1);

trans.bmAX = b - X.Avar;
trans.AC = trans.Amap(C);
trans.normb = max(1,norm(b));
trans.normC = max(1,ops_sdpnal(C,'norm'));
% [~, X, trans] = comp_res_fp(Z,blk,b,C,trans);


% trans = initial_ADMM_trans(trans,ADMMopts);
trans = initial_NEWT_trans(trans,NEWTopts);
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


% if (use_proximal)
%     sighatmin = 1e5; sighat = sighatmin; sighatfac = 10;
%     H = ones(m,1);
% else
%     trans.sighat = 0;
% end
%
% yhat = y;
% RpGradratio = 1;
%% One step ADMM (optional)
%------------------------------------------------------------------
% step 1: compute y
% y.var = (b-X.Avar)/trans.pmu + trans.AC - S.Avar;
% y.Avar = trans.ATmap(y.var);
% %------------------------------------------------------------------
% % step 2: compute X and S
% for k = 1:trans.nblock
%     W{k} = C{k}-(y.Avar{k}+X.var{k}/trans.pmu);
%     if strcmp(blk{k,1},'s')
%         W{k} = (W{k} + W{k}')/2;
%     end
% end
% % if (trans.hists.maxinf < 1e-1) && (trans.ADMM.recompeig < 5)
% %     trans.ADMM.can_use_eigs = 1; else; trans.ADMM.can_use_eigs = 0;
% % end
% [S.var,trans.ADMM.rankS,eighist] = blkprojSDP_sdpnal(blk,W);  %投影到S+
% Sm.var = ops_sdpnal(S.var,'-',W);
% Sm.Avar = trans.Amap(Sm.var);
% % trans.ADMM.recompeig = trans.ADMM.recompeig + max(eighist.recomp);
% 
% step = 1.618*trans.pmu;
% S.Avar = trans.Amap(S.var);
% for k = 1:trans.nblock
%     CS{k} = C{k} - S.var{k};
% end
% trans.AtymCSnrm = 0;
% trans.AtymCS = cell(trans.nblock,1);
% for k = 1:trans.nblock
%     trans.AtymCS{k} = y.Avar{k} - CS{k};
%     X.var{k} = X.var{k} + step*trans.AtymCS{k};%
%     trans.AtymCSnrm = trans.AtymCSnrm + norm(trans.AtymCS{k},'fro')^2;
% end
% trans.AtymCSnrm = sqrt(trans.AtymCSnrm);
% [X.Avar,X.Avarorg] = trans.Amap(X.var);
%% Compute projection
[FZ,S,trans] = compute_smooth_projection(blk,X,y,b,C,trans, smoothfun, d1smoothfun, d2smoothfun, smooth_par);
trans = set_param(trans);% 做牛顿步设置参数

count = 0; % cg restart
%% main iteration
for iter = 1:maxits
    %             FZ.FY=FZnew.FY;
    %             FZ.FX=FZnew.FX;
    FZ.FXres = ops_sdpnal(FZ.FX,'norm');
    FZ.FYres = ops_sdpnal(FZ.FY,'norm');
    FZ.res = FZ.FXres + FZ.FYres;

    zeta = zeta0 * min(1, (FZ.res ^ 2 + smooth_par^2)^(1 + 0.5));
    trans.zeta = zeta;

    if  use_proximal
        yhat = y.var;
        sighatmin = max(1,RpGradratio)*sighatmin; %% add: 27-Mar-2008
        sighatfac = max(1,RpGradratio^2)*sighatfac;
        sighat = max(sighatmin,sighatfac*trans.pmu);    %% add: 07-Apr-2008
        H = 1./max(1,sqrt(abs(yhat)));  %% better than H = ones(m,1);
        H2 = H.*H;
        %         FZ.par.yhat = yhat; FZ.par.sighat = sighat; FZ.par.H = H; FZ.par.H2 = H2;
    end

    if (use_proximal)
        trans.invdiagM = 1./full(trans.pmu*trans.diagAAt + H2/sighat);
    else
        trans.invdiagM = 1./full(trans.pmu*trans.diagAAt);
    end
    trans.NEWT.iter = trans.NEWT.iter + 1;
    trans.NEWT.totaliter = trans.NEWT.totaliter + 1;
    % parameter for the revised Newton system
    if false
        trans.NEWT.tau1 = 1e-4;
        trans.NEWT.tau2 = 1e-4;
    else
        trans.NEWT.tau1 = trans.NEWT.lambda*(FZ.FYres^trans.NEWT.sigPowy);
        trans.NEWT.tau2 = trans.NEWT.lambda*(FZ.FXres^trans.NEWT.sigPowx);
    end

    %         trans.NEWT.sig = 1 ;
    if cgopts.CG_adapt %更改CG的tol
        cgopts.CG_tol = max(min(0.1*FZ.FYres, 0.1),cgopts.cgtolmin)*0.001;% FZ.res
    end

    %         [d,trans] = gendirection(FZ,blk,At,trans,cgopts,NEWTopts,scale_data);
    [dx,dy,trans] = gendirectionxy(FZ,blk,At,trans,cgopts,NEWTopts,scale_data);
    for k = 1:trans.nblock
        if strcmp(blk{k,1},'s')
            dx{k} = (dx{k} + dx{k}')/2;
        end
    end
    %----------------------------------------------------------------------
    ynew = struct();
    ynew.var = y.var + dy;
    Xnew = struct();
    ynew.Avar = trans.ATmap(ynew.var);
    Xnew.var = ops_sdpnal(X.var,'+',dx) ;
    [AX,AXorg] = trans.Amap(Xnew.var);
    Xnew.Avar = AX;
    Xnew.Avarorg = AXorg;

    [FZnew,S,trans] = compute_smooth_projection(blk,Xnew,ynew,b,C,trans, smoothfun, d1smoothfun, d2smoothfun, smooth_par);

    trans.NEWT.nrmd = ops_sdpnal(dx,'norm')+ops_sdpnal(dy,'norm');
    trans.NEWT.FZd = -sum(FZnew.FY.*dy);
    for k = 1:trans.nblock
        trans.NEWT.FZd = trans.NEWT.FZd - sum(sum(FZnew.FX{k}.*dx{k}));
    end
    trans.NEWT.FZdorg = -sum(FZ.FY.*dy);
    for k = 1:trans.nblock
        trans.NEWT.FZdorg = trans.NEWT.FZdorg - sum(sum(FZ.FX{k}.*dx{k}));
    end
    % thetaFZdorg = trans.NEWT.FZdorg/FZ.res/trans.NEWT.nrmd;
    trans.NEWT.FZFZnew = sum(FZ.FY.*FZnew.FY);
    for k = 1:trans.nblock
        trans.NEWT.FZFZnew = trans.NEWT.FZFZnew + sum(sum(FZ.FX{k}.*FZnew.FX{k}));
    end
    % thetaFZFZnew = trans.NEWT.FZFZnew/FZnew.res/FZ.res;

    if FZ.res < 1e-3
        trans.NEWT.rhs = FZnew.res*trans.NEWT.nrmd^2;
    else
        trans.NEWT.rhs = trans.NEWT.nrmd^2;
    end
    trans.NEWT.ratio = trans.NEWT.FZd/trans.NEWT.rhs;

    % condition on performing a newton step 解的不好的时候
    if (trans.NEWT.iter == 1) && (FZnew.res >= FZ.res)
        trans.NEWT.lambda = trans.NEWT.lambda*trans.gfactor2;
        %             trans.NEWT.swt = trans.NEWT.swt + 1;
        continue;
    end

    if FZnew.res < 5e1*FZ.res || count <= 5 %5e1*FZ.res
        count = 0;
        y = ynew;
        FZ = FZnew;
        X = Xnew;
        trans = sigPow_CGmaxit_update(FZ.res,trans);
        cgopts.CG_maxit = trans.NEWT.CG_maxit;
    else
        %             trans.NEWT.tau1 = trans.NEWT.lambda*(FZ.FYres^trans.NEWT.sigPowy);
        %             trans.NEWT.tau2 = trans.NEWT.lambda*(FZ.FXres^trans.NEWT.sigPowx);
        trans.NEWT.lambda = trans.gfactor* trans.NEWT.lambda;
        count = count + 1;
        continue;
    end
 
    %------------------------------------------------------------------
    % check optimality
    trans = record_optimal_NEWT(C,b,trans,X,y,S,iter);

    if (max(trans.hists.pinforg(end),trans.hists.dinforg(end)*tolscale) < tol*retol)
        if opts.tolscale < 100
            [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,FZ,X,y,S,W,trans);
            if max([rec.pinf,rec.dinf*tolscale,rec.K1,rec.K1dual,rec.C1]) < tol*retol
                cstop = 1;
                X = Xnew; y = ynew; S = Snew;
                trans.rec = rec;
            else
                cstop = 0;
            end
        else
            cstop = 1;
            [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,FZ,X,y,S,W,trans);
            X = Xnew; y = ynew; S = Snew;
            trans.rec = rec;
        end
    end

    %% record log
    
    log_info = {'%5s', 'iter', num2str(iter);
                '%7s', 'ALpen', num2str(trans.pmu, '%2.1e');
                '%7s', 'smtpar', num2str(smooth_par, '%2.1e');
                '%8s', 'pobj', num2str(trans.hists.pobj(end), '%2.1e');
                '%8s', 'dobj', num2str(trans.hists.dobj(end), '%2.1e');
                '%7s', 'gap', num2str(trans.hists.gap(end), '%2.1e');
                '%7s', 'gaporg', num2str(trans.hists.gaporg(end), '%2.1e');
                '%7s', 'pinforg', num2str(trans.hists.pinforg(end), '%2.1e');
                '%7s', 'dinforg', num2str(trans.hists.dinforg(end), '%2.1e');
                '%7s', 'res', num2str(FZ.res, '%2.1e');
                '%7s', 'tau1', num2str(trans.NEWT.tau1, '%2.1e');
                '%7s', 'tau2', num2str(trans.NEWT.tau2, '%2.1e');
                '%7s', 'CGtol', num2str(trans.cgres(end), '%2.1e');
                '%5s', 'CGit', num2str(trans.cgiter, '%5d');
                '%7s', 'ratio', num2str(trans.NEWT.ratio, '%2.1e');
                '%7s', 'kappa', num2str(trans.NEWT.lambda, '%2.1e');
                '%5s', 'zeta', num2str(zeta, '%2.1e');
                '%7s', 'newtres', num2str(trans.newton_res, '%2.1e');};
                % '%7s', 'thetaFZdorg', num2str(thetaFZdorg, '%2.1e');
                % '%7s', 'thetaFZFZnew', num2str(thetaFZFZnew, '%2.1e');
                % '%8s', 'tolscale', num2str(tolscale, '%2.1e');
                % '%7s', 'retol', num2str(retol, '%2.1e')};

    %% define header for printing
    if record >= 1 && ~ exist('str_head2','var')
        str_head2 = "";
        for i = 1:size(log_info,1)
            str_head2 = str_head2 + sprintf(log_info{i,1}, log_info{i,2}) + " ";
        end
        fprintf('\n%s', str_head2);
    end

    %% print iteration information
    if record>=1 && (cstop || ...
            iter == 1 || iter==maxits || mod(iter,NEWTopts.print_itr)==0)
        if mod(iter,20*NEWTopts.print_itr) == 0 && iter ~= maxits && ~cstop
            fprintf('\n%s', str_head2);
        end

        str_info = "";
        for i = 1:size(log_info,1)
            str_info = str_info + sprintf(log_info{i,1}, log_info{i,3}) + " ";
        end
        fprintf('\n%s', str_info);

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
    if (iter == maxits)
        out.status = sprintf('reach the maximum iteration');
        [Xnew,ynew,Snew,rec] = recover_var_chk2nd(Z,FZ,X,y,S,W,trans);
        X = Xnew; y = ynew; S = Snew;
        trans.rec = rec;
        break;
    end


    %% update smooth_par
    smooth_par = smooth_par * zeta;
    if smooth_par > 0 && smooth_par < smooth_par_min
        smooth_par = 0;
        fprintf('\n  -- iter %d phaseII: smooth_par = 0', iter);
    end
    trans.smooth_par = smooth_par;



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
            pmup = trans.pmu;
            trans = mu_update_NEWT(trans,iter,muopts,opts);

            if trans.pmu ~= pmup
                if record >=1; fprintf('\n  -- AugLag penalty updated: %f',trans.pmu); end
                %                 [FZ, X, Z] = comp_res_fp_updmu(X, Z, blk, b, C, trans.pmu,pmup,trans);
            end
            %             trans = mu_update_ADMM(trans,iter,muopts);
    end


end

out = generate_outs(out,X,y,S,trans,opts,iter);

end


%% End main routine





function opts = default_opts(opts)
%% ------------------------------------------------------------------
if ~isfield(opts, 'record');      opts.record = 1;       end
if ~isfield(opts, 'scale_data');  opts.scale_data = 1;   end
if ~isfield(opts, 'pmu');         opts.pmu = 1;          end
if ~isfield(opts, 'sstruct');     opts.sstruct = 1;      end
if ~isfield(opts, 'maxits');      opts.maxits = 1000;   end
if ~isfield(opts,'tol');          opts.tol = 1e-6;       end
if ~isfield(opts,'tolscale');     opts.tolscale = 1;     end



%% -------------------------------------------------------------------------
% options for fixed-point algorithms
if ~isfield(opts, 'ADMMopts'); opts.ADMMopts = struct; end
ADMMopts = opts.ADMMopts;
if ~isfield(ADMMopts, 'print_itr'); ADMMopts.print_itr = 20;     end
if ~isfield(ADMMopts, 'rho');       ADMMopts.rho = 1.618;        end
if ~isfield(ADMMopts, 'imaxit');    ADMMopts.imaxit = 20000;     end
opts.ADMMopts = ADMMopts;

%% ------------------------------------------------------------------------
% options for the semismooth Newton algorithm
if ~isfield(opts, 'NEWTopts'); opts.NEWTopts = struct; end
NEWTopts = opts.NEWTopts;
if ~isfield(NEWTopts,'sigPow');   NEWTopts.sigPow = 1;       end
if ~isfield(NEWTopts,'resFac');   NEWTopts.resFac = 0.98;      end
if ~isfield(NEWTopts,'eta1');     NEWTopts.eta1 = 1e-4;        end
if ~isfield(NEWTopts,'eta2');     NEWTopts.eta2 = 0.9;         end
if ~isfield(NEWTopts,'gamma1');   NEWTopts.gamma1 = 0.5;       end %可调
if ~isfield(NEWTopts,'gamma2');   NEWTopts.gamma2 = 0.95;      end %可调
if ~isfield(NEWTopts,'gamma3');   NEWTopts.gamma3 = 10;         end %可调
if ~isfield(NEWTopts,'lambda');   NEWTopts.lambda = 1;         end %adjust mu
%     if ~isfield(NEWTopts,'sigma'),    NEWTopts.sigma = 0.01;       end
%     if ~isfield(NEWTopts,'NTstep');   NEWTopts.NTstep  = 5;        end
if ~isfield(NEWTopts, 'print_itr'); NEWTopts.print_itr = 1;    end
opts.NEWTopts = NEWTopts;


%% ------------------------------------------------------------------------
% parameter for CG
if ~isfield(opts, 'cgopts'); opts.cgopts = struct; end
cgopts = opts.cgopts;
if ~isfield(cgopts,'CG_maxit'); cgopts.CG_maxit = 500;   end %可调
if ~isfield(cgopts,'CG_tol');   cgopts.CG_tol = 1e-2;    end
if ~isfield(cgopts,'cgtolmin'); cgopts.cgtolmin = 1e-10;   end
if ~isfield(cgopts,'CG_adapt'); cgopts.CG_adapt = 1;     end
opts.cgopts = cgopts;


%% ------------------------------------------------------------------
% parameters for adjusting pmu
if ~isfield(opts, 'muopts'); opts.muopts = struct; end
muopts = opts.muopts;
if ~isfield(muopts, 'adp_mu');        muopts.adp_mu = 1;         end %1 or 0

if ~isfield(muopts,'NEWT');            muopts.NEWT = struct; end
% if ~isfield(muopts.NEWT, 'adpmu_cri');     muopts.NEWT.adpmu_cri = 1;      end %1 or 2
if ~isfield(muopts.NEWT, 'mu_min');        muopts.NEWT.mu_min = 1e-6;      end
if ~isfield(muopts.NEWT, 'mu_max');        muopts.NEWT.mu_max = 1e6;       end
if ~isfield(muopts.NEWT, 'mu_update_itr'); muopts.NEWT.mu_update_itr = opts.muopts.mu_update_itr; end %10 可调
if ~isfield(muopts.NEWT, 'mu_delta');      muopts.NEWT.mu_delta = 5e-1;    end % 可调
if ~isfield(muopts.NEWT, 'mu_fact');       muopts.NEWT.mu_fact = opts.muopts.mu_fact;      end %5/3 可调

% if ~isfield(muopts,'ADMM');          muopts.ADMM = struct; end
% if ~isfield(muopts.ADMM, 'ratio');   muopts.ADMM.ratio = 1.;      end %1 or 0
% if ~isfield(muopts.ADMM, 'mu_max');   muopts.ADMM.mu_max = 1e6;     end %1 or 2
% if ~isfield(muopts.ADMM, 'mu_min');   muopts.ADMM.mu_min = 1e-4;    end
% if ~isfield(muopts.ADMM, 'mu_delta');muopts.ADMM.mu_delta = 1.2;  end
% if ~isfield(muopts.ADMM, 'mu_fact'); muopts.ADMM.mu_fact = 1.8;  end

opts.muopts = muopts;
end

function [blk,At,C,b,trans] = preprocess_SDP(blk,At,C,b,opts,trans)

%% scale the data and cholesky decomposition
trans.mdim = length(b);
trans.nblock = size(blk,1);
DA = speye(trans.mdim);
bscale = 1;
Cscale = 1;

if opts.scale_data
    normA = zeros(trans.mdim,1);
    for k = 1:(trans.nblock)
        normA = normA + sum(At{k}.*At{k})';
    end
    normA = max(1,sqrt(normA));
    DA = spdiags(1./normA,0,trans.mdim,trans.mdim);
end
for k = 1:trans.nblock
    At{k} = At{k}*DA; %这里是乘 后面算org还要除去
end
b = DA*b;
scale.DA = DA;


[Lchol,AAt] = cholAAt_sdpnal(blk,At,trans.mdim);
Lchol.isidentity = false;
b = fwsolve(Lchol,b);

if (opts.scale_data==1)
    bscale = max(1,norm(b));
    Cscale = max(1,ops_sdpnal(C,'norm'));
end

b = b/bscale; %单位化
C = ops_sdpnal(C,'/',Cscale);
objscale = bscale*Cscale;
scale.bscale = bscale;
scale.Cscale = Cscale;
scale.objscale = objscale;
scale.scale_data = opts.scale_data;
trans.scale = scale;
trans.Lchol = Lchol;
end

function trans = initial_ADMM_trans(trans,ADMMopts)
trans.ADMM.rankS = zeros(trans.nblock,1);
trans.ADMM.recompeig = 0;
% trans.ADMM.can_use_eigs = 0;
trans.ADMM.swt = 0;
trans.ADMM.maxiter = ADMMopts.imaxit;
trans.ADMM.subiter = 0;
trans.ADMM.iter = 0;
trans.muADMM.prim_win = 0;
trans.muADMM.dual_win = 0;
end

function trans = initial_NEWT_trans(trans,NEWTopts)
trans.NEWT.lambda = NEWTopts.lambda;
trans.NEWT.sigPowx = NEWTopts.sigPow;
trans.NEWT.sigPowy = NEWTopts.sigPow;
trans.NEWT.iter = 0;
trans.NEWT.swt = 0;
trans.NEWT.subiter = 0;
trans.NEWT.maxiter = 20;
trans.NEWT.lastres = inf;
end

function trans = initial_hists_trans(trans)
trans.hists.maxinf = inf;
end


function [AX, AXorg] = AXmap(X,blk,At,Lchol)
AXorg = AXfun_sdpnal(blk,At,X); %原始的AX 在process_SDP函数那里A是经过除去对角的。
AX  = fwsolve(Lchol,AXorg);
end

function Aty = Atymap(y,blk,At,Lchol)
Aty = Atyfun_sdpnal(blk,At,bwsolve(Lchol,y));
end

function q = fwsolve(L,r)
if (L.isidentity)
    q = r;
else
    if strcmp(L.matfct_options,'chol')
        q = mextriang_sdpnal(L.R,r(L.perm),2);
    elseif strcmp(L.matfct_options,'spcholmatlab')
        q = mexfwsolve_sdpnal(L.R,r(L.perm,1));
    end
end
end

function q = bwsolve(L,r)

if (L.isidentity)
    q = r;
else
    if strcmp(L.matfct_options,'chol')
        q(L.perm,1) = mextriang_sdpnal(L.R,r,1);
    elseif strcmp(L.matfct_options,'spcholmatlab')
        q(L.perm,1) = mexbwsolve_sdpnal(L.Rt,r);
    end
end
end


function [FZ, X, trans] = comp_res_fp(Z,blk,b,C,trans)

pmu = trans.pmu;
[X.var,FZ.par] = project2_sdpnal(blk,Z.var);
[AX,AXorg] = trans.Amap(X.var);
X.Avar = AX;
X.Avarorg = AXorg;
mud = 2*AX - Z.Avar - pmu*trans.AC-b;
FZ.var = trans.ATmap(mud);
for k = 1:trans.nblock
    FZ.var{k} = FZ.var{k} + pmu*C{k} + Z.var{k} - X.var{k};
end
FZ.res = ops_sdpnal(FZ.var,'norm');
FZ.Avar = AX - b;
FZ.par.precond = trans.precond;
end

function [FZ,S,trans] = compute_projection(blk,X,y,b,C,trans)
tmp = ops_sdpnal(C,'-',ops_sdpnal(X.var,'/',trans.pmu));
W = ops_sdpnal(trans.ATmap(y.var),'-',tmp); % W = A' * y - C + X/pmu
[S0,FZ.par] = project2_sdpnal(blk,W);
FZ.par.precond = trans.precond; %加入预条件，为pre做准备

S.var = ops_sdpnal(S0,'-',W);
S.Avar = trans.Amap(S.var);
FZ.FY = -b + trans.pmu*trans.Amap(S0);
FZ.FX = ops_sdpnal(ops_sdpnal(X.var,'/',trans.pmu),'-',S0);
FZ.res = ops_sdpnal(FZ.FY,'norm')+ops_sdpnal(FZ.FX,'norm');
end


function [FZ,S,trans] = compute_smooth_projection(blk,X,y,b,C,trans, smoothfun, d1smoothfun, d2smoothfun, smooth_par)
    if smooth_par == 0
        [FZ,S,trans] = compute_projection(blk,X,y,b,C,trans);
        return;
    end
    tmp = ops_sdpnal(C,'-',ops_sdpnal(X.var,'/',trans.pmu));
    W = ops_sdpnal(trans.ATmap(y.var),'-',tmp); % W = A' * y - C + X/pmu
    [S0, FZ.par, d2S0] = smooth_projection(blk, W, smoothfun, d1smoothfun, d2smoothfun, smooth_par);
    FZ.par.precond = trans.precond; %加入预条件，为pre做准备
    
    S.var = ops_sdpnal(S0,'-',W);
    S.Avar = trans.Amap(S.var);
    FZ.FY = -b + trans.pmu*trans.Amap(S0);
    FZ.d2FY = trans.pmu * trans.Amap(d2S0);
    FZ.FX = ops_sdpnal(ops_sdpnal(X.var,'/',trans.pmu),'-',S0);
    FZ.d2FX = ops_sdpnal(0, '-', d2S0);
    FZ.res = ops_sdpnal(FZ.FY,'norm')+ops_sdpnal(FZ.FX,'norm');
end
    



function [dx,dy,trans]= gendirectionxy(Ftz,blk,At,trans,cgopts,NEWTopts,scale_data)

    %% [N1 N2;   * [dy;  = [rhsy;
    %   N3 N4]      dx]    rhsx]
nblock = trans.nblock;
Amap  = trans.Amap;
ATmap = trans.ATmap;
sig = trans.pmu; %注意这里与SSNSDP的不同之处
tau1 = trans.NEWT.tau1;
tau2 = trans.NEWT.tau2;
CG_maxit = cgopts.CG_maxit;
pcgTol = cgopts.CG_tol;
Lchol = trans.Lchol;
Fx = Ftz.FX;
Fy = Ftz.FY;
smooth_par = trans.smooth_par;
if smooth_par > 0
    d2Fx = Ftz.d2FX;
    d2Fy = Ftz.d2FY;
    rhsx = ops_sdpnal(ops_sdpnal(trans.zeta * smooth_par, '*', d2Fx), '-', Fx);
    rhsy = ops_sdpnal(ops_sdpnal(trans.zeta * smooth_par, '*', d2Fy), '-', Fy);
else
    rhsx = ops_sdpnal(0, '-', Fx);
    rhsy = ops_sdpnal(0, '-', Fy);
end
%Ftz.par: \Sigma


iHW = Ftz.par;
for k = 1:nblock
    if trans.smooth_par > 0
        iHW.Dsch11{k} = (sig*iHW.Dsch11{k})./(1+tau1*sig-iHW.Dsch11{k}); %iHW: \tilde{\Sigma}
    end
    iHW.Dsch2{k} = (sig*iHW.Dsch2{k})./(1+tau2*sig-iHW.Dsch2{k}); %iHW: \tilde{\Sigma}
end
ecoe1 = 1/tau2;
if smooth_par > 0
    N24rhsx = DPhi_smooth(blk,iHW,rhsx);  
else
    N24rhsx = DPhi_sdpnal(blk,iHW,rhsx,ecoe1);  
end
N24rhsx = Amap(N24rhsx);             %N_2(N_4)^{-1} Fx 
rhs1 = rhsy - N24rhsx;

iHWy = Ftz.par;
iHWy.sig = 1;
for k = 1:nblock
    if trans.smooth_par > 0
        iHWy.Dsch11{k} = (iHWy.Dsch11{k}.*(1+sig*tau1)*sig)./(1+sig*tau1-iHWy.Dsch11{k}); %iHWy: \hat{\Sigma}
    end
    iHWy.Dsch2{k} = (iHWy.Dsch2{k}.*(1+sig*tau2)*sig)./(1+sig*tau2-iHWy.Dsch2{k}); %iHWy: \hat{\Sigma}
end
ecoe2 = (1+sig*tau2)/tau2;
iHWy.epsilon = tau1;
L = struct();
L.invdiagM = trans.invdiagM;
if trans.smooth_par > 0
    matvec = @matvec_smooth;
else
    matvec = @matvec_sdpnaly;
end

[dy,relres,flag] = psqmr_sdpnal(matvec,blk,At,rhs1,ecoe2,iHWy,L,pcgTol,CG_maxit,Lchol);

% %%
% %% use direct method to solve linear system
% % lhs = Rt^{-1} * S * A * diag(iHWy.Dsch21) * At * Sinv * R^{-1} + epsilon * I
% % Sinv * y is equivalent to y(perm) = y
% % S * y is equivalent to y = y(perm)
% if ~isfield(Lchol, 'RtR') 
%     Lchol.RtR = Lchol.Rt * Lchol.R; 
%     Lchol.RtR(Lchol.perm, Lchol.perm) = Lchol.RtR;
% end

% for p =1: nblock
%     if p == 1
%         ADAt = At{p}' * spdiags(iHWy.Dsch2{p}, 0, length(iHWy.Dsch2{p}), length(iHWy.Dsch2{p})) * At{p};
%         %AorgDAt = trans.Atorg{p}' * spdiag(iHWy.Dsch2{p}) * trans.Atorg{p};
%     else
%         ADAt = ADAt + At{p}' * spdiags(iHWy.Dsch2{p}, 0, length(iHWy.Dsch2{p}), length(iHWy.Dsch2{p})) * At{p};
%         %AorgDAt = AorgDAt + trans.Atorg{p}' * spdiag(iHWy.Dsch2{p}) * trans.Atorg{p};
%     end
% end

% ADAt = iHWy.sig * ADAt + iHWy.epsilon * Lchol.RtR;

% rhs1_temp = zeros(size(rhs1));
% rhs1_temp(Lchol.perm) = Lchol.Rt * rhs1;
% dy = ADAt \ rhs1_temp;
% dy = Lchol.R * dy(Lchol.perm);

% %% check residual
% dy_temp = zeros(size(dy));
% if strcmp(Lchol.matfct_options,'chol')
%     dy_temp(Lchol.perm) = Lchol.R \ dy;
% elseif strcmp(Lchol.matfct_options,'spcholmatlab')
%     dy_temp(Lchol.perm) = mexbwsolve_sdpnal(Lchol.Rt, dy);
% end
% for p =1: nblock
%     if p == 1
%         lhs_temp = At{p}' * (reshape(iHWy.Dsch2{p}, [], 1) .* (At{p} * dy_temp)) ;
%     else
%         lhs_temp = lhs_temp + At{p}' * (reshape(iHWy.Dsch2{1}, [], 1) .* (At{p} * dy_temp)) ;
%     end
% end
% if strcmp(Lchol.matfct_options,'chol')
%     lhs_temp = Lchol.Rt \ lhs_temp(Lchol.perm);
% elseif strcmp(Lchol.matfct_options,'spcholmatlab')
%     lhs_temp = mexfwsolve_sdpnal(Lchol.R, lhs_temp(Lchol.perm));
% end
% lhs = iHWy.epsilon * dy + iHWy.sig * lhs_temp;
% relres = norm(rhs1 - lhs) / (1 + norm(rhs1));
% % fprintf("resdue of computing dy: %e\n", relres);
% flag = 1;



 %%

trans.cgres = relres;
trans.flag = flag;
trans.cgiter = length(relres);

rhs2 = ATmap(dy);
ecoe3 = 1;
if smooth_par > 0
    rhs2 = DPhi_smooth(blk,Ftz.par,rhs2);
else
    rhs2 = DPhi_sdpnal(blk,Ftz.par,rhs2,ecoe3);
end
for k = 1:nblock
    rhs2{k} = rhsx{k} + rhs2{k};          %  rhs2 = rhsx - N3 dy = rhsx + D * At * rhs2
end 
iHWx = Ftz.par;
for k = 1:nblock
    if smooth_par > 0
        iHWx.Dsch11{k} = (iHWx.Dsch11{k}*sig)./(1+sig*tau2-iHWx.Dsch11{k});  % \tilde{\Sigma}_{\tau_2}
    end
    iHWx.Dsch2{k} = (iHWx.Dsch2{k}*sig)./(1+sig*tau2-iHWx.Dsch2{k});  % \tilde{\Sigma}_{\tau_2}
end
if smooth_par > 0
    dx = DPhi_smooth(blk,iHWx,rhs2);  % dx =  (N4)^{-1} * rhs2
else
    dx = DPhi_sdpnal(blk,iHWx,rhs2,ecoe1); 
end
for k = 1:nblock
    dx{k} = 1 / (1+sig*tau2) * dx{k} + sig/(1+sig*tau2) * rhs2{k};
end

%% check residual
if smooth_par > 0
    D_lmut = @(x) DPhi_smooth(blk, Ftz.par, x);
else
    D_lmut = @(x) DPhi_sdpnal(blk, Ftz.par, x, 1);
end
Atdy = ATmap(dy);
N1dy = sig * Amap(D_lmut(Atdy)) + tau1 * dy;
N2dx = Amap(D_lmut(dx));
N3dy = ops_sdpnal(0, '-', D_lmut(Atdy));
N4dx = ops_sdpnal(ops_sdpnal((1 / sig + tau2), '*', dx), '-', ops_sdpnal(1 / sig, '*', D_lmut(dx)));
resy = ops_sdpnal(ops_sdpnal(N1dy, '+', N2dx), '-', rhsy);
resx = ops_sdpnal(ops_sdpnal(N3dy ,'+', N4dx), '-', rhsx);
trans.newton_res = sqrt(ops_sdpnal(resy, 'norm') ^ 2 + ops_sdpnal(resx, 'norm') ^ 2) / (1 + sqrt(ops_sdpnal(rhsy, 'norm') ^ 2 + ops_sdpnal(rhsx, 'norm') ^ 2));
end

function trans = record_optimal_NEWT(C,b,trans,X,y,S,iter)
% check optimality
% bTy  = full(b'*(trans.AC + (Z.Avar-X.Avar)/trans.pmu));
bTy  = full(b'*y.var);
trCX = full(ops_sdpnal(ops_sdpnal(C,'.*',X.var),'sum'));
dobj = trans.scale.objscale*bTy;
pobj = trans.scale.objscale*trCX;
gap = abs(trCX-bTy)/max(1,abs(trCX));
gaporg = abs(pobj-dobj)/max(1,abs(pobj));

bmAX = b-X.Avar;
bmAXnrm = norm(X.Avarorg./diag(trans.scale.DA)*trans.scale.bscale-trans.borg);
% pinforg = bmAXnrm/trans.normborg;
pinforg = bmAXnrm/(1+norm(trans.borg));
pinf = norm(bmAX)/trans.normb;

% S.var = ops_sdpnal(S.var,'*',trans.scale.Cscale);

y = y.var;
% S = S.var;
% scale = trans.scale;
y = trans.scale.Cscale*(trans.scale.DA*bwsolve(trans.Lchol,y));
% trans.AtymCSnrm = FZ.res/trans.pmu;
trans.AtymCSnrm = ops_sdpnal(ops_sdpnal(ops_sdpnal(Atyfun_sdpnal(trans.blkorg,trans.Atorg,y),'+',ops_sdpnal(S.var,'*',trans.scale.Cscale)),'-',trans.Corg),'norm');
dinf = trans.AtymCSnrm/trans.normC;
dinforg = trans.AtymCSnrm/(1+ops_sdpnal(trans.Corg,'norm'));

trans.hists.pobj(iter) = pobj;
trans.hists.dobj(iter) = dobj;
trans.hists.gaporg(iter)  = gaporg;
trans.hists.gap(iter)  = gap;
trans.hists.pinf(iter) = pinf;
trans.hists.pinforg(iter) = pinforg;
trans.hists.dinf(iter) = dinf;
trans.hists.dinforg(iter) = dinforg;
trans.hists.pvd(iter)  = pinf/dinf;
trans.hists.dvp(iter)  = dinf/pinf;
trans.hists.pvdorg(iter)  = pinforg/dinf;
trans.hists.dvporg(iter)  = dinf/pinforg;
trans.hists.cgiter(iter) = trans.cgiter;

maxinf = max(pinf,dinf);

trans.hists.isNEWT(iter) = 1;
trans.hists.maxinf = maxinf;

trans.bmAX = bmAX;
end

function [X,y,S,rec] = recover_var_chk2nd(Z,FZ,X,y,S,W,trans)
%     if strcmp(trans.last_iter, 'NEWT')
%         Zold = Z;
%         Z.var = Z.var - F.var;
%         Z.Avar = Z.Avar - FZ.Avar;
%         y = trans.AC + (Z.Avar - X.Avar)/trans.pmu;
%         S = X.var - Zold.var;
%         S = 1 / trans.pmu * S;
%         X = X.var;
%     else
y = y.var;
%         X = trans.pmu * (S.var - W);
X = X.var;
S = S.var;
%     end

scale = trans.scale;
X = scale.bscale * X;
y = scale.Cscale * (scale.DA*bwsolve(trans.Lchol,y));
S = scale.Cscale * S;


%% compute optimal index
AX = AXfun_sdpnal(trans.blkorg,trans.Atorg,X);
bmAXnrm = norm(AX-trans.borg);
rec.pinf = bmAXnrm/(1+norm(trans.borg));
AtymCSnrm = ops_sdpnal(ops_sdpnal(ops_sdpnal(Atyfun_sdpnal(trans.blkorg,trans.Atorg,y),'+',S),'-',trans.Corg),'norm');
rec.dinf = AtymCSnrm/(1+ops_sdpnal(trans.Corg,'norm'));
pobj = full(ops_sdpnal(ops_sdpnal(trans.Corg,'.*',X),'sum'));
dobj = full(trans.borg'*y);
rec.gap = abs(pobj - dobj)/(1+abs(pobj)+abs(dobj));
rec.pobj = pobj;
rec.dobj = dobj;
rec.totaltime = etime(clock,trans.tstart);
trXS = dot(X, S);

normX = ops_sdpnal(X,'norm');
normS = ops_sdpnal(S,'norm');
par = [];
[~,Xpar] = project2_sdpnal(trans.blkorg,X,par);
[~,Spar] = project2_sdpnal(trans.blkorg,S,par);
eigXnorm = 0;
eigSnorm = 0;
for k = 1:size(trans.blkorg,1)
    if strcmp(trans.blkorg{k,1},'s')
        rx = length(Xpar.posidx{k});
        rs = length(Spar.posidx{k});
        xnegidx = (rx+1):sum(trans.blkorg{k,2});
        snegidx = (rs+1):sum(trans.blkorg{k,2});
    elseif strcmp(trans.blkorg{k,1},'l')
        xnegidx = setdiff(1:trans.blkorg{k,2},Xpar.posidx{k});
        snegidx = setdiff(1:trans.blkorg{k,2},Spar.posidx{k});
    end
    eigXnorm = eigXnorm + norm(Xpar.dd{k}(xnegidx),'fro')^2;
    eigSnorm = eigSnorm + norm(Spar.dd{k}(snegidx),'fro')^2;
end
eigXnorm = sqrt(eigXnorm);
eigSnorm = sqrt(eigSnorm);
rec.K1 = eigXnorm/(1+normX);
rec.K1dual = eigSnorm/(1+normS);
rec.C1 = abs(trXS)/(1+normX+normS);
end


function out = generate_outs(out,X,y,S,trans,opts,iter)


rec = trans.rec;
if opts.record
    fprintf('\n-------------------------------------------------------------------------------------------------\n')
    fprintf('%14s %14s %10s %10s %10s %10s %10s %12s\n','pboj','dobj','gap', 'pinf','dinf','C1','K1','K1dual')
    fprintf('%14.8e %14.8e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n',...
        rec.pobj,rec.dobj,rec.gap,rec.pinf,rec.dinf,rec.C1,rec.K1,rec.K1dual)
    fprintf('\n-------------------------------------------------------------------------------------------------\n')
end

trans.hists.pobj(iter+2) = rec.pobj;
trans.hists.dobj(iter+2) = rec.dobj;
trans.hists.pinf(iter+2) = rec.pinf;
trans.hists.pinforg(iter+2) = rec.pinf;
trans.hists.dinf(iter+2) = rec.dinf;
trans.hists.dinforg(iter+2) = rec.dinf;
trans.hists.gap(iter+2) = rec.gap;
trans.hists.gaporg(iter+2) = rec.gap;
trans.hists.pvd(iter+2)  = rec.pinf/rec.dinf;
trans.hists.dvp(iter+2)  = rec.dinf/rec.pinf;
trans.hists.pvdorg(iter+2) = rec.pinf/rec.dinf;
trans.hists.dvporg(iter+2) = rec.dinf/rec.pinf;
%trans.hists.res(iter+2) = res;

out.pobj = rec.pobj;
out.dobj = rec.dobj;
out.gap  = rec.gap;
out.iter  = iter;
out.C1 = rec.C1;
out.K1 = rec.K1;
out.K1dual = rec.K1dual;
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

function trans = sigPow_CGmaxit_update(res,trans)
if res < 1e-4
    sigPowy      = trans.sigyl;
    sigPowx      = trans.sigxl;
else
    sigPowy      = trans.sigyu;
    sigPowx      = trans.sigxu;
end
CG_maxit = trans.cgmax;
% if min(trans.cgres(end),trans.cgres(end-1))<1e-6 && trans.cgiter < CG_maxit
%     CG_maxit = trans.cgmin;
% end

trans.NEWT.sigPowx = sigPowx;
trans.NEWT.sigPowy = sigPowy;
trans.NEWT.CG_maxit = CG_maxit;
end



function trans = mu_update_NEWT(trans,iter,muopts,opts)
smean = @geo_mean;
muNEWT = muopts.NEWT;
if mod(iter,muNEWT.mu_update_itr)==0
    sitr = iter-muNEWT.mu_update_itr+1;
    avg_pvd = smean(trans.hists.pvd(sitr:iter));
    avg_dvp = smean(trans.hists.dvp(sitr:iter));
    if avg_dvp > muNEWT.mu_delta
        trans.pmu = trans.pmu*(muNEWT.mu_fact);
    else
        trans.pmu = trans.pmu/(muNEWT.mu_fact);
    end
    trans.pmu = min(muNEWT.mu_max, max(muNEWT.mu_min, trans.pmu));
end
end




