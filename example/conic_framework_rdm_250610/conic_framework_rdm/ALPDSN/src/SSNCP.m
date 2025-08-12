%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-13 16:28:19
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%
function [X, out, y, S] = SSNCP(std_model, opts)

% SSNSDP: a semi-smooth Newton method for solving SDP.
%
%   (P) max <C,X> s.t. A(X)=b, X \in K,
%   (D) min b^Ty s.t.  C = At(y) + S,
%   where K is a positive semidefinite cone.
%
%%**********************************************************************
%

%% check input model
if ~isfield(std_model, 'At') || ~isfield(std_model, 'K') || ~isfield(std_model, 'C') || ~isfield(std_model, 'b')
    error("Input model is not valid");
end

%% make a copy
model = std_model;

At = model.At;
K = model.K;
b = model.b;
C = model.C;

%% Set the parameter
if nargin < 2; opts = []; end
opts = default_opts(opts);

if opts.log_path == ""
    fid = 1;
else
    fid = fopen(opts.log_path, 'w');
end

print_std_model_info(fid, model);

record = opts.record;
scale_data = opts.scale_data;
% trans.pmu = opts.pmu;% sigma
maxits = opts.maxits;
tol = opts.tol;
tolscale = opts.tolscale;

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
trans.K = K;
use_proximal = 0; % % better choice: use_proximal  = 1;
% modify_sig = 0;
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

trans.diagAAt = zeros(m, 1);

for p = 1:length(K)
    trans.diagAAt = trans.diagAAt + sum(At{p} .* At{p})';
end

trans.diagAAt = max(1e-4, full(trans.diagAAt));
trans.pmu = 1;

%% initial variable and operator
X.var = MatCell.zeros_like(C);
X.Avar = trans.Amap(X.var);
Z.var = MatCell.zeros_like(C);
Z.Avar = trans.Amap(Z.var);
W = Z.var;
y.var = zeros(trans.mdim, 1);
y.Avar = X.var;
S.var = MatCell.zeros_like(C);
S.Avar = trans.Amap(S.var);
% CS = cell(trans.nblock,1);

trans.bmAX = b - X.Avar;
trans.AC = trans.Amap(C);
trans.normb = max(1, norm(b));
trans.normC = max(1, norm(C));
% [~, X, trans] = comp_res_fp(Z,K,b,C,trans);

% trans = initial_ADMM_trans(trans,ADMMopts);
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

%% Compute projection
[F, S, trans] = compute_projection(K, X, y, b, C, trans);

trans = set_param(trans); % 做牛顿步设置参数

count = 0; % cg restart
%% main iteration
for iter = 1:maxits
    %             F.FY=FZnew.FY;
    %             F.FX=FZnew.FX;
    F.FXres = norm(F.FX);
    F.FYres = norm(F.FY);
    F.res = F.FXres + F.FYres;
    
    if use_proximal
        yhat = y.var;
        sighatmin = max(1, RpGradratio) * sighatmin; % % add: 27-Mar-2008
        sighatfac = max(1, RpGradratio ^ 2) * sighatfac;
        sighat = max(sighatmin, sighatfac * trans.pmu); % % add: 07-Apr-2008
        H = 1 ./ max(1, sqrt(abs(yhat))); % % better than H = ones(m,1);
        H2 = H .* H;
        %         F.par.yhat = yhat; F.par.sighat = sighat; F.par.H = H; F.par.H2 = H2;
    end
    
    if (use_proximal)
        trans.invdiagM = 1 ./ full(trans.pmu * trans.diagAAt + H2 / sighat);
    else
        trans.invdiagM = 1 ./ full(trans.pmu * trans.diagAAt);
    end
    
    trans.NEWT.iter = trans.NEWT.iter + 1;
    trans.NEWT.totaliter = trans.NEWT.totaliter + 1;
    % parameter for the revised Newton system
    trans.NEWT.tau1 = trans.NEWT.lambda * (F.FYres ^ trans.NEWT.sigPowy);
    trans.NEWT.tau2 = trans.NEWT.lambda * (F.FXres ^ trans.NEWT.sigPowx);
    %         trans.NEWT.sig = 1 ;
    if cgopts.CG_adapt %更改CG的tol
        cgopts.CG_tol = max(min(0.1 * F.FYres, 0.1), cgopts.cgtolmin) * 0.001; % F.res
    end
    
    %         [d,trans] = gendirection(F,K,At,trans,cgopts,NEWTopts,scale_data);
    % [dx, dy, trans] = gendirectionxy(F, K, At, trans, cgopts, 'iterative');
    [dx,dy,trans] = gendirectionxy(F, K, At, trans, cgopts, opts.method);
    for k = 1:trans.nblock
        if strcmp(K{p}.type, 's')
            dx{k} = (dx{k} + dx{k}') / 2;
        end
        
    end
    
    %----------------------------------------------------------------------
    ynew = struct();
    ynew.var = y.var + dy;
    Xnew = struct();
    ynew.Avar = trans.ATmap(ynew.var);
    Xnew.var = X.var + dx;
    [AX, AXorg] = trans.Amap(Xnew.var);
    Xnew.Avar = AX;
    Xnew.Avarorg = AXorg;
    
    [FZnew, S, trans] = compute_projection(K, Xnew, ynew, b, C, trans);
    
    trans.NEWT.nrmd = norm(dx) + norm(dy);
    trans.NEWT.FZd = -sum(FZnew.FY .* dy);
    
    for k = 1:trans.nblock
        trans.NEWT.FZd = trans.NEWT.FZd - sum(sum(FZnew.FX{k} .* dx{k}));
    end
    
    trans.NEWT.FZdorg = -sum(F.FY .* dy);
    
    for k = 1:trans.nblock
        trans.NEWT.FZdorg = trans.NEWT.FZdorg - sum(sum(F.FX{k} .* dx{k}));
    end
    
    thetaFZdorg = trans.NEWT.FZdorg / F.res / trans.NEWT.nrmd;
    trans.NEWT.FZFZnew = sum(F.FY .* FZnew.FY);
    
    for k = 1:trans.nblock
        trans.NEWT.FZFZnew = trans.NEWT.FZFZnew + sum(sum(F.FX{k} .* FZnew.FX{k}));
    end
    
    thetaFZFZnew = trans.NEWT.FZFZnew / FZnew.res / F.res;
    
    if F.res < 1e-3
        trans.NEWT.rhs = FZnew.res * trans.NEWT.nrmd ^ 2;
    else
        trans.NEWT.rhs = trans.NEWT.nrmd ^ 2;
    end
    
    trans.NEWT.ratio = trans.NEWT.FZd / trans.NEWT.rhs;
    
    % condition on performing a newton step 解的不好的时候
    if (trans.NEWT.iter == 1) && (FZnew.res >= F.res)
        trans.NEWT.lambda = trans.NEWT.lambda * trans.gfactor2;
        %             trans.NEWT.swt = trans.NEWT.swt + 1;
        continue;
    end
    
    if FZnew.res < 5e1 * F.res || count <= 5 %5e1*F.res
        count = 0;
        y = ynew;
        F = FZnew;
        X = Xnew;
        trans = sigPow_CGmaxit_update(F.res, trans);
        cgopts.CG_maxit = trans.NEWT.CG_maxit;
    else
        %             trans.NEWT.tau1 = trans.NEWT.lambda*(F.FYres^trans.NEWT.sigPowy);
        %             trans.NEWT.tau2 = trans.NEWT.lambda*(F.FXres^trans.NEWT.sigPowx);
        trans.NEWT.lambda = trans.gfactor * trans.NEWT.lambda;
        count = count +1;
        continue;
    end
    
    %------------------------------------------------------------------
    % check optimality
    trans = record_optimal_NEWT(C, b, trans, X, y, S, iter);
    
    if (max(trans.hists.pinforg(end), trans.hists.dinforg(end) * tolscale) < tol * retol)
        
        if opts.tolscale < 100
            [Xnew, ynew, Snew, rec] = recover_var_chk2nd(Z, F, X, y, S, W, trans);
            
            if max([rec.pinf, rec.dinf * tolscale, rec.K1, rec.K1dual, rec.C1]) < tol * retol
                cstop = 1;
                X = Xnew; y = ynew; S = Snew;
                trans.rec = rec;
            else
                cstop = 0;
            end
            
        else
            cstop = 1;
            [Xnew, ynew, Snew, rec] = recover_var_chk2nd(Z, F, X, y, S, W, trans);
            X = Xnew; y = ynew; S = Snew;
            trans.rec = rec;
        end
        
    end
    
    log_info = {'%5s', 'iter', num2str(iter);
        '%7s', 'mu', num2str(trans.pmu, '%2.1e');
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
        '%7s', 'nt_res', num2str(trans.newton_res, '%2.1e');
        '%5s', 'CGit', num2str(trans.cgiter, '%5d');
        % '%7s', 'ratio', num2str(trans.NEWT.ratio, '%2.1e')
        };
    
    %% define header for printing
    if record >= 1 && ~exist('str_head2', 'var')
        str_head2 = "";
        
        for i = 1:size(log_info, 1)
            str_head2 = str_head2 + sprintf(log_info{i, 1}, log_info{i, 2}) + " ";
        end
        
        fprintf(fid, '\n%s', str_head2);
    end
    
    %% print iteration information
    if record >= 1 && (cstop || ...
            iter == 1 || iter == maxits || mod(iter, NEWTopts.print_itr) == 0)
        
        if mod(iter, 20 * NEWTopts.print_itr) == 0 && iter ~= maxits && ~cstop
            fprintf(fid, '\n%s', str_head2);
        end
        
        str_info = "";
        
        for i = 1:size(log_info, 1)
            str_info = str_info + sprintf(log_info{i, 1}, log_info{i, 3}) + " ";
        end
        
        fprintf(fid, '\n%s', str_info);
        
    end
    
    %% 更新kappa 这里用的lambdba表示
    if trans.NEWT.ratio >= NEWTopts.eta2
        trans.NEWT.lambda = max(NEWTopts.gamma1 * trans.NEWT.lambda, 1e-16); % ratio越大说明越好，大的话要减小kappa
    elseif trans.NEWT.ratio >= NEWTopts.eta1
        trans.NEWT.lambda = NEWTopts.gamma2 * trans.NEWT.lambda;
    else
        trans.NEWT.lambda = NEWTopts.gamma3 * trans.NEWT.lambda;
    end
    
    if cstop
        out.status = sprintf('max(pinf,dinf) < %3.2e', tol);
        break;
    end
    
    if (iter == maxits)
        out.status = sprintf('reach the maximum iteration');
        [Xnew, ynew, Snew, rec] = recover_var_chk2nd(Z, F, X, y, S, W, trans);
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
            tolscale = tolscale / 2;
            if record == 1; fprintf(fid, 'change the tol\n'); end
        elseif trans.NEWT.totaliter == 750
            retol = 6;
            tolscale = tolscale / 2;
            if record == 1; fprintf(fid, 'change the tol\n'); end
        end
        
    end
    
    if muopts.adp_mu == 1
        pmup = trans.pmu;
        trans = mu_update_NEWT(trans, iter, muopts, opts);
        
        if trans.pmu ~= pmup
            if record >= 1; fprintf(fid, '\n  -- mu updated: %f\n', trans.pmu); end
            %                 [F, X, Z] = comp_res_fp_updmu(X, Z, K, b, C, trans.pmu,pmup,trans);
        end
        
        %             trans = mu_update_ADMM(trans,iter,muopts);
    end
    
end

out = generate_outs(fid, out, X, y, S, trans, opts, iter);

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
if ~isfield(NEWTopts, 'sigPow'); NEWTopts.sigPow = 1; end
if ~isfield(NEWTopts, 'resFac'); NEWTopts.resFac = 0.98; end
if ~isfield(NEWTopts, 'eta1'); NEWTopts.eta1 = 1e-4; end
if ~isfield(NEWTopts, 'eta2'); NEWTopts.eta2 = 0.9; end
if ~isfield(NEWTopts, 'gamma1'); NEWTopts.gamma1 = 0.5; end %可调
if ~isfield(NEWTopts, 'gamma2'); NEWTopts.gamma2 = 0.95; end %可调
if ~isfield(NEWTopts, 'gamma3'); NEWTopts.gamma3 = 10; end %可调
if ~isfield(NEWTopts, 'lambda'); NEWTopts.lambda = 1; end %adjust mu
%     if ~isfield(NEWTopts,'sigma'),    NEWTopts.sigma = 0.01;       end
%     if ~isfield(NEWTopts,'NTstep');   NEWTopts.NTstep  = 5;        end
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
if ~isfield(cgopts, 'cgtolmin'); cgopts.cgtolmin = 1e-10; end
if ~isfield(cgopts, 'CG_adapt'); cgopts.CG_adapt = 1; end
opts.cgopts = cgopts;

%% ------------------------------------------------------------------
% parameters for adjusting pmu
if ~isfield(opts, 'muopts'); opts.muopts = struct; end
muopts = opts.muopts;
if ~isfield(muopts, 'adp_mu'); muopts.adp_mu = 1; end %1 or 0

if ~isfield(muopts, 'NEWT'); muopts.NEWT = struct; end
% if ~isfield(muopts.NEWT, 'adpmu_cri');     muopts.NEWT.adpmu_cri = 1;      end %1 or 2
if ~isfield(muopts.NEWT, 'mu_min'); muopts.NEWT.mu_min = 1e-6; end
if ~isfield(muopts.NEWT, 'mu_max'); muopts.NEWT.mu_max = 1e6; end
if ~isfield(muopts.NEWT, 'mu_update_itr'); muopts.NEWT.mu_update_itr = opts.muopts.mu_update_itr; end %10 可调
if ~isfield(muopts.NEWT, 'mu_delta'); muopts.NEWT.mu_delta = 5e-1; end % 可调
if ~isfield(muopts.NEWT, 'mu_fact'); muopts.NEWT.mu_fact = opts.muopts.mu_fact; end %5/3 可调

% if ~isfield(muopts,'ADMM');          muopts.ADMM = struct; end
% if ~isfield(muopts.ADMM, 'ratio');   muopts.ADMM.ratio = 1.;      end %1 or 0
% if ~isfield(muopts.ADMM, 'mu_max');   muopts.ADMM.mu_max = 1e6;     end %1 or 2
% if ~isfield(muopts.ADMM, 'mu_min');   muopts.ADMM.mu_min = 1e-4;    end
% if ~isfield(muopts.ADMM, 'mu_delta');muopts.ADMM.mu_delta = 1.2;  end
% if ~isfield(muopts.ADMM, 'mu_fact'); muopts.ADMM.mu_fact = 1.8;  end

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
    At{k} = At{k}*DA; %这里是乘 后面算org还要除去
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

Lchol.isidentity = false;
b = fwsolve(Lchol,b);

if (opts.scale_data==1)
    bscale = max(1,norm(b));
    Cscale = max(1,norm(C));
end

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


function trans = initial_ADMM_trans(trans, ADMMopts)
trans.ADMM.rankS = zeros(trans.nblock, 1);
trans.ADMM.recompeig = 0;
% trans.ADMM.can_use_eigs = 0;
trans.ADMM.swt = 0;
trans.ADMM.maxiter = ADMMopts.imaxit;
trans.ADMM.subiter = 0;
trans.ADMM.iter = 0;
trans.muADMM.prim_win = 0;
trans.muADMM.dual_win = 0;
end

function trans = initial_NEWT_trans(trans, NEWTopts)
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

function [AX, AXorg] = AXmap(X, K, At, Lchol)
AXorg = AXfun(K, At, X); %原始的AX 在process_SDP函数那里A是经过除去对角的。
AX = fwsolve(Lchol, AXorg);
end

function Aty = Atymap(y, K, At, Lchol)
Aty = Atyfun(K, At, bwsolve(Lchol, y));
end

function [F, X, trans] = comp_res_fp(Z, K, b, C, trans)

pmu = trans.pmu;
[X.var, F.par] = project(K, Z.var);
[AX, AXorg] = trans.Amap(X.var);
X.Avar = AX;
X.Avarorg = AXorg;
mud = 2 * AX - Z.Avar - pmu * trans.AC - b;
F.var = trans.ATmap(mud) + pmu * C + Z.var - X.var;

F.res = norm(F.var);
F.Avar = AX - b;
F.par.precond = trans.precond;
end

function [F, S, trans] = compute_projection(K, X, y, b, C, trans)
tmp = C - X.var / trans.pmu;
W = trans.ATmap(y.var) - tmp;
[S0, F.par] = project(K, W);
F.par.precond = trans.precond; %加入预条件，为pre做准备

S.var = S0 - W;
S.Avar = trans.Amap(S.var);
F.FY = -b + trans.pmu * trans.Amap(S0);
F.FX = X.var / trans.pmu - S0;
F.res = norm(F.FY) + norm(F.FX);
end



function trans = record_optimal_NEWT(C, b, trans, X, y, S, iter)
% check optimality
% bTy  = full(b'*(trans.AC + (Z.Avar-X.Avar)/trans.pmu));
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

% S.var = ops_sdpnal(S.var,'*',trans.scale.Cscale);

y = y.var;
% S = S.var;
% scale = trans.scale;
y = trans.scale.Cscale * (trans.scale.DA * bwsolve(trans.Lchol, y));
% trans.AtymCSnrm = F.res/trans.pmu;
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

function [X, y, S, rec] = recover_var_chk2nd(Z, F, X, y, S, W, trans)
%     if strcmp(trans.last_iter, 'NEWT')
%         Zold = Z;
%         Z.var = ops_sdpnal(Z.var,'-',F.var);
%         Z.Avar = Z.Avar - F.Avar;
%         y = trans.AC + (Z.Avar - X.Avar)/trans.pmu;
%         S = ops_sdpnal(X.var,'-',Zold.var);
%         S = ops_sdpnal(S,'/',trans.pmu);
%         X = X.var;
%     else
y = y.var;
%         X = ops_sdpnal(ops_sdpnal(S.var,'-',W),'*',trans.pmu);
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
out.gap = rec.gap;
out.iter = iter;
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

function trans = sigPow_CGmaxit_update(res, trans)

if res < 1e-4
    sigPowy = trans.sigyl;
    sigPowx = trans.sigxl;
else
    sigPowy = trans.sigyu;
    sigPowx = trans.sigxu;
end

CG_maxit = trans.cgmax;
% if min(trans.cgres(end),trans.cgres(end-1))<1e-6 && trans.cgiter < CG_maxit
%     CG_maxit = trans.cgmin;
% end

trans.NEWT.sigPowx = sigPowx;
trans.NEWT.sigPowy = sigPowy;
trans.NEWT.CG_maxit = CG_maxit;
end

function trans = mu_update_NEWT(trans, iter, muopts, opts)
smean = @geo_mean;
muNEWT = muopts.NEWT;

if mod(iter, muNEWT.mu_update_itr) == 0
    sitr = iter - muNEWT.mu_update_itr + 1;
    avg_pvd = smean(trans.hists.pvd(sitr:iter));
    avg_dvp = smean(trans.hists.dvp(sitr:iter));
    
    if avg_dvp > muNEWT.mu_delta
        trans.pmu = trans.pmu * (muNEWT.mu_fact);
    else
        trans.pmu = trans.pmu / (muNEWT.mu_fact);
    end
    
    trans.pmu = min(muNEWT.mu_max, max(muNEWT.mu_min, trans.pmu));
end

end
