function [xopt, out] = SSNCVX(x0,pblk,Bt,f,Q,C,l,u,At,lb,ub,opts,y,z,v,r)
params.pflag = 1;
Ayes = 1;
method = 'iterative';
% 
nnls = 5;
params.fap = 0;
if isfield(opts,'tol');        tol           = opts.tol; end
if isfield(opts,'nnls');      nnls           = opts.nnls; end
if ~isfield(opts, 'nu');             opts.nu = 0.9999;  end
if ~isfield(opts, 'adaplambda');     opts.adaplambda = 0; end
if ~isfield(opts, 'lambda');     opts.lambda = 1; end
if ~isfield(opts, 'sigma');     opts.sigma = 1; end
if ~isfield(opts, 'resmin');     opts.resmin = 1e10; end
if isfield(opts, 'fap');     params.fap = 1; end

if isfield(opts,'method'); method    = opts.method; end

params.method = method;
params.Axeqb = 0;
params.fnonsmooth = 0;
params.nblock = length(pblk);
params.Ayes = Ayes;

params.l = l;
params.u = u;
params.lb = lb;
params.ub = ub;

%% Detect
if isfield(opts,'sva_path')
    params.save_path = opts.save_path;
    params.saveflag = 1;
else
    params.saveflag = 0;
end


if isempty(pblk)
    error('pblk is not exist');
else
    for j = 1:length(pblk)
        if ~isfield(pblk{j},'type')
            error('type of pblk does not exist');
        end
        if isfield(pblk{j},'cofficient')
            pblk{j}.cofficient = 1;
        end
        if isfield(pblk{j}.type,'fused')
            params.Binput = pblk{j}.Binput;
        end
    end
end
params.pblk = pblk;
params.DPhiP = @(iHW, d) DPhi_all(iHW, d, pblk);


%% inentify the type of the problem
if ~isempty(f)
    for i = 1:length(f)
        if ~isfield(f{i},'type')
            error('f exist but the type of f does not exist');
        else
            if strcmp(f{i}.type,'square') || strcmp(f{i}.type,'exp') || strcmp(f{i}.type,'logdet') || strcmp(f{i}.type,'log') || strcmp(f{i}.type,'logsumexp')
                params.fnonsmooth = 0;
            else
                params.fnonsmooth = 1;
            end
        end
        if isfield(f,'cofficient')
            f{i}.cofficient = 1;
        end
    end
    params.fflag = 1;
    params.f = f;
else
    params.fflag = 0;
end


if ~isempty(Q)
    params.Qflag = 1;
    if ismatrix(Q) && ~iscell(Q)
        params.Qmap = @(X) Q*X;
        params.Qcmap = @(X)  Q*X{1};
    elseif ismatrix(Q) && iscell(Q)
        params.Qcmap = @(X) Qmap(X, Q);
        params.Qmap = @(X) Qmap(X, Q');
    elseif class(Amap) == 'function_handle'
        params.Qmap = Q(X);
    end
else
    params.Qflag = 0;
    params.Qmap = @(X) 0;
    params.Qcmap = @(X) 0;
end




if xor(isempty(l),isempty(u))
    error('one box constraint is identified but the other is not');
else
    if ~(isempty(l) && isempty(u))
        params.boxflag = 1;
        params.idmap = @(X) idmap(X);
        params.isL = 1;
        if norm(l - u) == 0
            warning('the lb and ub of x is the same');
        end
    else
        params.isL = 0;
        params.boxflag = 0;
        params.idmap = @(X) 0;
    end
end

if xor(isempty(l),isempty(u))
    error('one box constraint is identified but the other is not');
else
    if ~(isempty(lb) && isempty(ub))
        params.Aboxflag = 1;
        if isempty(At)
            warning('Abox constraint is not exist, At is set as identity')
            At = [];
        end
        if norm(lb - ub,'fro') == 0
            params.Axeqb = 1;
            params.b = lb;
        end
    else
        params.Aboxflag = 0;
        At = [];
    end
end

if isempty(C)
    for i = 1:length(pblk)
        C{i} = zeros([pblk{i}.size,1]);
    end
end



%% Preprocess
params.Lchol = struct();
params.Lchol.isidentity = 1;
if nargin < 2; opts = []; end
opts = default_opts(opts);
NEWTopts = opts.NEWTopts;

if opts.log_path == ""
    fid = 1;
else
    fid = fopen(opts.log_path, 'w');
end
muopts = opts.muopts;
cgopts = opts.cgopts;
record = opts.record;
maxits = opts.maxits;
params.C = C;
if ~isempty(At)
    if ~strcmp(pblk{1}.type,'s')
        params.Amap = @(X) Amap(X, At);
        params.ATmap = @(y) ATmap(y, At);
    elseif norm(lb-ub) == 0
        b =lb;
        params.borg = b;
        [At,b,C, params] = preprocess_SDP(At,b,C, opts, params);
        params.Amap = @(X) AXmap(X, opts.K, At, params.Lchol);
        params.ATmap = @(y) Atymap(y, opts.K, At, params.Lchol);
        params.lb = b;
        params.ub = b;
        params.K = opts.K;
        params.b = b;
        params.At = At;
        params.C = C;
    end


else
    params.Amap = @(X) 0;
    params.ATmap = @(y) 0;
end
params.prox_p = @(X,params) prox_p(pblk, X, params);
params.pvalue = @(X,params) pvalue(pblk, X, params);
params.pdual_value = @(X,params) pdual_value(pblk, X, params);
if isfield(params,'f')
    params.fvalue = @(X,params) pvalue(params.f, X, params);
    params.fdualvalue = @(X,params) pdual_value(params.f, X, params);
else
    params.fvalue = @(X,params) 0;
    params.fdualvalue = @(X,params) 0;
end


if ~isempty(f) && params.fnonsmooth
    params.prox_f = @(X,params) prox_p(f, X, params);
elseif ~isempty(f) && ~params.fnonsmooth
    params.dual_gradientf = @(X) dual_gradientf(f,X,params);
    params.fhess = @(X) 1;
end

if isempty(Bt) && isempty(f)
    params.Bmap = @(X) Zeromap(X);
    params.BTmap = @(y) 0;
    for i = 1:params.nblock
        params.Bt{i} = speye(sum(pblk{i}.size));
    end
elseif isempty(Bt) && ~isempty(f)
    params.Bmap = @(X) idmap(X);
    params.BTmap = @(y) idmap(y);
    for i = 1:params.nblock
        params.Bt{i,1} = speye(pblk{i}.size);
    end
elseif isa(Bt,'double')

    
    params.Bmap = @(X) Amap(X, Bt);
    params.BTmap = @(X) ATmap(X, Bt);
    params.Bt = {Bt};
    params.Btd = Bt;
    params.B = Bt';
elseif  isa(Bt,'cell')
    params.Bmap = @(X) Amap(X, Bt);
    params.BTmap = @(X) ATmap(X, Bt);
elseif isa(Bt,"struct")
    if isfield(Bt,'Bmap') && isfield(Bt,'BTmap')
        params.Bmap = @(X) Bt.Bmap(X);
        params.BTmap = @(X) Bt.BTmap(X);
    else
        params.Bmap = @(X) idmap(X);
        params.BTmap = @(y) idmap(y);
    end
end
params.mdim = size(At,2);

params.tic = tic;
for i =  1:length(pblk)
    if ~strcmp(pblk{i}.type,'s')
        if iscell(At)
        params.At{i} = At{i};
        else
        params.At{i} = At;
        end
    end
    params.Qt{i} = Q;
end

params.proj_AX = @(y) projBox( y, params.lb,params.ub);
if isfield(params,'xL')
    % params.proj_X = @(X) projBox( X, params.xL,params.xU);
    params.P1box = @(X) projBox( X, params.xL,params.xU);
else
% params.proj_X = @(X) projBox( X, params.l,params.u);
params.P1box = @(X) projBox(X,l,u);
end
% params.P1box = @(X) projBox(X,l,u);
params.P2box = @(y) projBox(y,lb,ub);
params.sigma = opts.sigma;

params.cgmin = opts.cgmin;
params.cgmed = opts.cgmed;
params.cgmax = opts.cgmax;



params.sigyl = opts.sigyl;
params.sigym = opts.sigym;
params.sigyu = opts.sigyu;
params.sigzl = opts.sigzl ;
params.sigzm =    opts.sigzm ;
params.sigzu =    opts.sigzu ;
params.sigrl =    opts.sigrl ;
params.sigrm =   opts.sigrm;
params.sigru =    opts.sigru ;
params.sigvl =    opts.sigvl ;
params.sigvm =   opts.sigvm;
params.sigvu =    opts.sigvu ;

params.sigx1l = opts.sigx1l;
params.sigx1m = opts.sigx1m;
params.sigx1u = opts.sigx1u;
params.sigx2l = opts.sigx2l;
params.sigx2m = opts.sigx2m;
params.sigx2u = opts.sigx2u;
params.sigx3l = opts.sigx3l;
params.sigx3m = opts.sigx3m;
params.sigx3u = opts.sigx3u;
params.sigx4l = opts.sigx4l;
params.sigx4m = opts.sigx4m;
params.sigx4u = opts.sigx4u;
str = [opts.basename, newline];

params = initial_NEWT_params(params, NEWTopts);
params.NEWT.lambda = opts.lambda;
%% initial variable and operator
if params.Aboxflag == 1
    params.mdim = opts.m;
    y.var = zeros(params.mdim, 1);
else
    y.var = {zeros(1,1)};
end
params.y0 = y.var;
if ~isempty(f)
    if isempty(Bt)
        if strcmp(params.f{1}.type,'exp')
        Z.var = -1e-7*MatCell.eyes_like(full(C));
        elseif  strcmp(params.f{1}.type,'logdet')
            Z.var = 1e-7*MatCell.eyes_like(full(C));
        else
        Z.var = 0*MatCell.eyes_like(full(C));
        end
        X2.var = MatCell.zeros_like(Z.var);
    elseif isa(Bt,"struct")
        Z.var = {zeros([Bt.out_size,1])};
        Z.var = MatCell.zeros_like(Z.var);
        X2.var = MatCell.zeros_like(C);
    elseif isa(Bt,"double")
        Z.var = {zeros([size(Bt,2),1])};
        Z.var = MatCell.zeros_like(Z.var);
        X2.var = MatCell.zeros_like(Z.var);
        % X2.Avar = params.Amap(X2.var);
    elseif isa(Bt,"cell")
        for i = 1:length(Bt)
            Z.var{i,1} = zeros([size(Bt{i},2),1]);
        end
        Z.var = MatCell.zeros_like(Z.var );
        X2.var = MatCell.zeros_like(Z.var);
    end
    params.z0 = Z.var;
else
    for i = 1:params.nblock
        Z.var{i,1} = zeros(1,1);
        params.z0{i,1} = zeros(1,1);
    end
    X2.var = MatCell.zeros_like(C);
end

if params.boxflag == 1
    R.var = MatCell.zeros_like(C);
else
    for i = 1:params.nblock
    R.var{i,1} = zeros(1,1);
    end
end
params.r0 = R.var;
if params.Qflag == 1
    V.var = MatCell.zeros_like(C);
    V.Avar = params.Amap(V.var);
else
    for i = 1:params.nblock
    V.var{i,1} = zeros(1,1);
    end
end
params.v0 = V.var;
% params.hists.pobj = 1;
% params.hists.dobj = 1;
% params.hists.gap = 1;

S.var = MatCell.zeros_like(C);
S.Avar = params.Amap(S.var);

X1.var = {zeros(size(y.var))};

X3.var = MatCell.zeros_like(C);
X3.Avar = params.Amap(X3.var);
X4.var = 0*MatCell.eyes_like(C);
X4.Avar = params.Amap(X4.var);

params.xsize = size(X4.var);
params.cgall = 0;
params.AC = params.Amap(C);
if isfield(params,'b')
    params.normb = max(1, norm(params.b));
end
params.normC = max(1, norm(C));
params.FYtmp = 0;
params.FXtmp = 0;
params.NEWT.iter = 0;
params.res = 1;
cstop = 0;
params.cgtol = opts.cgtol;
params.cgratio = opts.cgratio;
params.resratio = opts.resratio;
params.tstart = tic;
[init] = init_parameters(params);
xzi_x = zeros(nnls, 1);
[F,~,~,params] = compute_gradient(y,Z,R,V,X1,X2,X3,X4,params);
nFxz = max(norm(F.Fzres, 'fro'), norm(F.Fx4res, 'fro'));
xzi_x(1) = nFxz;
count = 0;
%% main solver
for iter = 1:opts.maxits


    params.NEWT.iter = params.NEWT.iter + 1;

    if cgopts.CG_adapt
        cgopts.CG_tol = max(min(0.1 * F.res, 0.1)*params.cgratio, params.cgtol)  ; % F.res
    end

    params.NEWT.tau1 = init.tau1*params.NEWT.lambda*(F.Fyres^params.NEWT.sigPowy);
    % params.NEWT.tau2 = init.tau2*params.NEWT.lambda*(F.Fzres^params.NEWT.sigPowz) + 1e-7;
    params.NEWT.tau2 =  params.NEWT.lambda*init.tau2*max(min((F.Fzres^params.NEWT.sigPowz),opts.resmin),1e-5);
    params.NEWT.tau3 =  init.tau3*params.NEWT.lambda*(F.Frres^params.NEWT.sigPowr)+ 1e-7;
    params.NEWT.tau4 = init.tau4*params.NEWT.lambda*(F.Fvres)^params.NEWT.sigPowv + 1e-7;

    params.NEWT.taux1 = init.taux1*params.NEWT.lambda*(F.Fx1res^params.NEWT.sigPowx1)+1e-7;
    params.NEWT.taux2 = init.taux2*params.NEWT.lambda*(F.Fx2res^params.NEWT.sigPowx2)+ 1e-7;
    params.NEWT.taux3 = init.taux3*params.NEWT.lambda*(F.Fx3res^params.NEWT.sigPowx3)+ 1e-7;
    params.NEWT.taux4 = params.NEWT.lambda*init.taux4*max(min((F.Fx4res)^params.NEWT.sigPowx4,opts.resmin),1e-5);

    if iter == 13
        % break;
        1;
    end

    [out,params] = getdirectionall(F,params,cgopts);
    params.cgall = params.cgall + params.cgiter;



    ynew.var = y.var + out.dy;
    ynew.Avar = params.ATmap(ynew.var);

    Znew.var = Z.var + out.dz;
    Rnew.var = R.var + out.dr;
    Vnew.var = V.var + out.dv;
    X1new.var = X1.var + out.dx1;
    X2new.var = X2.var + out.dx2;

    if isfield(params.D,'D3')
        X3new.var = X3.var + out.dx3;
    else
        X3new.var = [];
    end

    X4new.var = X4.var + out.dx4;
    [Fnew,Snew,X4new,params] = compute_gradient(ynew,Znew,Rnew,Vnew,X1new,X2new,X3new,X4new,params);

    


    

    [AX] = params.Amap(X4new.var);
    X4new.Avar = AX;


    %----------------------------------------------------------------------

    params.NEWT.nrmd = norm(out.dy)^2 + norm(out.dz)^2 + norm(out.dr)^2 + norm(out.dv)^2 +  norm(out.dx4)^2;
    if params.Aboxflag == 1 && params.Axeqb == 0
        params.NEWT.nrmd = params.NEWT.nrmd + norm(out.dx1)^2;
    end

    if params.fnonsmooth
        params.NEWT.nrmd = params.NEWT.nrmd + norm(out.dx2)^2;
    end

    if params.boxflag == 1
        params.NEWT.nrmd = params.NEWT.nrmd + norm(out.dx3)^2;
    end

    params.NEWT.nrmd = sqrt(params.NEWT.nrmd);
    params.NEWT.Fd = -dot_ssn(Fnew.Fy, out.dy);

    for k = 1:params.nblock
        params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fz{k} , out.dz{k});
        params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fv{k}, out.dv{k});
        params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fr{k}, out.dr{k});

        if params.Aboxflag == 1 && params.Axeqb == 0
            params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fx1{k}, out.dx1{k});
        end
        if params.fnonsmooth
            params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fx2{k}, out.dx2{k});
        end
        if params.boxflag == 1
            params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fx3{k}, out.dx3{k});
        end
        params.NEWT.Fd = params.NEWT.Fd - dot_ssn(Fnew.Fx4{k}, out.dx4{k});
    end


    if F.res < params.resratio
        params.NEWT.rhs = 1e-1*sqrt(Fnew.res) * params.NEWT.nrmd ^ 2;
    else
        params.NEWT.rhs = params.NEWT.nrmd ^ 2;
    end

    params.NEWT.ratio = params.NEWT.Fd / params.NEWT.rhs;
    if opts.adaplambda
        if params.NEWT.ratio >= NEWTopts.eta2
            params.NEWT.lambda = max(NEWTopts.gamma1 * params.NEWT.lambda, 1e-16);
        elseif params.NEWT.ratio >= NEWTopts.eta1
            params.NEWT.lambda = NEWTopts.gamma2 * params.NEWT.lambda;
        else
            params.NEWT.lambda = NEWTopts.gamma3 * params.NEWT.lambda;
        end
    end

    if ~isempty(f) && strcmp(params.f{1}.type,'square') && (strcmp(params.pblk{1}.type,'l1') || strcmp(params.pblk{1}.type,'fused'))
        nFu = max(Fnew.Fzres, Fnew.Fx4res);
        iter_ns = 1;
        while nFu > 0.999 * max(xzi_x) && iter_ns < 2
            Znew.var = Z.var + opts.linratio^iter_ns*out.dz;
            X4new.var = X4.var + opts.linratio^iter_ns*out.dx4;
            [Fnew,Snew,X4new,params] = compute_gradient(ynew,Znew,Rnew,Vnew,X1new,X2new,X3new,X4new,params);
            nFu1 = Fnew.Fzres; nFu2 = Fnew.Fx4res;
            nFu = max(nFu1, nFu2);
            iter_ns = iter_ns +1;
        end
        xzi_x(mod(iter-1, nnls) + 1) = nFu;

        if (nFu < opts.nu * max(xzi_x)  || count >5) || iter>0
            params.NEWT.lambda = params.NEWT.lambda * opts.lfactor;
        else
            % lambda = lambda * 1.1;
            % count = count+1;
            % status = 'o';
        end
    end




    if Fnew.res < 5e2 * F.res || 1
        count = 0;
        y = ynew;
        F = Fnew;
        Z = Znew;
        S = Snew;
        V = Vnew;
        R = Rnew;
        if params.Aboxflag == 1 && params.Axeqb == 0
            X1 = X1new;
        end

        if params.fnonsmooth
            X2 = X2new;
        end

        if params.boxflag == 1
            X3 = X3new;
        end

        X4 = X4new;
        if iter > 2
            params = sigPow_CGmaxit_updateplus(F.res,params,cgopts);
            cgopts.CG_maxit = params.NEWT.CG_maxit;
        end
    elseif count > 5
        params.NEWT.lambda = iter;
        count = 0;
    else
        params.NEWT.lambda = params.gfactor * params.NEWT.lambda;
        count = count +1;

        %% print...   to be implemented

        continue;
    end
    %------------------------------------------------------------------
    % check optimality
    % if iter == 1 || xor(~isempty(f) && strcmp(params.f{1}.type,'square') && (strcmp(params.pblk{1}.type,'l1') || strcmp(params.pblk{1}.type,'fused')), F.Fx4res > 1e-5)
    params = record_optimal_NEWT(y,Z,R,V,S,X4,params,iter);
    % end
    %% define header for printing
    log_info = {'%5s', 'iter', num2str(iter);
        '%7s', 'sigma', num2str(params.sigma, '%2.1e');
        '%8s', 'pobj', num2str(params.hists.pobj(end), '%2.1e');
        '%8s', 'dobj', num2str(params.hists.dobj(end), '%2.1e');
        '%7s', 'gap', num2str(params.hists.gap(end), '%2.1e');
        % '%7s', 'gaporg', num2str(params.hists.gaporg(end), '%2.1e');
        '%7s', 'pinf', num2str(params.hists.pinf(end), '%2.1e');
        '%7s', 'dinf', num2str(params.hists.dinf(end), '%2.1e');
        '%7s', 'F', num2str(F.res, '%2.1e');
        '%7s', 'tau1', num2str(params.NEWT.tau1, '%2.1e');
        '%7s', 'tau2', num2str(params.NEWT.tau2, '%2.1e');
        '%7s', 'tau3', num2str(params.NEWT.tau3, '%2.1e');
        '%7s', 'tau4', num2str(params.NEWT.tau4, '%2.1e');
        '%7s', 'CGres', num2str(params.cgres(end), '%2.1e');
        %                     '%7s', 'nt_res', num2str(params.newton_res, '%2.1e');
        '%5s', 'CGit', num2str(params.cgiter, '%5d');
        %         '%7s', 'ratio', num2str(params.NEWT.ratio, '%2.1e')
        };


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


    if (max([params.hists.pinf(end),params.hists.K1,params.hists.K2,params.hists.P1, params.hists.dinf(end)]) < tol) && iter > 10
        if opts.tolscale < 100
            [ Xnew2, ynew2, Znew2,rec] = recover_var_chk2nd(y,Z,R,V,S,X4,params,iter);
            if max([rec.pinf,rec.dinf]) < tol
                cstop = 1;
                y = ynew2;  Z= Znew2;
                X= Xnew2;
                params.rec = rec;
            else
                cstop = 0;
            end
        else
            cstop = 1;
            [ Xnew2, ynew2, Znew2,rec] = recover_var_chk2nd(y,Z,V,R,S,X4,params);
            y = ynew2; R = Rnew2; Z= Znew2; V = Vnew2;
            X = Xnew2;
            params.rec = rec;
        end

    end


    if cstop
        out.status = sprintf('max(pinf,dinf) < %3.2e', tol);
        break;
    end

    if (iter == maxits) || toc(params.tic) > 10000
        out.status = sprintf('reach the maximum iteration');
        [ynew, Znew, Vnew, Rnew, X1new, X2new, X3new, X4new,rec] = recover_var_chk2ndplus(y,Z,V,R,X1,X2,X3,X4,params);
        y = ynew; R = Rnew; Z= Znew; V = Vnew;
        X1 = X1new; X2 = X2new; X3 = X3new; X4 = X4new;
        params.rec = rec;
        break;
    end



    if muopts.adp_mu == 1
        pmup = params.sigma;
        params = mu_update_NEWT(params, iter, muopts, opts);

        if params.sigma ~= pmup
            if record >= 1; fprintf(fid, '\n  -- mu updated: %f\n', params.sigma); end
        end

    end

end
params.str = str;
if params.saveflag == 1
    fid = fopen(params.save_path,'a+');
    fprintf(fid,'%s',params.str);
else
    fid = 0;
end
out = generate_outs(fid, out, params, opts, iter);
% out = [];

out.cgall = params.cgall;
xopt = X4;
end



%%
function params = record_optimal_NEWT(y,Z,R,V,S,X4,params,iter)
trCX = full(dot_ssn(full(params.C), X4.var));
QX4 = params.Qcmap(X4.var);
Qv = params.Qcmap(V.var);
xQx = 0.5*dot_ssn(QX4,X4.var);
if params.fap == 1
    for i = 1:params.nblock
        R.var{i}(abs(R.var{i})<5e-7) = 0;
        XP{i,1} = max(-R.var{i}*params.scale.Cscale,0);
        XN{i,1} = min(-R.var{i}*params.scale.Cscale,0);
    end    
else
    for i = 1:params.nblock
        XP{i,1} = max(-R.var{i},0);
        XN{i,1} = min(-R.var{i},0);
    end
end

if params.Aboxflag == 1
    if params.fap == 1
        ytmp = params.scale.Cscale*(params.scale.DA*bwsolve(params.Lchol,y.var));
    else
        ytmp = y.var;
    end
    yP = max(-ytmp,0);
    yN = min(-ytmp,0);
else
    yP = 0;
    yN = 0;
end


tmp2 = params.Bmap(X4.var);
if isfield(params,'fvalue')
    if params.fap == 1
    pobj =  params.scale.objscale*(xQx + trCX + params.fvalue(tmp2,params) + params.pvalue(X4.var,params));
    else
    pobj =  (xQx + trCX + params.fvalue(tmp2,params) + params.pvalue(X4.var,params));
    end
else
    pobj = xQx + trCX  + params.pvalue(X4.var,params);
end
if params.boxflag == 1
    if params.fap == 1
    tmpl{1} = max(params.l{1},-1e5);
    tmpu{1} = min(params.u{1},1e5);
    else
    tmpl = params.l;
    tmpu = params.u;
    end
    dualtmp1 = dot_ssn(XN, tmpl) + dot_ssn(XP,tmpu); % sum(XP.*params.u) ;
else
    dualtmp1 = 0;
end

if params.Aboxflag == 1
    if params.fap == 1
    dualtmp2 = sum(yN.*params.borg) + sum(yP.*params.borg);
    else
    dualtmp2 = sum(yN.*params.lb) + sum(yP.*params.ub);
    end
else
    dualtmp2 = 0;
end

if isfield(params,'f')
    dobj = -(dualtmp2 + dualtmp1 + 0.5*dot_ssn(Qv,X4.var) + params.fdualvalue(-Z.var,params) + params.pdual_value(-S.var,params) );
    % dobj = params.scale.objscale*dobj;
else
    dobj = -(dualtmp2 + dualtmp1 + 0.5*dot_ssn(X4.var{1},Qv)) + params.pdual_value(-S.var,params);
end

gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));

normX = norm(X4.var);
normS = norm(S.var);
normR = norm(R.var);
normZ = norm(Z.var);

pinf = norm(X4.Avar - params.proj_AX(X4.Avar - y.var))/(1 + normX);
y = y.var;
dinf = norm(params.ATmap(y) + params.BTmap(Z.var) + params.idmap(R.var) +  S.var - params.Qcmap(V.var) -params.C)/(1 + norm(params.C));

sigmatmp = params.sigma;
params.sigma = 1;
etaK1 = 0;

for i = 1:params.nblock
    if isfield(params.pblk{i},'shift')
        etaK1 = norm(X4.var - params.pblk{i}.shift - params.prox_p(X4.var -params.pblk{i}.shift - S.var,params))/(1 + normX + normS);
    else
        etaK1 = norm(X4.var - params.prox_p(X4.var- S.var,params))/(1 + normX + normS);
    end
end
params.sigma = sigmatmp;


QX4nrm = norm(QX4);
Qvnrm = norm(Qv);

etaK2 = norm(QX4 - Qv)/(1 + QX4nrm +  Qvnrm);
etaP1 = 0;
if params.fflag && params.fnonsmooth
    tmp = params.Bmap(X4.var);
    for p =1:params.nblock
        if isfield(params.f{p},'shift')
            tmp = tmp - params.f{p}.shift ;
        end
    end
    etaP1 = norm(params.prox_f(tmp - Z.var,params) - tmp);
elseif params.fflag && ~params.fnonsmooth
    etaP1 = norm(params.dual_gradientf(Z.var) - params.Bmap(X4.var ))/(1 + normX + normZ);
    params.hists.P1 = etaP1;
end

if params.boxflag == 1
    etaP2 = norm(params.P1box(X4.var - R.var) - X4.var)/(1 + normX + normR);
else
    etaP2 = 0;
end

params.hists.pobj(iter) = pobj;
params.hists.dobj(iter) = dobj;
params.hists.gap(iter) = gap;
% params.hists.gap(iter) = gap;
params.hists.pinf(iter) = pinf;
% params.hists.pinforg(iter) = pinforg;
params.hists.dinf(iter) = dinf;
% params.hists.dinforg(iter) = dinforg;
params.hists.pvd(iter) = pinf / dinf;
params.hists.dvp(iter) = dinf / pinf;
% params.hists.pvd(iter) = pinforg / dinf;
% params.hists.dvg(iter) = dinf / pinforg;
params.hists.cgiter(iter) = params.cgiter;
params.hists.K1 = etaK1;
params.hists.K2 = etaK2;

params.hists.P2 = etaP2;
params.hists.P1 = etaP1;
% params.hists.P1org = etaP1org;
% params.hists.P1dual = etaP1dual;
% params.hists.P1dualorg = etaP1dualorg;
maxinf = max(pinf, dinf);

params.hists.isNEWT(iter) = 1;
params.hists.maxinf = maxinf;

% params.bmAX = bmAX;
end


function [X, y,  Z,rec] = recover_var_chk2nd(y,Z,R,V,S,X,params,iter)

trCX = full(dot_ssn(full(params.C), X.var));
QX4 = params.Qcmap(X.var);
Qv = params.Qcmap(V.var);
xQx = 0.5*dot_ssn(QX4,X.var);
for i = 1:params.nblock
    XP{i,1} = max(-R.var{i},0);
    XN{i,1} = min(-R.var{i},0);
end
if params.Aboxflag == 1
    yP = max(-y.var,0);
    yN = min(-y.var,0);
else
    yP = 0;
    yN = 0;
end
% dobj = params.scale.objscale * bTy;
% pobj = params.scale.objscale * trCX ;

tmp2 = params.Bmap(X.var);
if isfield(params,'fvalue')
    pobj = xQx + trCX + params.fvalue(tmp2,params) + params.pvalue(X.var,params);
else
    pobj = xQx + trCX  + params.pvalue(X4.var,params);
end
% gap = abs(trCX - bTy) / max(1, abs(trCX));
% gaporg = abs(pobj - dobj) / max(1, abs(pobj));
if params.boxflag == 1
    dualtmp1 = dot_ssn(XN,params.l) + dot_ssn(yN,params.lb);
    dualtmp2 = dot_ssn(XP,params.u) + dot_ssn(yP,params.ub);
else
    dualtmp1 = 0;
    dualtmp2 = 0;
end

if isfield(params,'f')
    dobj = -(dualtmp2 + dualtmp1 + 0.5*dot_ssn(Qv,X.var) + params.fdualvalue(-Z.var,params) + params.pdual_value(-S.var,params) );
else
    dobj = -(dualtmp2 + dualtmp1 + 0.5*dot_ssn(X.var{1},Qv)) + params.pdual_value(-S.var,params);
end

gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
% bmAX = b - X.Avar;
% bmAXnrm = norm(X.Avarorg ./ diag(params.scale.DA) * params.scale.bscale - params.borg);

% pinforg = bmAXnrm / (1 + norm(params.borg));
% pinf = norm(bmAX) / params.normb;

normX = norm(X.var);
normS = norm(S.var);
normR = norm(R.var);
normZ = norm(Z.var);

pinf = norm(X.Avar - params.proj_AX(X.Avar - y.var))/(1 + normX);
y = y.var;
dinf = norm(params.ATmap(y) + params.BTmap(Z.var) + params.idmap(R.var) +  S.var - params.Qcmap(V.var) -params.C)/(1 + norm(params.C));

sigmatmp = params.sigma;
params.sigma = 1;
etaK1 = norm(X.var - params.prox_p(X.var- S.var,params))/(1 + normX + normS);
params.sigma = sigmatmp;


QX4nrm = norm(QX4);
Qvnrm = norm(Qv);

etaK2 = norm(QX4 - Qv)/(1 + QX4nrm +  Qvnrm);
etaP1 = 0;
if params.fflag && params.fnonsmooth
    tmp = params.Bmap(X.var);
    etaP1 = norm(params.prox_f(tmp - Z.var,params) - tmp);
    params.hists.P1 = etaP1;
elseif params.fflag && ~params.fnonsmooth
    etaP1 = norm(-params.dual_gradientf(Z.var) - params.Bmap(X.var))/(1 + normX + normZ);
    params.hists.P1 = etaP1;
end

if params.boxflag == 1
    etaP2 = norm(params.P1box(X.var - R.var) - X.var)/(1 + normX + normR);
else
    etaP2 = 0;
end

% PX = norm(params.projectP(X.var) - X.var);
% PZb = norm(params.projectP(- Zb.var));
% etaP1 = PX/params.normC;
% etaP1dual = PZb/params.normC;
% etaP1org = PX/(1+norm(params.Corg));
% etaP1dualorg = PZb/(1+norm(params.Corg));

% params.hists.pobj(iter) = pobj;
% params.hists.dobj(iter) = dobj;
% params.hists.gaporg(iter) = gaporg;
% params.hists.gap(iter) = gap;
params.hists.pinf(iter) = pinf;
% params.hists.pinforg(iter) = pinforg;
params.hists.dinf(iter) = dinf;
% params.hists.dinforg(iter) = dinforg;
params.hists.pvd(iter) = pinf / dinf;
params.hists.dvp(iter) = dinf / pinf;
% params.hists.pvd(iter) = pinforg / dinf;
% params.hists.dvg(iter) = dinf / pinforg;
params.hists.cgiter(iter) = params.cgiter;
params.hists.K1 = etaK1;
params.hists.K2 = etaK2;

params.hists.P2 = etaP2;
% params.hists.P1org = etaP1org;
% params.hists.P1dual = etaP1dual;
% params.hists.P1dualorg = etaP1dualorg;
maxinf = max(pinf, dinf);

params.hists.isNEWT(iter) = 1;
params.hists.maxinf = maxinf;

% params.bmAX = bmAX;

rec.pinf = pinf;
rec.dinf = dinf;
rec.gap = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
rec.pobj = pobj;
rec.dobj = dobj;
rec.totaltime = toc(params.tstart);
rec.P1 = etaP1;
rec.P2 = etaP2;

rec.K1 = etaK1;
rec.K2 = etaK2;


end

function params = sigPow_CGmaxit_updateplus(res,params,cgopts)
if res < 1e-4
    sigPowy      = params.sigyl;
    sigPowz     = params.sigzl;
    sigPowv     = params.sigvl;
    sigPowr     = params.sigrl;

    sigPowx1      =params.sigx1l;
    sigPowx2      =params.sigx2l;
    sigPowx3      =params.sigx3l;
    sigPowx4      =params.sigx4l;
elseif res < 1e-2
    sigPowy      = params.sigym;
    sigPowz     = params.sigzm;
    sigPowv     = params.sigvm;
    sigPowr     = params.sigrm;

    sigPowx1      =params.sigx1m;
    sigPowx2      =params.sigx2m;
    sigPowx3      =params.sigx3m;
    sigPowx4      =params.sigx4m;
else
    sigPowy      = params.sigyu;
    sigPowz     = params.sigzu;
    sigPowv     = params.sigvu;
    sigPowr     = params.sigru;

    sigPowx1      =params.sigx1u;
    sigPowx2      =params.sigx2u;
    sigPowx3      =params.sigx3u;
    sigPowx4      =params.sigx4u;
end



% CG_maxit = params.cgmax;
if params.NEWT.CG_maxit == params.cgmin  && min(params.cgres(end),params.cgres(end-1))>cgopts.CG_tol && params.cgiter > params.cgmin
    params.NEWT.CG_maxit = params.cgmed;
elseif params.NEWT.CG_maxit == params.cgmed  && params.cgres(end)>max(cgopts.CG_tol,1e-6) && params.cgiter > params.cgmed
    params.NEWT.CG_maxit = params.cgmax;
elseif params.NEWT.CG_maxit == params.cgmax  && params.cgres(end)>max(cgopts.CG_tol,1e-6) && params.cgiter > params.cgmax
    params.NEWT.CG_maxit = params.cgmin;
end


params.NEWT.sigPowy = sigPowy;
params.NEWT.sigPowz = sigPowz;
params.NEWT.sigPowr = sigPowr;
params.NEWT.sigPowv = sigPowv;

params.NEWT.sigPowx1 = sigPowx1;
params.NEWT.sigPowx2 = sigPowx2;
params.NEWT.sigPowx3 = sigPowx3;
params.NEWT.sigPowx4 = sigPowx4;
end


function out = generate_outs(fid, out, params, opts, iter)

rec = params.rec;

if params.saveflag
    fprintf(fid, '\n-------------------------------------------------------------------------------------------------\n') ;
    fprintf(fid, '%14s %14s %10s %10s %10s %10s %10s %12s\n', 'pboj', 'dobj', 'gap', 'pinf', 'dinf', 'C1', 'K1', 'K2');
    fprintf(fid, '%14.8e %14.8e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n', ...
        rec.pobj, rec.dobj, rec.gap, rec.pinf, rec.dinf, rec.C1, rec.K1, rec.K2) ;
    fprintf(fid, '\n-------------------------------------------------------------------------------------------------\n') ;
end

params.hists.pobj(iter + 2) = rec.pobj;
params.hists.dobj(iter + 2) = rec.dobj;
params.hists.pinf(iter + 2) = rec.pinf;
params.hists.pinforg(iter + 2) = rec.pinf;
params.hists.dinf(iter + 2) = rec.dinf;
params.hists.dinforg(iter + 2) = rec.dinf;
params.hists.gap(iter + 2) = rec.gap;
params.hists.gaporg(iter + 2) = rec.gap;
params.hists.pvd(iter + 2) = rec.pinf / rec.dinf;
params.hists.dvp(iter + 2) = rec.dinf / rec.pinf;
params.hists.pvdorg(iter + 2) = rec.pinf / rec.dinf;
params.hists.dvporg(iter + 2) = rec.dinf / rec.pinf;
%params.hists.res(iter+2) = res;

out.pobj = rec.pobj;
out.dobj = rec.dobj;
out.gap  = rec.gap;
out.iter  = iter;
out.P1 = rec.P1;

out.K1 = rec.K1;
out.K2 = rec.K2;

% if params.isL
%     out.C2 = rec.C2;
%     out.P1 = rec.P1;
%     out.P1dual = rec.P1dual;
% else
%     out.Knorm = rec.Knorm;
% end
out.pinf = rec.pinf;
out.dinf = rec.dinf;
out.rec = rec;
out.hists = params.hists;
out.totaltime = rec.totaltime;

end

function params = initial_NEWT_params(params, NEWTopts)
params.NEWT.lambda = NEWTopts.lambda;

params.NEWT.sigPowy = params.sigyu;
params.NEWT.sigPowz = params.sigzu;
params.NEWT.sigPowr = params.sigru;
params.NEWT.sigPowv = params.sigvu;

params.NEWT.sigPowx1 = params.sigx1u;
params.NEWT.sigPowx2 = params.sigx2u;
params.NEWT.sigPowx3 = params.sigx3u;
params.NEWT.sigPowx4 = params.sigx4u;

params.NEWT.iter = 0;
params.NEWT.swt = 0;
params.NEWT.CG_maxit = params.cgmin;
params.NEWT.subiter = 0;
params.NEWT.maxiter = 20;
params.NEWT.lastres = inf;
end

function params = mu_update_NEWT(params, iter, muopts, opts)
smean = @geo_mean;
muNEWT = muopts.NEWT;

if mod(iter, muNEWT.mu_update_itr) == 0
    sitr = iter - muNEWT.mu_update_itr + 1;
    avg_pvd = smean(params.hists.pvd(sitr:iter));
    avg_dvp = smean(params.hists.dvp(sitr:iter));
    if iter < 10000000000000
        avg_dvp = smean(params.hists.dvp(sitr:iter));
    else
        avg_dvp = smean(params.hists.dvp(sitr:iter));
        muNEWT.mu_fact = 1.2;
        muNEWT.mu_delta = 5;
    end
    if iter == 10000000000000
        params.sigma = 1;
    end
    if avg_dvp > muNEWT.mu_delta
        params.sigma = params.sigma * (muNEWT.mu_fact);
    else
        params.sigma = params.sigma / (muNEWT.mu_fact);
    end

    params.sigma = min(muNEWT.mu_max, max(muNEWT.mu_min, params.sigma));
end

end

function result = Qmap(X, Q)

if iscell(X)
    result = cell(length(X),1);
    for i =1:length(X)
        result{i} =  Q{i} * X{i};
    end
else
    result = Q * X;
end
end

function result = Amap(X, At)
if iscell(X) && iscell(At)
    result = 0;
    for i =1:length(X)
        result = result + At{i}' * X{i};
    end
elseif iscell(X) && ~iscell(At)
    result = 0;
    for i = 1:length(X)
        result = result + At' * X{i};
    end
elseif ~iscell(X) && iscell(At)
    result = 0;
    for i = 1:length(At)
        result = result + At{i}' * X;
    end
elseif ~iscell(X) && ~iscell(At)
    result = At' * X;
end
end




function result = ATmap(X, Q)

if iscell(X) && iscell(Q)
    result = cell(length(X),1);
    for i =1:length(X)
        result{i,1} =  Q{i} * X{i};
    end
elseif iscell(X) && ~iscell(Q)
    result = cell(length(X),1);
    for i =1:length(X)
        result{i,1} =  Q * X{i};
    end
elseif ~iscell(X) && iscell(Q)
    result = cell(length(Q),1);
    for i =1:length(Q)
        result{i,1} =  Q{i} * X;
    end
else
    result = Q * X;
end
end

function result = idmap(X)
result =  X;
end

function out = Zeromap(X)
if iscell(X)
    out= cell(length(X),1);
    for i =1:length(X)
        out{i} = 0;
    end
else
    out = 0;
end
end

function [AX, AXorg] = AXmap(X, K, At, Lchol)
AXorg = AXfun(K,At,X);
AX = fwsolve(Lchol, AXorg);
end

function Aty = Atymap(y, K, At, Lchol)

Aty = Atyfun(K, At, bwsolve(Lchol, y));
end

function [At,b,C,params] = preprocess_SDP(At,b,C,opts,params)

% K = model.K;
% At = model.At;
% b = model.b;
% C = model.C;


%% scale the data and cholesky decomposition
params.mdim = length(b);
params.nblock = length(params.pblk);
DA = speye(params.mdim);
bscale = 1;
Cscale = 1;

if opts.scale_data
    % normA: norm of each row of A
    normA = zeros(params.mdim,1);
    for k = 1:(params.nblock)
        normA = normA + sum(At{k}.*At{k})';
    end
    normA = max(1,sqrt(normA));
    DA = spdiag(1./normA);
end
for k = 1:params.nblock
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

b = b / bscale; %单位化
C = 1 / Cscale * C;
if params.isL
    params.xL = params.l/bscale;
    params.xU = params.u/bscale;
end
% model_new.K = K;
% model_new.At = At;
% model_new.b = b;
% model_new.C = C;

objscale = bscale*Cscale;
scale.bscale = bscale;
scale.Cscale = Cscale;
scale.objscale = objscale;
scale.scale_data = opts.scale_data;
params.scale = scale;
params.Lchol = Lchol;
end