addpath(genpath('../'));
clear
problemtype = 'Lasso';
datadir = '../data/Lasso';
fname{1} = 'uci_CT';
for i =1
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    load([probname,'.mat'])

    [m,n] = size(A);
    Bt = A';




    

    %% opts setting
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    opts.resratio = 1e-4;

    x0 = zeros(m,1);
    At = A';
    [m ,n]=size(A);

    %% pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.topk = 5;
    pblk{1}.size = n;
    Bt = eye(n);
    pblk{1}.coefficient = 1;

    %% f setting
    f{1} = struct;
    f{1}.type = 'l2con2';
    f{1}.user_defined = 1;
    f{1}.size = n;
    f{1}.coefficient = 1;
    f{1}.pobj = @(X) 0;
    f{1}.dobj = @(X,pblk) dobj(X,pblk);
    f{1}.prox = @(X,coefficient,sigma) prox(X,coefficient,sigma);
    f{1}.Matvec = @(X,par) Matvec(X,par);
    f{1}.DPhi = @(iHW,d) DPhi(iHW,d);
    f{1}.DPhi2 = @(D2,tx2,taux2) DPhi2(D2,tx2,taux2);
    f{1}.DPhi3 = @(D2,D2y,dz,Fx2,taux2) DPhi3(D2,D2y,dz,Fx2,taux2);
    seed = 2024;
    rng(seed)
    b = 10*ones(n,1);
    f{1}.shift = b;



    opts.cgmin = 50;
    opts.cgmed = 300;
    opts.method = 'iterative';
    opts.cgmax = 300;
    %% solve
     [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);
     x02{1,1} = x0;
     x02{2,1} = x0;
     1;
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     Bt2{1,1} = Bt;
     Bt2{2,1} = Bt;
     f2{1,1} = f{1};
     f2{2,1} = f{1};


end
out.totaltime


function [value, genD] = prox(X,coefficient,sigma) 
            nrmx = norm(X);
            lambda = coefficient*sigma;
        if nrmx < lambda
            value = X;
            genD.nrmx = 1;
            genD.type = 1;
            genD.coefficient = coefficient;
            genD.sigma = sigma;
            genD.coe1 = 1;
            genD.coe2 = 0;
        else
            value = lambda*X/nrmx;
            genD.coefficient = coefficient;
            genD.sigma = sigma;
            genD.nrmx = nrmx;
            genD.coe1 = coefficient*sigma/nrmx;
            genD.coe2 = -coefficient*sigma/nrmx;
            genD.X = X/nrmx;
            genD.type = 2;
        end
end

function [value] = dobj(X,pblk)
if isfield(pblk,'shift')
    value =  dot_ssn(pblk.shift,X) + pblk.coefficient*norm(X,2);
else
    value = pblk.coefficient*norm(X,2);
end
end

function out = DPhi(iHW,d)
if iHW.type == 1
    if iscell(d)
        out= d{1};
    else
        out = d;
    end
else
    if iscell(d)
        out = iHW.coe1*d{1} + iHW.coe2*iHW.X'*d{1}*iHW.X;
    else
        out = iHW.coe1*d + iHW.coe2*iHW.X'*d*iHW.X;
    end
end
end

function  tmp = Matvec(Aty,par,p)

        if iscell(Aty)
            Aty = Aty{1};
        end
        if par.D4.type == 1
            tmp = zeros(size(Aty));
        else
            tmp = par.D11{p,1} *Aty + par.D12{p,1} *par.D4.X'*Aty*par.D4.X;
        end
end

function tmp2 = DPhi2(D2,tx2,taux2)
if isfield( D2,'type') &&  D2.type == 2
    tmpcoe = taux2 + (D2.nrmx - D2.sigma*D2.coefficient)/D2.sigma*D2.nrmx;
    tmpcoe = tmpcoe^(-1);
    tmpcoe2 = D2.sigma*D2.coefficient/D2.nrmx;
    D2.coe1 = tmpcoe2*tmpcoe;
    D2.coe2 = -tmpcoe2*tmpcoe;
    tmp2 = D2.coe1*tx2 + D2.coe2*D2.X'*tx2*D2.X;
else
    tmp2 = tx2;
end
end

function dx2 = DPhi3(D2,D2y,dz,Fx2,taux2)
iHWx = D2;
if iHWx.type == 1
    dx2 = dz;
else
    dx2 = iHWx.coe1*dz + iHWx.coe2*iHWx.X'*dz*iHWx.X;
end
dx2 = (dx2 - Fx2);
tmpcoe = taux2 + (D2y.nrmx - D2y.sigma*D2y.coefficient)/D2y.sigma*D2y.nrmx;
tmpcoe = tmpcoe^(-1);
tmpcoe2 = D2y.coefficient/D2y.nrmx;
iHWx.coe1 = tmpcoe;
iHWx.coe2 = tmpcoe*tmpcoe2/(1 - tmpcoe*tmpcoe2);
if iHWx.type == 1
    dx2 = dz;
else
    dx2 = iHWx.coe1*dx2 + iHWx.coe2*iHWx.X'*dx2*iHWx.X;
end
end