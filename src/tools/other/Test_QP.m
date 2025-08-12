addpath(genpath('../'));
clear
%% QP
problemtype = 'QP';
datadir = '../data/QP';
fname{1} = 'AUG2DCQP';
fname{2} = 'CONT-101';
fname{3} = 'CONT-201';
fname{4} = 'CVXQP3_L';
fname{5} = 'LISWET5';
fname{6} = 'LISWET11';
fname{7} = 'POWELL20';
fname{8} = 'Q25FV47';
fname{9} = 'QPILOTNO';
fname{10} = 'QSHIP12L';
fname{11} = 'STADAT1';
fname{12} = 'UBH1';
idxMaros = 1:12;
fname{13} = 'lipa50a';
fname{14} = 'sko56';
fname{15} = 'tai50a';
fname{16} = 'wil50';
fname{17} = 'esc64a';
idxQAP = 13:17;
fname{18} = 'be100.1';
fname{19} = 'bqp100-1';
fname{20} = 'gka5b';
fname{21} = 'Portfolio5';
fname{22} = 'Portfolio8';
seed = 2024;
rng(seed);
for i = 14% 1  3 11 12
    %% load data
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    if exist([probname,'.mat'])
        load([probname,'.mat'])
    else
        fprintf('\n Can not find the file in UCIdata');
        fprintf('\n ');
        return
    end
        [m,n] = size(AE);
        A =AE;
        b = A*ones(size(A,2),1);
        C = cell(1);
        C{1} = randn(size(A,2),1);
        Amap  = @(x) A*x;
        ATmap = @(x) A'*x;
        opts.Amap = @(x) Amap(x);
        opts.ATmap = @(x) ATmap(x);
        lambda = 1;
    if exist([probname,'.mat'])
        load([probname,'.mat'])
    else
        fprintf('\n Can not find the file in UCIdata');
        fprintf('\n ');
        return
    end


    %%

    stoptol = 1e-4;%stopping tol
    tol = 1e-8;
    % opts.maxits =  10000
    % opts.maxtime = 10000;
    % opts.solver = 1; 
    % opts.record = 1;
    % opts.gtol = tol;
    % opts.xtol = 1e-9;
    % opts.ftol = 1e-9;
    % opts.solver_init = [];
    % opts.opts_init.record = 0;
    % opts.opts_init.tau   = 1e-10;
    % opts.Amap = @(x) Amap(x);
    % opts.ATmap = @(x) ATmap(x);
    % opts.opts_init.maxit = 100;
    % opts.opts_init.gtol  = opts.gtol*1e0;
    % opts.A = A;
    % opts.AT = A';




    % opts.Ascale = 1;
    % opts.is_CG = 0;
    % opts.resmin = 0.05;
    % opts.tol = 1e-6;

    opts.b = b;
    % opts.lambda = lambda;

    opts.m = m;

    % opts.sigma = 1;
    opts.linratio = 1;
    opts.lfactor = 1;

    % opts.sigyl = 1;
    % opts.sigym = 1;
    % opts.sigyu = 1;
    % opts.sigzl = 1;
    % opts.sigzm = 1;
    % opts.sigzu = 1;
    % opts.sigrl = 1;
    % opts.sigrm = 1;
    % opts.sigru = 1;
    % opts.sigvl = 1;
    % opts.sigvm = 1;
    % opts.sigvu = 1;
    % 
    % opts.sigx1l = 0.5;
    % opts.sigx1m = 0.5;
    % opts.sigx1u = 0.5;
    % opts.sigx2l = 0.5;
    % opts.sigx2m = 0.5;
    % opts.sigx2u = 0.5;
    % opts.sigx3l = 0.5;
    % opts.sigx3m = 0.5;
    % opts.sigx3u = 0.5;
    % opts.sigx4l = 0.5;
    % opts.sigx4m = 0.5;
    % opts.sigx4u = 0.5;
    % opts.cgratio = 0.1;
    % opts.gamma1 = 0.5;
    % opts.gamma2 = 0.9;
    % opts.gamma3 = 5;
    % opts.cgtol = 1e-7;
    opts.muopts.mu_update_itr = 10;
    opts.muopts.mu_fact = 1;
    % opts.resratio = 1e-4;
    % opts.adaplambda = 1;
    opts.basename = 'test';
    x0 = zeros(m,1);

    At = {A'};

    opts.t_adap = 0;
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l2';
    pblk{1}.topk = 1000;
    pblk{1}.size = n;
    pblk{1}.coefficient = 1;
    Bt = speye(n);
    if exist('L','var')
        l{1} = 0*ones(size(L));
        u{1} = 20000*ones(size(L));
    else
        l = [];
        u = [];
    end

    if strcmp(problemtype,'QP')
        lb = b;
        ub = b;
        try
            Q= full(Q.Qmat);
        end
        f{1} = struct;
        f{1}.type = 'square';
        f{1}.size = n;
        f{1}.coefficient = 0.5;
    end
    % opts.cgmin = 50;
    % opts.cgmed = 300;
    % opts.method = 'iterative';
    % opts.cgmax = 300;
    % C{1} = randn(size(C{1}));
    % C{1,1} = [C{1,1},C{1,1}];
    % x0 = [x0,x0];

    % profile off
    % profile on
    opts.tol = 1e-5;
    [xopt out] = SSNCVX(x0,pblk,[],f,Q,C,l,u,At,lb,ub,opts);
    % profile off
    % profile on
    [xopt, out] = SSNCVX(x0,pblk,[],f,[],C,l,u,At,lb,ub,opts);
    0.5*norm(xopt.var{1},2)^2 + norm(xopt.var{1}-5,1) + dot_ssn(xopt.var{1},C)
    % profile off
    % profile viewer
    % [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],At,[],[],opts);
    % profile off
    % profile viewer
    1;
    % profile off
    % profile viewer
end
out.totaltime

% function data = prepare_cache_ARNT(data, A,opts)
% y1 = data.YP; z1 = data.XP;
% data.Bx = -opts.ATmap(y1);
% if A == 1
%     data.BTz = -opts.Amap(z1);
% end
% end
