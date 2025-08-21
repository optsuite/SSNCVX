addpath(genpath('../'));
clear
problemtype = 'QP';
datadir = '../data/QP/QP_Mat/QPDATA_Mat';
fname = {    'AUG2D','AUG2DC','AUG2DCQP','AUG2DQP','BOYD1','BOYD2','CONT-100','CONT-101', ...
    'CONT-200','CONT-201','CONT-300','DTOC3','LISWET1','LISWET10','LISWET11','LISWET12',...
    'LISWET2','LISWET3','LISWET4','LISWET5','LISWET6','LISWET7','LISWET8','LISWET9',...
    'UBH1'};
seed = 2025;
rng(seed);
total = [];
totalobj = [];
totalpinf = [];
totaldinf = [];
totalK1 = [];
for i = 1:length(fname)% 1  3 11 12
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    load('QSTAIR.mat')
    % if exist([probname,'.mat'])
    %     load([probname,'.mat'])
    % else
    %     fprintf('\n Can not find the file in UCIdata');
    %     fprintf('\n ');
    %     return
    % end
        A = model.A;
        n = length(model.lb);
        A = ones(n,1)';
        b = 1;
        m = length(b);
        Q = model.Q;
  

        C{1} = model.c;

    % k = 10;
    % n = 2048; p = ceil(n*0.1); % k = 1:10
    % F = sprandn(n, p, 0.1); D = sparse(diag(sqrt(p)*rand(n,1)));
    % % Q = cov(F') + D;
    % Q = full(cov(F') + D);
    % A = ones(n,1)';
    % b = 1;
    % m = 1;
    % C = cell(1);
    % C{1} = randn(n,1); gamma = 1.0;
    % lambda = 1;



    %%

    tol = 1e-5;
    opts.maxits =  10000;
    opts.maxtime = 10000;
    opts.solver = 1; 
    opts.record = 1;


    opts.Amap = @(x) Amap(x);
    opts.ATmap = @(x) ATmap(x);
    opts.A = A;
    opts.AT = A';

    opts.maxit = 100000;
    opts.resmin = 0.05;
    opts.tol = 1e-6;

    opts.m = m;


    opts.linratio = 0.5;
    opts.lfactor = 0.98;


    opts.sigx1l = 0.5;
    opts.sigx1m = 0.5;
    opts.sigx1u = 0.5;
    opts.sigx2l = 0.5;
    opts.sigx2m = 0.5;
    opts.sigx2u = 0.5;
    opts.sigx3l = 0.5;
    opts.sigx3m = 0.5;
    opts.sigx3u = 0.5;
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    opts.cgratio = 0.1;
    opts.gamma1 = 0.5;
    opts.gamma2 = 0.9;
    opts.gamma3 = 5;
    opts.cgtol = 1e-7;
    opts.muopts.mu_update_itr = 10;
    opts.muopts.mu_fact = 1;
    opts.resratio = 1e-4;
    opts.adaplambda = 1;
    opts.basename = 'test';
    x0 = zeros(m,1);

    At = {A'};


    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'box';
    pblk{1}.topk = 1000;
    pblk{1}.size = n;
    pblk{1}.l = 0*ones(n,1);
    pblk{1}.u = inf*ones(n,1);
    pblk{1}.coefficient = 1;
    pblk{1}.coefficient2 = 2;
    Bt = speye(n);


    lb = b;
    ub = b;
    
    l = {model.lb-1000};
    u = {model.ub+1000};
    opts.cgmin = 50;
    opts.cgmed = 300;
    opts.method = 'iterative';
    opts.cgmax = 300;

    [xopt, out] = SSNCVX(x0,pblk,[],[],Q,C,l,u,At,lb,ub,opts);

    total = [total out.totaltime];
    totalobj = [totalobj out.pobj];
    totalpinf = [totalpinf out.pinf];
    totaldinf = [totaldinf out.dinf];
    totalK1 = [totalK1 out.K1];
    0.5*norm(xopt.var{1},2)^2 + norm(xopt.var{1}-5,1) + dot_ssn(xopt.var{1},C)

end



