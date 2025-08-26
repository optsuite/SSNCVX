%% Test_QP: test the QP problem
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
problemtype = 'QP';
datadir = '../data/QP/';
fname = {'AUG2D','AUG2DC','AUG2DCQP','AUG2DQP'};
seed = 2025;
rng(seed);
for i = 1:length(fname)% 1  3 11 12
    %% One block problem
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    if exist([probname,'.mat'])
        load([probname,'.mat'])
    else
        fprintf('\n Can not find the file in UCIdata');
        fprintf('\n ');
        return
    end
        A = model.A;
        n = length(model.lb);
        A = ones(n,1)';
        b = 1;
        m = 1;
        Q = model.Q;
  
        C = cell(1);
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

    % opts setting
    opts.solver = 1; 
    opts.record = 1;
    opts.resmin = 0.05;
    opts.m = m;

    opts.linratio = 0.5;
    opts.lfactor = 0.98;
    opts.cgtol = 1e-7;
    opts.resratio = 1e-4;
    opts.adaplambda = 1;
    opts.basename = 'test';
    x0 = zeros(m,1);

    At = {A'};

    % pblk setting
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

    % opts setting
    opts.cgmin = 50;
    opts.cgmed = 300;
    opts.method = 'iterative';
    opts.cgmax = 300;

    [xopt, out] = SSNCVX(x0,pblk,[],[],Q,C,[],[],At,lb,ub,opts);

    0.5*norm(xopt.var{1},2)^2 + norm(xopt.var{1}-5,1) + dot_ssn(xopt.var{1},C)

end



