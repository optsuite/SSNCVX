%% Test_multi_block: test the problem that has multi function block 
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
problemtype = 'QP';
datadir = '../data/QP';
fname{1} = 'AUG2DCQP';
fname{2} = 'CONT-101';
fname{3} = 'CONT-201';
fname{4} = 'CVXQP3_L';
seed = 2024;
rng(seed);
for i = 16%
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
    [m,n] = size(AE);
    A =AE;
    b = A*ones(size(A,2),1);
    C = cell(1);
    C{1} = randn(size(A,2),1);

    % opts setting
    opts.resmin = 0.05;
    opts.m = length(b);
    opts.sigx3l = 0.5;
    opts.sigx3m = 0.5;
    opts.sigx3u = 0.5;
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;


    % pblk setting
    x0 = zeros(m,1);
    At = {A'};
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l2';
    pblk{1}.topk = 1000;
    pblk{1}.size = n;
    pblk{1}.coefficient = 1;

    Bt = speye(n);
    l{1} = 0*ones(size(L));
    u{1} = 20000*ones(size(L));

    lb = b;
    ub = b;
    try
        Q= full(Q.Qmat);
    end
    % f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = n;
    f{1}.coefficient = 0.5;

    % Solve
    [xopt, out] = SSNCVX(x0,pblk,[],f,Q,C,l,u,At,lb,ub,opts);

    %% Two block problem
    opts.fap = 0;
    f2{1,1} = f{1};
    f2{2,1} = f{1};
    pblk2{1,1} = pblk{1};
    pblk2{2,1} = pblk{1};
    At2{1,1} = At{1};
    At2{2,1} = At{1};
    C2{1,1} = C{1};
    C2{2,1} = C{1};
    l2{1,1} = l{1};
    l2{2,1} = l{1};
    u2{1,1} = u{1};
    u2{2,1} = u{1};
    lb = lb*2;
    ub = ub*2;
    Q2{1,1} = Q;
    Q2{2,1} = Q;

    [xopt, out] = SSNCVX([],pblk2,[],f2,[],C2,l2,u2,At2,lb,ub,opts);
end
out.totaltime
