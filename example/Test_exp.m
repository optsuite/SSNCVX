%% Test_exp: test the Logistic regression problem
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
%% Problem
problemtype = 'QP';
datadir = '../data/QP';
fname{1} = 'AUG2DCQP';
fname{2} = 'CONT-101';
fname{3} = 'CONT-201';
fname{4} = 'CVXQP3_L';
fname{5} = 'LISWET5';

seed = 2024;
rng(seed);
for i = 1
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
    A = AE;
    b = A*ones(size(A,2),1);
    C = cell(1);
    C{1} = randn(size(A,2),1);
    x0 = zeros(m,1);
    % opts setting
    opts = struct();
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    At = {A'};
    [m ,n] = size(A);

    % pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.size = n;
    pblk{1}.coefficient = 1;

    % f setting
    f{1} = struct;
    f{1}.type = 'exp';
    f{1}.size = n;
    f{1}.coefficient = 1;
    % Solve
    [xopt, out] = SSNCVX([],pblk,[],f,[],[],[],[],[],[],[],opts);

    %% Two block problem
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     f2{1,1} = f{1};
     f2{2,1} = f{1};

    [xopt, out] = SSNCVX([],pblk2,[],f2,[],[],[],[],[],[],[],opts);
end
out.totaltime