%% Test_fap: test the fap problem
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
clear;
test_dir = '..';
addpath(genpath(test_dir));

dir_data = addpath_data();
dataset = "fap";
probnames = fapprobs;

file_len = length(probnames);
for i = 1
    %% One block problem
    probname = probnames{i};
    model = SDPdata2model(dataset,probname,dir_data);
    A = model.At{1};
    b = model.b;
    C = model.C;
    L = model.L;
    U = model.U;
    lambda = 1;
    K = model.K;
    At = model.At;

    % opts setting
    opts.sigx4l = 0.4;
    opts.sigx4m = 0.4;
    opts.sigx4u = 0.5;
    opts.cgratio = 0.01;
    opts.lambda = 5;
    opts.muopts.mu_update_itr = 12;
    opts.muopts.mu_fact = 1.5;
    [m ,n] = size(A);

    % pblk setting
    pblk{1} = struct;
    pblk{1}.type = 's';
    pblk{1}.size = size(C{1},1);
    pblk{1}.coefficient = 1;
    opts.K = K;
    opts.m = length(b);
    opts.fap = 1;
    l = L;
    u = U;
    lb = b;
    ub = b;

    % solve
    [xopt, out] = SSNCVX([],pblk,[],[],[],C,l,u,At,lb,ub,opts);

     %% Two block problem
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     At2{1,1} = At{1};
     At2{2,1} = At{1};
     C2{1,1} = C{1};
     C2{2,1} = C{1};
     lb = lb*2;
     ub = ub*2;
     l2{1,1} = l{1};
     l2{2,1} = l{1};
     u2{1,1} = u{1};
     u2{2,1} = u{1};
     opts.K{1,1} = K{1};
     opts.K{2,1} = K{1};
     [xopt, out] = SSNCVX([],pblk2,[],[],[],C2,l2,u2,At2,lb,ub,opts);

end
out.totaltime

