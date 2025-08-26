%% Test_SPCA: test the sparse PCA problem
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
clear;
% clc;
test_dir = '..';
addpath(genpath(test_dir));
dir_data = '../data';

dataset = "SPCA";
probnames = sPCAprobs;
table_str = [];
timegeo = [];
file_len = length(probnames);
for i = [ 1] 
    %% One block problem
    probname = probnames{i};
    model = SDPdata2model(dataset,probname,dir_data);
    n = size(model.C,1);
    Acell{1} = speye(n);
    blk = cell(1,2);
    blk{1,1} = 's';
    blk{1,2} = n;
    K = cell(1);
    K{1} = struct;
    K{1}.type = 's';
    K{1}.size = n;
    At = svec(blk, Acell, 1);
    A = At';
    b = 1;
    nrmC = norm(model.C);
    C = {-model.C/nrmC};

    % opts setting
    opts.adaplambda = 1;

    [m ,n]=size(A);
    % pblk setting   
    pblk{1} = struct;
    pblk{1}.type = 's';
    pblk{1}.size = size(C{1},1);
    pblk{1}.coefficient = 1;

    lb = b;
    ub = b;
    opts.K = K;
    opts.m = length(b);
    % f setting   
    f{1} = struct;
    f{1}.type = 'l1';
    f{1}.size = n;
    f{1}.coefficient = 1/nrmC;
    
    % solve
    [xopt, out] = SSNCVX([],pblk,[],f,[],C,[],[],At,lb,ub,opts);

    %% Two block problem
    pblk2{1,1} = pblk{1};
    pblk2{2,1} = pblk{1};
    At2{1,1} = At{1};
    At2{2,1} = At{1};
    C2{1,1} = C{1};
    C2{2,1} = C{1};
    lb = lb*2;
    ub = ub*2;
    f2{1,1} = f{1};
    f2{2,1} = f{1};
    opts.K{1,1} = K{1};
    opts.K{2,1} = K{1};
    [xopt, out] = SSNCVX([],pblk2,[],f2,[],C2,[],[],At2,lb,ub,opts);

end



