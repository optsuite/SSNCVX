%% Test_l1con: test the problem that has l1 constraint
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
%% Lasso
problemtype = 'Lasso';
datadir = '../data/Lasso';
fname{1} = 'uci_CT';
for i =1
    %% One block problem
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    load([probname,'.mat'])

    [m,n] = size(A);
    Bt = A';
    % opts setting
    opts.sigx2l = 0.5;
    opts.sigx2m = 0.5;
    opts.sigx2u = 0.5;
    opts.resratio = 1e-4;
    x0 = zeros(m,1);
    At = A';
    
    % pblk setting
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.topk = 5;
    pblk{1}.size = n;
    Bt = eye(n);
    pblk{1}.coefficient = 1;

    % f setting
    f{1} = struct;
    f{1}.type = 'l1con';
    f{1}.size = n;
    f{1}.coefficient = 1;
    seed = 2024;
    rng(seed)
    b = 10*ones(n,1);
    f{1}.shift = b;

    % Solve
    [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);

     %% Two block porblem
     x02{1,1} = x0;
     x02{2,1} = x0;
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     Bt2{1,1} = Bt;
     Bt2{2,1} = Bt;
     f2{1,1} = f{1};
     f2{2,1} = f{1};
     [xopt, out] = SSNCVX(x02,pblk2,Bt2,f2,[],[],[],[],[],[],[],opts);
     norm(Bt'*xopt.var{1} + b,2)
end
out.totaltime

