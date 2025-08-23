%% Test_l2con: test the problem that has l2 constraint
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
fname{2} = 'log1p.E2006.train';
fname{3} = 'E2006.test';
fname{4} = 'log1p.E2006.test';
fname{5} = 'pyrim_scale_expanded5';
fname{6} = 'triazines_scale_expanded4';
fname{7} = 'abalone_scale_expanded7';
fname{8} = 'bodyfat_scale_expanded7';
fname{9} = 'housing_scale_expanded7';
fname{10} = 'mpg_scale_expanded7';
fname{11} = 'space_ga_scale_expanded9';
fname{12} = 'E2006.train';
for i =1
    %% One block problem
    probname = [datadir,filesep,fname{i}];
    fprintf('\n Problem name: %s \n', fname{i});
    load([probname,'.mat'])
    [m,n] = size(A);
    Bt = A';

    % opts setting
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
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
    f{1}.type = 'l2con';
    f{1}.size = n;
    f{1}.coefficient = 1;
    seed = 2024;
    rng(seed)
    b = 10*ones(n,1);
    f{1}.shift = b;

    % solve
    [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);
    %% Two block problem
    
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

