%% Test_linftycon: test the problem that has linfty norm constraint
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
datadir = '../data/Lasso';
fname{1} = 'uci_CT';
for i =1% 1  3 11 12
    %% One block problem
    probname = [datadir,filesep,fname{i}];
    load([probname,'.mat'])
    [m,n] = size(A);

    % opts setting
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    x0 = zeros(m,1);
    At = A';
    [m ,n]=size(A);

    % pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.size = n;
    Bt = eye(n);
    pblk{1}.coefficient = 1;

    % f setting
    f{1} = struct;
    f{1}.type = 'linftycon';
    f{1}.size = n;
    f{1}.coefficient = 1;
    b = 10*ones(n,1);
    f{1}.shift = b;
    
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
     norm(Bt'*xopt.var{1} + b,'inf')
end
out.totaltime

