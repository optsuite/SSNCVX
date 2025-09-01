%% Test_l1l2: test the problem that has l1l2 norm
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
problemtype = 'Lasso';
datadir = '../data/Lasso';
fname{1} = 'uci_CT';
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
    [m,n] = size(A);
    Bt = A';
    crho = 1e-4;
    lambdamax=norm(Bt*b,'inf');
    lambda=crho*lambdamax;
    stoptol = 1e-4;
    tol = 1e-6;
    % opts setting
    opts.maxits =  10000;
    opts.maxtime = 10000;
    opts.Amap = @(x) Amap(x);
    opts.ATmap = @(x) ATmap(x);
    x0 = zeros(m,1);
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    opts.resmin = 1e-4;

    x0 = zeros(m,1);
    [m ,n]=size(A);
    % pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'l1l2';
    pblk{1}.size = [n,2];
    pblk{1}.coefficient = lambda;

    % f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = n;
    f{1}.coefficient = 0.5;
    f{1}.shift = [-b,-b];
    
     % solve
     [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);
     x02{1,1} = x0;
     x02{2,1} = x0;
     
     %% Two block problem
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     Bt2{1,1} = Bt;
     Bt2{2,1} = Bt;
     f2{1,1} = f{1};
     f2{2,1} = f{1};
     [xopt, out] = SSNCVX(x02,pblk2,Bt2,f2,[],[],[],[],[],[],[],opts);

end
out.totaltime

