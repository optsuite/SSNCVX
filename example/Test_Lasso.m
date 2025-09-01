%% Test_Lasso: test the lasso problem
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
fname{2} = 'E2006.test';
for i = 1% 1  3 11 12
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
    crho = 1e-3;%
    lambdamax=norm(Bt*b,'inf');
    lambda=crho*lambdamax;

    % opts setting
    opts.lambda = 10;
    opts.sigma = 1;
    opts.sigzl = 0.8;
    opts.sigx4l = 0.6;
    opts.nu = 5;
    opts.sigzm = opts.sigzl;
    opts.sigzu = opts.sigzl;
    opts.sigx4m = opts.sigx4l;
    opts.sigx4u = opts.sigx4l;
    x0 = zeros(m,1);
    At = A';
    % pblk setting
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.topk = 10;
    pblk{1}.size = n;
    opts.resmin = 0.07;
    pblk{1}.coefficient = lambda;

    % f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = n;
    f{1}.coefficient = 0.5;
    f{1}.shift = -b;
    opts.method = 'direct';

    % Solve
    [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);


end
out.totaltime
