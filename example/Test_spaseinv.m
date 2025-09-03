%% Test_sparseinv: test the sparse covariance selection problem
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
    A =AE;
    b = A*ones(size(A,2),1);
    C = cell(1);
    C{1} = randn(size(A,2),1);

    % opts setting
    opts.resmin = 0.05;
    opts.cgratio = 0.1;
    opts.cgtol = 1e-7;
    opts.resratio = 1e-4;

    x0 = zeros(100,100);
    seed = 2025;
    rng(seed);
    n = 100;
    A=randn(n,n);
    mean_1=mean(A,1);
    A=A-repmat(mean_1,n,1);
    norm_2=sqrt( sum(A.^2,1) );
    A=A./repmat(norm_2,n,1);
    C = {A*A'};

    % pblk setting
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.topk = 1000;
    pblk{1}.size = n;
    pblk{1}.coefficient = 1;
    lb = b;
    ub = b;
    try
        Q= full(Q.Qmat);
    end

    % f setting
    f{1} = struct;
    f{1}.type = 'logdet';
    f{1}.size = n;
    f{1}.coefficient = 1;
    f{1}.obj = @(X) log(det(X));
    f{1}.dual_obj = @(X) n - det(X);
    f{1}.D2 = @(X) - X;
    
    [xopt, out] = SSNCVX(x0,pblk,[],f,[],C,[],[],[],[],[],opts);

    %% Two block problem
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     C2{1,1} = C{1};
     C2{2,1} = C{1};
     f2{1,1} = f{1};
     f2{2,1} = f{1};


     [xopt, out] = SSNCVX([],pblk2,[],f2,[],C2,[],[],[],[],[],opts);

    -0.5*log(det(xopt.var{1})) + norm(reshape(xopt.var{1},[],1),1) + xopt.var{1}(:)'*C{1}(:);

end
out.totaltime

