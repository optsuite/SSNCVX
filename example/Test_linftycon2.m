addpath(genpath('../'));
clear
datadir = '../data/Lasso';
fname{1} = 'uci_CT';
seed = 2024;
rng(seed)
for i = 1
    probname = [datadir,filesep,fname{i}];
    load([probname,'.mat'])
    [m,n] = size(A);


    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    x0 = zeros(m,1);
    At = A';

    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'linftycon';
    pblk{1}.topk = 5;
    pblk{1}.size = n;
    Bt = eye(n);
    pblk{1}.coefficient = 1;

    b = 10*rand(n,1);
    pblk{1}.shift = b;
    f{1} = struct;
    f{1}.type = 'l1';
    f{1}.size = n;
    f{1}.coefficient = 1;



    [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);
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
