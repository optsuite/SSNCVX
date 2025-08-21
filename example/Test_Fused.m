addpath(genpath('../'));
clear
%% Fused Lasso problem
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
for i = 11
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



    %% opts setting
    opts.lambda = 20;
    opts.resmin = 500;
    opts.sigma = 5;
    opts.linratio = 0.5;
    opts.lfactor = 0.98;
    % opts = deflalt_opts(crho,i);
    opts.sigma = 1/opts.sigma;
    opts.sigzl = 1;
    opts.sigx4l = 1;
    opts.sigzm = opts.sigzl;
    opts.sigzu = opts.sigzl;
    opts.sigx4m = opts.sigx4l;
    opts.sigx4u = opts.sigx4l;
    opts.method = 'direct';


    x0 = zeros(m,1);

    %% pblk setting
    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'fused';
    pblk{1}.size = n;
    pblk{1}.coefficient = lambda;
    pblk{1}.coefficient2 = lambda;
    [Bmap,BTmap] = FLBmap(n);
    pblk{1}.Binput = struct();
    pblk{1}.Binput.Bmap = Bmap;
    pblk{1}.Binput.BTmap = BTmap;

    %% f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = n;
    f{1}.coefficient = 0.5;
    f{1}.shift = b;

    
    %% Solve
     [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);

end
out.totaltime
