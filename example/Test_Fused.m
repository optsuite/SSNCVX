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
for i = 12
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
    opts.sigma = 30;
    opts.linratio = 0.5;
    opts.lfactor = 0.98;
    opts = deflalt_opts(crho,i);
    opts.sigma = 1/opts.sigma;
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
    pblk{1}.topk = 10;
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



    

     [xopt, out] = SSNCVX(x0,pblk,Bt,f,[],[],[],[],[],[],[],opts);

end
out.totaltime

function opts = deflalt_opts(crho,i)

    if crho == 1e-3
        if i == 1 opts.sigzl = 0.63; opts.sigx4l = 0.6; opts.sigma = 1;opts.nu = 5;opts.resmin = 7;opts.t_adap = 0; end %1e-3
        if i == 2 opts.sigzl = 0.6; opts.sigx4l = 0.6;opts.resmin = 100;opts.sigma = 30; opts.lfactor = 0.8;opts.linratio = 0.4;end %？
        if i == 3 opts.sigzl = 1; opts.sigx4l = 1; opts.sigma = 1;opts.nu = 5;opts.resmin = 0.05;opts.t_adap = 0; end%1e-3
        if i == 4 opts.sigzl = 1; opts.sigx4l = 0.8;opts.resmin = 500;opts.sigma = 70; opts.lfactor = 0.8;opts.linratio = 0.6;end %？
        if i == 5 opts.sigzl = 1; opts.sigx4l = 1;opts.resmin = 0.1; opts.t_adap = 0;opts.sigma = 1;end % ?
        if i == 6 opts.sigzl = 1; opts.sigx4l = 1;opts.resmin = 500;opts.sigma = 10;end
        if i == 7 opts.sigzl = 1; opts.sigx4l = 1;opts.resmin = 2;opts.sigma = 10; opts.lfactor = 0.8;end
        if i == 8 opts.tau1 = 0.8; opts.tau2 = 0.8;opts.resmin = 300;opts.sigma = 1; opts.lfactor = 0.9;end
        if i == 9 opts.sigzl = 0.7; opts.sigx4l = 0.7;opts.resmin = 300;opts.sigma = 5;opts.lfactor = 0.9; opts.linratio = 0.8;end%？
        if i == 11 opts.sigzl = 1; opts.sigx4l = 1; opts.sigma = 1;opts.nu = 5;opts.resmin = 0.05;opts.t_adap = 0; end%1e-3
        if i == 12 opts.sigzl = 1; opts.sigx4l = 1; opts.sigma = 5;opts.nu = 5;opts.resmin = 0.05;opts.t_adap = 0; end%1e-3
        if i == 10 opts.sigzl = 0.8; opts.sigx4l = 0.8;opts.resmin = 0.02;opts.t_adap = 0; opts.sigma = 0.5;end
    elseif crho == 1e-4
        if i == 1 opts.sigzl = 0.5; opts.sigx4l = 0.5; opts.sigma = 50;opts.nu = 5;opts.resmin = 5;opts.t_adap = 0;opts.lfactor = 0.8;opts.adaplambda = 1; end %1e-3
        if i == 2 opts.sigzl = 0.6; opts.sigx4l = 0.6;opts.resmin = 700;opts.sigma = 600; opts.lfactor = 0.98;opts.linratio = 0.8;end %？
        if i == 3 opts.sigzl = 1; opts.sigx4l = 1; opts.sigma = 1;opts.nu = 5;opts.resmin = 0.05;opts.t_adap = 0; end%1e-3
        if i == 4 opts.sigzl = 1.2; opts.sigx4l = 0.8;opts.resmin = 700;opts.sigma = 500; opts.lfactor = 0.98;opts.linratio = 0.8;end %？
        if i == 5 opts.sigzl = 0.6; opts.sigx4l = 0.6;opts.resmin = 0.1; opts.t_adap = 0;opts.sigma = 0.01;opts.adaplambda = 1;opts.lfactor = 0.7; end % ?
        if i == 6 opts.sigzl = 1; opts.sigx4l = 1;opts.resmin = 500;opts.sigma = 1;opts.adaplambda = 1;   end
        if i == 7 opts.sigzl = 1; opts.sigx4l = 1;opts.resmin = 2;opts.sigma = 10; opts.lfactor = 0.8;end
        if i == 8 opts.sigzl = 0.8; opts.sigx4l = 0.8;opts.resmin = 300;opts.sigma = 1;opts.adaplambda = 1; end
        if i == 9 opts.sigzl = 0.6; opts.sigx4l = 0.6;opts.resmin = 300;opts.sigma = 0.1;opts.lfactor = 0.9; opts.linratio = 0.8;opts.adaplambda = 1;end%？
        if i == 11 opts.sigzl = 1; opts.sigx4l = 1; opts.sigma = 0.1;opts.nu = 5;opts.resmin = 5;opts.t_adap = 0; end%1e-3
        if i == 12 opts.sigzl = 0.6; opts.sigx4l = 0.6; opts.sigma = 5;opts.nu = 5;opts.resmin = 0.5;opts.t_adap = 0;opts.adaplambda = 1; end%1e-3
        if i == 10 opts.sigzl = 0.8; opts.sigx4l = 0.8;opts.resmin = 0.02;opts.t_adap = 0; opts.sigma = 0.5;end
    end
end