clear;
root_dir = '../..';
test_dir = '..';
addpath([root_dir]);
addpath(genpath(test_dir));
dir_results = "../results";
dir_data = addpath_data();
dir_logs = 'logs';
dataset = "theta";
instances = "all";
probnames = thetaprobs;


save_root = strcat(dir_results,'/' ,dataset);

table_str = [];
timegeo = [];
file_len = length(probnames);
for i = 10
    probname = probnames{i};
    model = SDPdata2model(dataset,probname,dir_data);
    A = model.At{1};
    b = model.b;
    C = model.C;
    Amap  = @(x) A*x;
    ATmap = @(x) A'*x;
 
    opts.gtol = 1e-6;
    K = model.K;
    At = model.At;


    

    opts.sigx4l = 0.3;
    opts.sigx4m = 0.4;
    opts.sigx4u = 0.4;
    opts.cgratio = 0.1;
    opts.muopts.mu_fact = 1.5;
    opts.adaplambda = 1;


    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 's';
    pblk{1}.topk = 10;
    pblk{1}.size = size(C{1},1);
    pblk{1}.coefficient = 1;


    lb = b;
    ub = b;
    opts.cgmin = 50; 
    opts.cgmed = 300;
    opts.cgmax = 300; 
    opts.K = K;
    opts.m = length(b);

     [xopt, out] = SSNCVX([],pblk,[],[],[],C,[],[],At,lb,ub,opts);

     1;
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     At2{1,1} = At{1};
     At2{2,1} = At{1};
     C2{1,1} = C{1};
     C2{2,1} = C{1};
     lb = lb*2;
     ub = ub*2;
     opts.K{1,1} = K{1};
     opts.K{2,1} = K{1};
     [xopt, out] = SSNCVX([],pblk2,[],[],[],C2,[],[],At2,lb,ub,opts);

end
out.totaltime



% function data = prepare_cache_ARNT(data, A,opts)
% y1 = data.YP; z1 = data.XP;
% data.Bx = -opts.ATmap(y1);
% if A == 1
%     data.BTz = -opts.Amap(z1);
% end
% end
