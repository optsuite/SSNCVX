clear;
clc;
root_dir = '../..';
test_dir = '..';
addpath([root_dir]);
addpath(genpath(test_dir));
dir_results = "../results";

dir_data = "../data/"
dir_logs = 'logs';

dataset = "DIMACS";
probnames = ["nb", "nb_L1", "nb_L2", "nb_L2_bessel", "nql180", "nql30", "nql60", "qssp180", "qssp30", "qssp60", "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled", "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
save_root = strcat(dir_results,'/' ,dataset);


table_str = [];
timegeo = [];
file_len = length(probnames);
for i = 3
    %% load data
    probname = probnames{i};
    model = SOCPdata2model(dataset,probname,dir_data);
  
    A = model.At{1};
    n = size(model.C{1},1);
    b = model.b;
    C = model.C;
    Amap  = @(x) A*x;
    ATmap = @(x) A'*x;
    K = model.K;
    At = model.At;

    


    
    %% opts setting
    opts.method = 'direct';
    opts.sigx4l = 0.4;
    opts.sigx4m = 0.4;
    opts.sigx4u = 0.4;
    opts.resratio = 1e-4;


    opts.t_adap = 0;
    [m ,n]=size(A);

    %% pblk setting
    for i = 1:length(K)
        pblk{i,1} = struct;
        pblk{i,1}.type = K{i}.type;
        pblk{i,1}.size = K{i}.size;
        pblk{i,1}.coefficient = 1;
    end

    lb = b;
    ub = b;
    opts.K = K;
    opts.m = length(b);

     [xopt, out] = SSNCVX([],pblk,[],[],[],C,[],[],At,lb,ub,opts);

end
out.totaltime



