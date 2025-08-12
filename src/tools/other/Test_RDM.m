clear;
% clc;
root_dir = '../..';
test_dir = '..';
addpath([root_dir]);
addpath(genpath(test_dir));
dir_results = "../results";
dir_data = addpath_data();
dir_logs = 'logs';

if ~exist("dataset", "var")
    dataset = "RDM";
end
instances = "all";

if  numel(instances) == 1 && instances == "all"

    if dataset == "theta"
        probnames = thetaprobs;
    elseif dataset == "qap"
        probnames = qapprobs;
    elseif dataset == "DIMACS"
        probnames = ["nb", "nb_L1", "nb_L2", "nb_L2_bessel", "nql180", "nql30", "nql60", "qssp180", "qssp30", "qssp60", "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled", "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
    elseif dataset == "CBLIB"
        probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph",   "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", "beam30", "beam7", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", "db-plate-yield-line"];
    elseif dataset == "sdplib"
        probnames = ["arch0", "arch2", "arch4", "arch8", "control1", "control10", "control11", "control2", "control3", "control4", "control5", "control6", "control7", "control8", "control9", "equalG11", "equalG51", "gpp100", "gpp124-1", "gpp124-2", "gpp124-3", "gpp124-4", "gpp250-1", "gpp250-2", "gpp250-3", "gpp250-4", "gpp500-1", "gpp500-2", "gpp500-3", "gpp500-4", "hinf1", "hinf10", "hinf11", "hinf12", "hinf13", "hinf14", "hinf15", "hinf2", "hinf3", "hinf4", "hinf5", "hinf6", "hinf7", "hinf8", "hinf9", "infd1", "infd2", "infp1", "infp2", "maxG11", "maxG32", "maxG51", "maxG55", "maxG60", "mcp100", "mcp124-1", "mcp124-2", "mcp124-3", "mcp124-4", "mcp250-1", "mcp250-2", "mcp250-3", "mcp250-4", "mcp500-1", "mcp500-2", "mcp500-3", "mcp500-4", "qap10", "qap5", "qap6", "qap7", "qap8", "qap9", "qpG11", "qpG51", "ss30", "theta1", "theta2", "theta3", "theta4", "theta5", "theta6", "thetaG11", "thetaG51", "truss1", "truss2", "truss3", "truss4", "truss5", "truss6", "truss7", "truss8"];
    elseif dataset == "theta+"
        probnames = thetaprobs;
    elseif dataset == "biq"
        probnames = biqprobs;
    elseif dataset == "fap"
        probnames = fapprobs;
    elseif dataset == "RDM"
        probnames = rdmproblem;
    elseif dataset == "R1TA"
        probnames = R1TAprobs;
    elseif dataset == "rcp"
        probnames = rcpprobs;
    end

else
    probnames = instances;
end

save_root = strcat(dir_results,'/' ,dataset);

table_str = [];
timegeo = [];
file_len = length(probnames);
for i = 1
    probname = probnames{i};
    model = SDPdata2model(dataset,probname,dir_data);
       kk = length(model.K);
         modeltmp = model;
            lenk = 0;
            model.K = {};
            for ii = 1:kk
                tmpC2 = modeltmp.C{ii};
                Attmp = modeltmp.At{ii};
                cone2 = modeltmp.K{ii};
                blktmp = modeltmp.blk{ii,2};
                cumsum_size = [0 cumsum(cone2.size)];
                sqcumsum_size = [0 cumsum((cone2.size.^2 + cone2.size)/2  )];
                for jj = 1 : length(modeltmp.K{ii}.size)
                    cone2.size = modeltmp.K{ii}.size(jj);
                    if strcmp(modeltmp.K{ii}.type,'s')
                        model.C{lenk+jj,1} = tmpC2(cumsum_size(jj)+1:cumsum_size(jj+1),cumsum_size(jj)+1:cumsum_size(jj+1));
                        model.At{lenk+jj,1} = Attmp(sqcumsum_size(jj)+1:sqcumsum_size(jj+1),:);
                        model.blk{lenk+jj,1} = 's';
                        model.blk{lenk+jj,2} = blktmp(jj);
                        model.K{lenk+jj,1} = cone2;
                    else
                        model.C{lenk+jj,1} = tmpC2;
                        model.At{lenk+jj,1} = Attmp;
                        model.K{lenk+jj,1} = cone2;
                        model.blk{lenk+jj,1} = 'l';
                        model.blk{lenk+jj,2} = blktmp(jj);
                    end
                end
                lenk = lenk + length(modeltmp.K{ii}.size);
            end
    A = model.At{1};
    b = model.b;
    C = model.C;
    Amap  = @(x) A*x;
    ATmap = @(x) A'*x;
    opts.Amap = @(x) Amap(x);
    opts.ATmap = @(x) ATmap(x);
    lambda = 1;

    stoptol = 1e-4;%stopping tol
    tol = 1e-6;
 
    opts.maxits =  10000;
    opts.maxtime = 10000;
    opts.solver = 1; 
    opts.record = 1;
    opts.gtol = tol;
    opts.xtol = 1e-9;
    opts.ftol = 1e-9;
    opts.solver_init = [];
    opts.opts_init.record = 0;
    opts.opts_init.tau   = 1e-10;%1e-3
    K = model.K;
    At = model.At;
    opts.Amap = @(X) AXmap(X, K, At, Lchol);
    opts.ATmap = @(y) Atymap(y, K, At, Lchol);
    opts.opts_init.maxit = 100;
    opts.opts_init.gtol  = opts.gtol*1e0;
    opts.A = A;
    opts.AT = A';
    
    opts.maxit = 100000;
    opts.sub_maxit = 20;
    opts.adap_subiter = 1;
    opts.t_adap = 1;
    
    opts.Ascale = 1;
    opts.is_CG = 0;
    opts.resmin = 0.05;
    opts.tol = 1e-7;

    opts.b = b;

    opts.stoptol = 1e-2;
    opts.tau1 = 0.8; opts.tau2 = 0.8;opts.resmin = 500;opts.sigma = 30;
    opts.linratio = 0.5;
    opts.lfactor = 0.98;


    
    opts.sigyl = 1;
    opts.sigym = 1;
    opts.sigyu = 1;
    opts.sigzl = 1;
    opts.sigzm = 1;
    opts.sigzu = 1;
    opts.sigrl = 1;
    opts.sigrm = 1;
    opts.sigru = 1;
    opts.sigvl = 1;
    opts.sigvm = 1;
    opts.sigvu = 1;

    opts.sigx1l = 0.5;
    opts.sigx1m = 0.5;
    opts.sigx1u = 0.5;
    opts.sigx2l = 0.5;
    opts.sigx2m = 0.5;
    opts.sigx2u = 0.5;
    opts.sigx3l = 0.5;
    opts.sigx3m = 0.5;
    opts.sigx3u = 0.5;
    opts.sigx4l = 0.3;
    opts.sigx4m = 0.4;
    opts.sigx4u = 0.4;
    opts.cgratio = 0.001;
    opts.gamma1 = 1;
    opts.gamma2 = 1;
    opts.gamma3 = 1;
    opts.cgtol = 1e-7;
    opts.muopts.mu_update_itr = 15;
    opts.muopts.mu_fact = 1.1;
    opts.resratio = 1e-4;
    opts.basename = 'test';


    opts.t_adap = 0;
    [m ,n]=size(A);
    for i = 1:length(model.blk)
    pblk{i,1} = struct;
    pblk{i,1}.type = model.blk{i,1};
    pblk{i,1}.size = model.blk{i,2};
    pblk{i,1}.coefficient = 1;
    end
    % l = L;
    if exist('L','var')
    l = -2*ones(size(L));
    u = 20*ones(size(L));
    else
    l = [];
    u = [];
    end
    problemtype = 'SDP';
    lb = b;
    ub = b;
    opts.cgmin = 50; 
    opts.cgmed = 700;
    opts.cgmax = 700; 
    opts.K = K;
    opts.m = length(b);
    % profile off
    % profile on
    % [xopt out] = SSNCVX(x0,pblk,[],f,Q,C,l,u,At,lb,ub,opts);
    % profile off
        % profile on
     % [xopt, out] = SSNCVX(x0,pblk,[],f,[],C,l,u,At,lb,ub,opts);
     % profile off
      % profile viewer
     [xopt, out] = SSNCVX([],pblk,[],[],[],C,[],[],At,lb,ub,opts);
     % profile off
     % profile viewer
     1;
    % profile off
    % profile viewer
end
out.totaltime



% function data = prepare_cache_ARNT(data, A,opts)
% y1 = data.YP; z1 = data.XP;
% data.Bx = -opts.ATmap(y1);
% if A == 1
%     data.BTz = -opts.Amap(z1);
% end
% end
