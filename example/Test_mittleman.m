clear; clc;
root_dir = '../..';
test_dir = '..';
addpath([root_dir]);
startup(root_dir);
addpath("../../sdpt3fun/Mexfun/");

dir_results = "../results";
dir_logs = 'logs';

dataset = "pass";
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
    elseif dataset == "rdm"
        probnames = rdmproblem;
    elseif dataset == "R1TA"
        probnames = R1TAprobs;
    elseif dataset == "rcp"
        probnames = rcpprobs;
    elseif dataset == "sedumi"
        probnames = sedumi_problem;
    elseif dataset == "gen" || dataset == "pass"
        probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2a", ...
        "2013_firL2L1alph", "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firLinf",...
         "2013_wbNRL", "beam7", "beam30", "chainsing-50000-1", "chainsing-50000-2",...
          "chainsing-50000-3", "db-joint-soerensen", "db-plate-yield-line"];
    end
else
    probnames = instances;
end

save_root = strcat(dir_results,'/' ,dataset);
dir_data = '../data/';
table_str = [];
timegeo = [];
file_len = length(probnames);

for i = [6]
    probname = probnames{i};
    [model, model_original] = data2model(dataset,probname,dir_data);
    
    for j = 1:length(model.At)
        At{j,1} = model.At{j};
        pblk{j,1} = struct;
        pblk{j,1}.type = model.K{j}.type;
        pblk{j,1}.size = model.K{j}.size;
        pblk{j,1}.coefficient = 1;
        C{j,1} = model.C{j};
    end
    b = model.b;
    opts.m = length(b);




    kk = length(model.K);
    lenk = 0;
    modeltmp = model;

    table_str = [table_str char(probname)];
    L = [];
    U = [];

    opts.log_path = "";
    opts.barrier = 0;
    opts.mu = 1e-1;
    opts.tol = 1e-6;
    opts.record = 1;
    opts.save_path = "../../record_weit/archive";
    opts.basename = probname;
    opts.cgmin = 1000;
    opts.cgmed = 500; 
    opts.cgmax = 500;
    opts.sigxl = 0.1;
    opts.sigxm = 0.1;
    opts.sigxu = 0.1;
    opts.xfactor = 1;
    opts.project_option = 0;
    opts.use_chol = 0;
    opts.use_AAtchol = 0;
    opts.scale_A = 1;
    opts.scale_bc_flag = 1;
    opts.gamma1 = 0.5;
    opts.gamma2 = 0.7;
    opts.gamma3 = 5;
    opts.cgtol = 1e-6;
    opts.sigyl = 1;
    opts.sigym = 1; 
    opts.sigyu = 1;
    opts.sigzl = 0.5;
    opts.sigzm = 0.5;
    opts.sigzu = 0.8;
    opts.sigql = 0.1;
    opts.sigqm = 0.1;
    opts.sigqu = 0.3;
    opts.cgratio = 0.05;
    opts.gfactor = 5;
    opts.gfactor2 = 10;
    opts.resratio = 1e-4;

    opts.random_init_x = 1e-3; % 初始点的缩放系数，这里表示乘1e-3
    opts.gap_pinf_elimination = 0; % 是否使用使用投影手段把gap和pinf降到0，只需调dinf和knorm，这个功能matlab测试后发现迭代步数反而增加，先不使用

    opts.mintau1 = 0; % tau1最小值，小于此值时，被设置为此值
    opts.mintau2 = 0; % tau2最小值，小于此值时，被设置为此值

    opts.method = 'direct'; % SOCP都使用直接法
    opts.swapDirectStep = 0; % 在第几步时强制切换为直接法，设为0
    opts.system_opt = 2;
    opts.direct_solver = 'ldl';
    opts.socp_formula = 1;
    opts.sys2_sparse_strategy = 0;

    opts.linesearch = 0; % 是否使用线搜索，0或1
    opts.stepsize = 1; % 线搜索初始步长
    opts.linesearch_minstep = 0.6; % 线搜索最小步长
    opts.linesearch_rho = 0.9; % 线搜索步长缩放因子
    opts.linesearch_threshold_dinf = -1; % dinf大于此值时不开启线搜索，-1表示不使用阈值
    opts.ls_const = 1;

    opts.muopts.init = 1; % sigma初始值
    opts.muopts.mu_update_itr = 10; % sigma每12步更新一次
    opts.muopts.mu_fact = 1.5; % sigma每次更新、乘上或除以1.5。
    opts.muopts.NEWT.mu_min = 1e-6; % sigma不会小于1e-6
    opts.muopts.NEWT.mu_max = 1e6; % sigma不会大于1e6
    opts.muopts.NEWT.strategy = 'monotony'; % 使用单调策略，即每次都乘1.5，而非自动选择乘上或除以1.5

    opts.smooth = 1; % 是否使用光滑化，0或1
    opts.smooth_mu = 1e-12;
    opts.smooth_linesearch_update_mu = 0; % 是否使用线搜索更新mu，0或1。
    opts.smooth_linesearch_mu_fact = 0.6;
    opts.smooth_ratio_update_mu = 0; % 是否使用固定比例更新mu，0表示不使用，若使用则为一个0~1间的数，这里表示每次迭代时mu会乘上0.6. 
    opts.smooth_gap_update_mu = 0; % 是否使用gap更新mu，0或1
    opts.smooth_gap_update_mu_pow = 2; % mu_next = gaporg^trans.smooth_gap_update_mu_pow * trans.smooth_gap_update_mu_coeff;
    opts.smooth_gap_update_mu_coeff = 1; % 同上
    opts.smooth_threshold = 1e-17;

    opts.adaptive_mu = 0;
    opts.param_nm = 1;
    opts.param_am = 1;
    opts.max_non_monotone = 0;
    opts.iter_am = 10;
    opts.get_param = 1;

    if opts.get_param
        opts = get_param(opts, i);
    end

    lb = b;
    ub = b;

    if opts.record >=1
        fprintf("*************************************************************************************\n");
        fprintf("Solving problem %s from %s dataset ...\n", probname, dataset);
    end

    % [X,out,y,S,Zb] = SSNCPplusall_mittlemanmit(model,model_original,opts);
    [xopt, out] = SSNCVX([],pblk,[],[],[],C,[],[],At,lb,ub,opts);

    iter = out.iter;
    time = out.totaltime
    fprintf("iter: %d, time: %.2e\n", iter, time);
    % append to file

end

function [opts] = get_param(opts, i)
    if i == 1 % 2013_dsNRL
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.01;
        opts.linesearch = 1;
        opts.max_non_monotone = 20;
        opts.muopts.mu_fact = 0.5;
        opts.muopts.mu_update_itr = 1000;
        opts.param_am = 0.8;
        opts.resratio = 0.1;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-06;
        opts.smooth_ratio_update_mu = 1;
        opts.smooth_threshold = 1.00e-13;
        opts.socp_formula = 3;
    elseif i == 2 % 2013_firL1
        opts.gamma1 = 0.6913901;
        opts.gamma2 = 0.340769006;
        opts.linesearch = 1;
        opts.muopts.mu_fact = 22;
        opts.muopts.init = 100;
        opts.muopts.mu_update_itr = 12;
        opts.resratio = 2.593002745;
        opts.scale_A = 1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
    elseif i == 3 % 2013_firL1Linfalph
        opts.scale_A = 1;
        opts.scale_bc_flag = 1;
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.gamma3 = 1.5;
        opts.iter_am = 12;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.1;
        opts.max_non_monotone = 8;
        opts.muopts.mu_fact = 0.7;
        opts.muopts.init = 1;
        opts.muopts.mu_update_itr = 10000;
        opts.param_am = 100;
        opts.resratio = 1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-06;
        opts.smooth_threshold = 1.00e-15;
        opts.socp_formula = 3;
    elseif i == 4 % 2013_firL1Linfeps
        opts.scale_A = 0;
        opts.scale_bc_flag = 1;
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.gamma3 = 1;
        opts.linesearch = 1;
        opts.linesearch_min_step = 1.00e-01;
        opts.max_non_monotone = 4;
        opts.muopts.mu_fact = 0.6;
        opts.muopts.init = 100;
        opts.muopts.mu_update_itr = 10000;
        opts.muopts.mu_fact = 0.5;
        opts.param_am = 100;
        opts.resratio = 1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-05;
        opts.smooth_threshold = 1.00e-15;
        opts.socp_formula = 3;
    elseif i == 5 % 2013_firL2a
        opts.linesearch = 1;
        opts.scale_A = 0;
    elseif i == 6 % 2013_firL2L1alph
        opts.muopts.init = 0.01;
        opts.resratio = 1;
        opts.scale_A = 0;
        opts.smooth_mu = 1e-15;
    elseif i == 7 % 2013_firL2L1eps
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.8;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.639444217;
        opts.ls_const = 1.932571717;
        opts.max_non_monotone = 5;
        opts.muopts.mu_fact = 0.5;
        opts.muopts.init = 1.00e+04;
        opts.muopts.NEWT.mu_max = 1.00e+09;
        opts.muopts.mu_update_itr = 100;
        opts.resratio = 1;
        opts.scale_A = 1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
        opts.smooth = 1;
        opts.smooth_mu = 1.00e-10;
        opts.socp_formula = 3;
        opts.stepsize = 1;
    elseif i == 8 % 2013_firL2Linfalph
        opts.gamma1 = 0.4;
        opts.gamma2 = 0.75;
        opts.linesearch = 1;
        opts.muopts.mu_fact = 10;
        opts.muopts.init = 100;
        opts.muopts.NEWT.mu_max = 1.00e+09;
        opts.muopts.mu_update_itr = 12;
        opts.resratio = 0.1;
        opts.scale_A = 0;
        opts.scale_bc_flag = 1;
        opts.sigxl = 0.5;
        opts.sigxu = 0.5;
        opts.sigyl = 1;
        opts.sigyu = 1;
        opts.stepsize = 1;
    elseif i == 9 % 2013_firL2Linfeps
        opts.linesearch = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1e-05;
        opts.smooth_threshold = 1e-15;
        opts.resratio = 1;
        opts.smooth_ratio_update_mu = 1;
        opts.linesearch_min_step = 0.1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
        opts.sigma = 10;
        opts.scale_A = 0;
        opts.scale_bc_flag = 0;
    
        opts.muopts.mu_update_itr = 1000;
        opts.gamma1 = 0.3;
        opts.gamma3 = 1;
    
        opts.adaptive_mu = 1;
        opts.smooth_linesearch_mu_fact = 0.46;
        opts.max_non_monotone = 50;
        opts.param_am = 100;
        opts.refinement_tol = 8;
        opts.refine_max_iter = 3;
        % opts.gamma1 = 0.3;
        % opts.gamma2 = 0.85;
        % opts.linesearch = 1;
        % opts.muopts.mu_fact = 10;
        % opts.muopts.init = 100;
        % opts.muopts.NEWT.mu_max = 1.00e+09;
        % opts.muopts.mu_update_itr = 50;
        % opts.resratio = 0.001;
        % opts.scale_A = 1;
        % opts.sigxl = 0.3;
        % opts.sigxu = 0.4;
        % opts.smooth = 1;
        % opts.socp_formula = 3;
        % opts.stepsize = 1;
    elseif i == 10 % 2013_firLinf
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.333856951;
        opts.gamma2 = 0.153438045;
        opts.iter_am = 14;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.1;
        opts.max_non_monotone = 5;
        opts.muopts.mu_fact = 0.5;
        opts.muopts.init = 405;
        opts.muopts.mu_update_itr = 160;
        opts.param_am = 100;
        opts.resratio = 0.001;
        opts.scale_A = 1;
        opts.sigxl = 0.197601116;
        opts.sigxu = 0.130925183;
        opts.sigyl = 0.889558081;
        opts.sigyu = 0.83881669;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-06;
        opts.smooth_threshold = 0;
        opts.socp_formula = 3;
        opts.stepsize = 1;
    elseif i == 11 % 2013_wbNRL
        opts.gamma1 = 0.01;
        opts.resratio = 1;
        opts.scale_A = 1;
        opts.smooth = 1;
        opts.smooth_mu = 1.00e-05;
        opts.smooth_ratio_update_mu = 0.5;
        opts.socp_formula = 3;
    elseif i == 12 % beam7
        opts.scale_A = 0;
        opts.scale_bc_flag = 1;
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.linesearch = 1;
        opts.linesearch_min_step = 1.00e-06;
        opts.max_non_monotone = 5;
        opts.muopts.mu_fact = 0.3;
        opts.muopts.init = 0.001;
        opts.muopts.mu_update_itr = 10000;
        opts.param_am = 10;
        opts.resratio = 1;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 4.00e-06;
        opts.smooth_threshold = 1.00e-16;
        opts.socp_formula = 3;
    elseif i == 13 % beam30
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.606795764;
        opts.gamma2 = 0.68804605;
        opts.linesearch = 1;
        opts.linesearch_min_step = 1.00e-06;
        opts.max_non_monotone = 5;
        opts.muopts.mu_fact = 0.2;
        opts.muopts.init = 0.001;
        opts.muopts.mu_update_itr = 555559;
        opts.param_am = 0.1;
        opts.resratio = 0.562805684;
        opts.scale_A = 1;
        opts.sigxl = 0.771996973;
        opts.sigxu = 0.607934249;
        opts.sigyl = 0.355832806;
        opts.sigyu = 0.662545248;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 0.0003;
        opts.smooth_threshold = 5.00e-15;
        opts.socp_formula = 3;
    elseif i == 14 % chainsing-50000-1
        opts.gamma1 = 0.518102261;
        opts.gamma2 = 0.252996919;
        opts.linesearch = 1;
        opts.muopts.mu_fact = 4.978932538;
        opts.resratio = 0.1;
        opts.scale_A = 1;
        opts.sigxl = 0.253799078;
        opts.sigxu = 0.494959889;
        opts.sigyl = 0.9;
        opts.sigyu = 0.543022212;
        opts.smooth = 1;
        opts.socp_formula = 3;
    elseif i == 15 % chainsing-50000-2
        opts.gamma1 = 0.4;
        opts.gamma2 = 0.75;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.3;
        opts.linesearch_rho = 0.8;
        opts.muopts.mu_fact = 10;
        opts.muopts.init = 10;
        opts.muopts.NEWT.mu_max = 400;
        opts.muopts.NEWT.mu_min = 1;
        opts.muopts.mu_update_itr = 12;
        opts.scale_A = 0;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 1;
        opts.sigyu = 1;
        opts.smooth = 1;
        opts.smooth_gap_update_mu_coeff = 1.00e-07;
        opts.smooth_gap_update_mu_pow = 1.5;
        opts.smooth_mu = 1.00e-05;
        opts.smooth_ratio_update_mu = 0.450591289;
        opts.smooth_threshold = 0;
        opts.socp_formula = 3;
        opts.stepsize = 0.9;
    elseif i == 16 % chainsing-50000-3
        opts.gamma1 = 0.146474688;
        opts.gamma2 = 0.350508534;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.6;
        opts.muopts.mu_fact = 2.087021963;
        opts.muopts.NEWT.mu_max = 1000;
        opts.muopts.NEWT.mu_min = 1;
        opts.resratio = 0.726995821;
        opts.scale_A = 0;
        opts.scale_bc_flag = 0;
        opts.sigxl = 0.1;
        opts.sigxu = 0.362283251;
        opts.sigyl = 0.473310698;
        opts.sigyu = 0.154057395;
        opts.smooth_gap_update_mu = 1;
        opts.smooth_gap_update_mu_coeff = 1.00e-05;
        opts.smooth_gap_update_mu_pow = 1.5;
        opts.smooth_mu = 1.00e-05;
        opts.smooth_threshold = 0;
    elseif i == 17 % db-joint-soerensen
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.linesearch = 1;
        opts.linesearch_min_step = 0.1;
        opts.max_non_monotone = 5;
        opts.muopts.mu_fact = 0.5;
        opts.muopts.init = 3;
        opts.muopts.NEWT.mu_max = 10000;
        opts.muopts.NEWT.mu_min = 0.0001;
        opts.muopts.mu_update_itr = 1000;
        opts.param_am = 1000;
        opts.param_nm = 5;
        opts.resratio = 0.5;
        opts.scale_A = 1;
        opts.sigxl = 0.1;
        opts.sigxu = 0.1;
        opts.sigyl = 0.5;
        opts.sigyu = 0.5;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-06;
        opts.smooth_threshold = 4.00e-16;
        opts.socp_formula = 3;
    elseif i == 18 % db-plate-yield-line
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.linesearch = 1;
        opts.linesearch_min_step = 1.00e-06;
        opts.max_non_monotone = 15;
        opts.muopts.mu_fact = 0.6;
        opts.muopts.init = 0.01;
        opts.muopts.mu_update_itr = 1000;
        opts.param_am = 1000;
        opts.resratio = 0.5;
        opts.scale_A = 1;
        opts.sigxl = 0.1;
        opts.sigxu = 0.1;
        opts.sigyl = 0.5;
        opts.sigyu = 0.5;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 0.001;
        opts.smooth_threshold = 1.00e-14;
        opts.socp_formula = 3;
        opts.sys2_sparse_strategy = 3;
    end
end