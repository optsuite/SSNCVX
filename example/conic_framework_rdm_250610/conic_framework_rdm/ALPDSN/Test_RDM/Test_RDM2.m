
file_dir = '.';
root_dir = '..';
addpath(file_dir);
%     addpath(genpath([root_dir '/SDPNAL+v1.0']));
addpath(genpath([root_dir '/backup']));
addpath(genpath(pwd));
data_dir = [root_dir '/sdp_data/'];

    files = rdmproblem;


file_len = length(files);
prob_range = 1:file_len;
%     big_index = [1,3,7,9,13,15,19,21,23,25,27,38,39,40,41];
big_index = [];
prob_range(big_index) = [];
file_len = length(prob_range);
table_str = '';

save_root = strcat(root_dir,'/results/RDM/');
if ~exist(save_root,'dir')
    mkdir(save_root)
end

save_root_mat = strcat(save_root,'mat/');
if ~exist(save_root_mat,'dir')
    mkdir(save_root_mat)
end

save_root_eig = strcat(save_root,'eig/');
if ~exist(save_root_eig,'dir')
    mkdir(save_root_eig)
end

save_root_res = strcat(save_root,'res/');
if ~exist(save_root_res,'dir')
    mkdir(save_root_res)
end

timegeo = [];
for gamma2 = [0.9  ] % 0.85 0.95 0.8
    %     for gfactor = [1.5 2 2.5 3 4]
    %     for gfactor2 = [5 8 10 15]
    %     for mu_update_itr = [15 12 10]%15 12 10 8
    %     for mu_factor =[1.1 1.2 1.3 1.5] %1.1 1.2  1.3 1.5
    for i = 1%1:file_len % 4 16 [49 50 177 178 185 263]
        %         ssigxl = num2str(sigxl); %测试sigmax
        %         ssigxu = num2str(sigxu);
        %         ssigyl = num2str(sigyl); %测试sigmax
        %         ssigyu = num2str(sigyu);
        %         scgmax = num2str(cgmax);%测试cg迭代次数
        %         scgmin = num2str(cgmin);
        %                       sgamma1 = num2str(gamma1);
        sgamma2 = num2str(gamma2);
        %         mu_update_itr = 15;
        %         mu_factor = 1.2;
        %         smu_update_itr = num2str(mu_update_itr); %测试sigmay
        %         smu_factor = num2str(mu_factor);
        %         smu_update_itr = strcat('new',smu_update_itr);
        save_root = strcat(root_dir,'/results/RDM/rerdm2x0.3_0.4_cg10gamma1_0.4');
        %             save_root = [save_root,'_',sgamma1,'_',sgamma2,'/'];
        save_root = [save_root,'_',sgamma2,'/'];
        %         save_root = strcat(save_root,scgmin,scgmax,'/');
        %         save_root = [save_root,scgmin,scgmax,'/'];
        %         save_root = [save_root,'_',ssigyl,'_',ssigyu,'/'];
        %            save_root = [save_root,'_',sgamma1,'/'];
        %         save_root = [save_root,'_',scgmax ,'_',scgmin,'/'];
        if ~exist(save_root,'dir')
            mkdir(save_root)
        end

        save_root_mat = strcat(save_root,'mat/');
        if ~exist(save_root_mat,'dir')
            mkdir(save_root_mat)
        end

        save_root_eig = strcat(save_root,'eig/');
        if ~exist(save_root_eig,'dir')
            mkdir(save_root_eig)
        end

        save_root_res = strcat(save_root,'res/');
        if ~exist(save_root_res,'dir')
            mkdir(save_root_res)
        end


        relpath = files{prob_range(i)};
        [~, basename1, basename2] = fileparts(relpath);
        basename = [basename1, basename2];
        file = [data_dir, 'RDM/', relpath, '.mat'];
        load(file);
        %         blk = blk(19,:);
        %         At = At(19,:);
        %         C = C(19,:);
        [b,At,cnz] = data_process(blk,At,b);
        table_str = [table_str basename];
        L = [];
        U = [];
        OPTIONS.tol = 1e-6;
        OPTIONS.printlevel = 0;
        OPTIONS.stopoption = 0;
        opts = [];
        opts.ADMMopts.imaxit = 0;
        opts.tol = 1e-6;
        opts.record = 1;
        opts.basename = basename;
        opts.cgmin = 50;% best:100
        opts.cgmed = 300;% best:100
        opts.cgmax = 300;% best:700
        opts.sigxl = 0.3; % best:0.3 modify
        opts.sigxm = 0.4; % best:0.4 modify
        opts.sigxu = 0.4; % 0.3
        opts.sigyl = 1; % best:0.2
        opts.sigym = 1;
        opts.sigyu = 1;
        opts.save_root = save_root;
        opts.gfactor = 5;
        opts.gfactor2 = 12;

        opts.gamma2 = 0.9;
        opts.gamma1 = 0.4;% 0.3 0.35 0.4
        opts.gamma3 = 5;
        opts.muopts.mu_fact = 1.5; % 1.3
        opts.muopts.mu_update_itr = 10; % 15 modify
        profile on
        [X,out,y,S] = SSNCP2RDM(blk,At,C,b,opts);%RDM 全部投影 RDM2 部分投影
        profile off
        profsave
%         profile viewer
        %         L = {0};
        %         U = {inf};
        %         [X,out,y,S] = SSNCPplus(blk,At,C,b,L,U,opts);
        %         [obj,X,~,y,S,Z,~,~,out,runhist] = ...
        %             sdpnalplus(blk,At,C,b,L,U,[],[],[],OPTIONS);
        iter = out.iter;
        time = out.totaltime
        pinf = out.pinf;
        dinf = out.dinf;
        relgap = out.gap;
        pobj = out.pobj;
        dobj = out.dobj;
        %         pinf = norm(AXfun_sdpnal(blk,At,X) - b)/(1 + norm(b));
        %         dinf = ops(ops(ops(ops(Atyfun_sdpnal(blk,At,y),'+',S),'+',Z), ...
        %             '-',C),'norm')/(1 + ops(C,'norm'));
        %         pobj = full(ops(ops(C,'.*',X),'sum'));
        %         dobj = dualobj(b, y, Z, L, U);
        %         relgap = abs(pobj - dobj)/(1+abs(pobj)+abs(dobj));

        X_item = X{1};
        [n,~] = size(X_item);
        m = length(y);

        if test_type == 1 || test_type == -1 || n<3000
            save_path = strcat(save_root_mat,basename,'_',mat2str(test_id),'.mat');
            save(save_path,'X','y','S','out');
        end

        X_eig = eig(X_item);
        S_eig = eig(S{1});
        save_path = strcat(save_root_eig,basename,'_eig_',mat2str(test_id),'.mat');
        save(save_path,'X_eig','S_eig');
        timegeo = [timegeo time];
        % X_rank = sum(X_eig-max(X_eig)*rank_thres>0);
        % S_rank = sum(S_eig-max(S_eig)*rank_thres>0);
        X_rank = sum(X_eig-max(X_eig)*relgap>0);
        S_rank = sum(S_eig-max(S_eig)*relgap>0);

        etaK1 = out.K1; etaK2 = 0;
        etaC1 = out.C1; etaC2 = 0;

        table_str = [table_str, sprintf('& %d & %d & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e & %d & %.1f & %d & %d', ...
            m,n, pobj, pinf, dinf, relgap, etaK1, etaK2, etaC1, etaC2, iter, time, X_rank, S_rank)];
        table_str = [table_str ' \\ \hline' newline];

        save_path = strcat(save_root_res,basename,'_res_',mat2str(test_id),'.mat');
        save(save_path,'m','n','pobj','pinf','dinf', 'relgap', 'etaK1', 'etaK2', 'etaC1', 'etaC2', 'iter', 'time', 'X_rank', 'S_rank');


    end
    %     fprintf('时间的几何平均值为:%.2e,增长因子1为%.1e,增长因子2为%.1e\n',geo_mean(timegeo),gfactor,gfactor2)
    %         fprintf('时间的几何平均值为:%.2e,cg迭代次数最小为%d,cg迭代次数最大为%d\n',geo_mean(timegeo),cgmin,cgmax)
    %         fprintf('时间的几何平均值为:%.2e,cg每次迭代次数%.2e,sigma增大指数%.2e\n',geo_mean(timegeo),mu_update_itr,mu_factor)
    fprintf('\nx_0.3_0.4_cg10gamma10.4时间的几何平均值为:%.2e,cg每次迭代次数%.2e\n',geo_mean(timegeo),gamma2)
    timegeo  = [];
end
% end

disp(newline);
disp(table_str);

save_path = strcat(save_root,'test',mat2str(test_type),'_',mat2str(test_id),'.txt');
fid = fopen(save_path,'w+');
fprintf(fid,'%s',table_str);

