function prob = convert_RDM_to_mosek(At,blk,C,b)

%% 已知 blk At b C 将该四元素转化为Mosek的问题形式prob
% blk: n*2 cell 第一列表示变量类型：'l'或's'，其中'l'表示非负向量，'s'表示半正定矩阵；
%               第二列表示每个变量对应的维度（若维度为一个数组，则表示该正定矩阵由若干个半正定小块组成，每个小块的维度对应数组中的元素，并且每小块排列在对角线上）
% At: 系数矩阵的转置
% b: 约束的右端项
% C: 目标函数中的系数矩阵和系数向量
%
% Remark: 具体blk,At,b,C的含义可参考SDPT3的文档 (https://bpb-us-w2.wpmucdn.com/blog.nus.edu.sg/dist/b/11036/files/2019/10/guide4-0-draft.pdf)
%         prob中每个变量的含义参考Mosek的文档，包括但不限于第6.7节 (https://docs.mosek.com/latest/toolbox/index.html)
%
%% code
[num_var,~] = size(blk);   % 一共的变量个数
c_cell = [];   % 目标函数中'l'类型变量的系数
bardim = [];   % 's'类型变量的维数
num_PSD = 0;   % 's'类型变量的个数
At_PSD={};   % 's'类型变量在约束中的稀疏系数矩阵

% 's'类型变量在约束中的稀疏系数矩阵转化为行列值三元素
r = {};   % 行
c = {};   % 列
v = {};   % 值

a_cell = [];   % 'l'类型变量在约束中的稀疏系数矩阵

blx = [];   % 'l'类型变量bound下界
bux = [];   % 'l'类型变量bound上界

for i = 1:num_var   % 遍历每个变量
    if blk{i,1} == 'l'   % 若变量类型为'l'
        c_cell = [c_cell;C(i)];   % 存储系数向量
        a_cell = [a_cell;At(i)];   % 存储系数矩阵
        blx = [blx,zeros(1,blk{i,2})];   % 下界
        bux = [bux,inf*ones(1,blk{i,2})];   % 上界
    end
    if blk{i,1} == 's'   % 若变量类型为's'
        if size(blk{i,2},2) == 1   % 若该变量是一个大块
            num_PSD = num_PSD + 1;
            [r1,c1,v1] = find(tril(C{i,1}));   % 存储目标函数中稀疏的系数矩阵的下三角的三元素：行列值
            r{num_PSD,1} = r1;
            c{num_PSD,1} = c1;
            v{num_PSD,1} = v1;
            At_PSD{num_PSD,1} = At{i};   % 存储约束中的系数矩阵
        else
            temp_blk_arr = blk{i,2};   % 若该变量是由若干个小块组成
            temp_At = At{i};
            len_ = size(temp_blk_arr,2);   % 多少个小块
            tt_C = 0;
            tt_At = 0;
            for ii = 1:len_
                start_C = tt_C + 1;
                stop_C = tt_C + temp_blk_arr(ii);
                temp_C = C{i,1}(start_C: stop_C, start_C: stop_C);   % 提取出每个小块的C中的对应的矩阵
                [r1,c1,v1] = find(tril(temp_C));   % 将下三角转化为对应的三元素
                r{num_PSD+ii,1} = r1;
                c{num_PSD+ii,1} = c1;
                v{num_PSD+ii,1} = v1;
                len_At = temp_blk_arr(ii)*(temp_blk_arr(ii)+1)/2;   % 每个小块At中的系数矩阵的长度
                start_At = tt_At + 1;
                stop_At = tt_At + len_At;
                At_PSD{num_PSD+ii,1} = temp_At(start_At:stop_At,:);
                tt_C = tt_C + temp_blk_arr(ii);
                tt_At = tt_At + len_At;
            end
            num_PSD = num_PSD + len_;
        end
        bardim = [bardim,blk{i,2}];
    end
end

% The scalar part, as in linear optimization examples
prob.c = full(cell2mat(c_cell)');    % 目标函数中'l'类型变量的系数 从cell转化为向量
if size(b,1)==1
    b = b';
end
prob.blc = b;   % 约束的右端项下界
prob.buc = b;   % 约束的右端项上界
prob.a = cell2mat(a_cell)';   % 'l'类型变量在约束中的稀疏系数矩阵
if isempty(prob.a)
    prob.a = zeros(size(b,1), 0);
end
prob.blx = blx;
prob.bux = bux;

% Dimensions of PSD variables
prob.bardim = bardim;   % 所有's'类型变量的维数

% Coefficients in the objective
barc_subj = [];   % 's'类型的变量
barc_subk = [];   % 行
barc_subl = [];   % 列
barc_val= [];     % 值 上述四个变量维度要一致
for i = 1:num_PSD   % 遍历's'类型的所有变量 第i个变量
    barc_subj = [barc_subj,repmat(i,1,length(v{i,1}))];
    barc_subk = [barc_subk, r{i,1}'];
    barc_subl = [barc_subl, c{i,1}'];
    barc_val = [barc_val, v{i,1}'];
end
prob.barc.subj = barc_subj;
prob.barc.subk = barc_subk;
prob.barc.subl = barc_subl;
prob.barc.val = barc_val;

% Coefficients in the constraint
num_constraint = size(b,1);   % 约束的个数
bara_subi = cell(1, num_constraint);   % 约束
bara_subj = cell(1, num_constraint);   % 's'类型的变量
bara_subk = cell(1, num_constraint);   % 行
bara_subl = cell(1, num_constraint);   % 列
bara_val = cell(1, num_constraint);    % 值
arr=zeros(num_PSD,1);
for i = 1:num_PSD   % 遍历变量 第i个变量
    tmp2 = sqrt(2)*triu(ones(bardim(i)),1) + speye(bardim(i),bardim(i));   % 原始数据中的At与Mosek的稀疏矩阵有所不一样
    tmp2 = tmp2(:);
    dd = tmp2(find(tmp2));
    temp_subk = [];
    temp_subl = [];

    temp_bara_subi = [];
    temp_bara_subj = [];
    temp_bara_val = [];
    for j = 1:num_constraint      % 遍历约束
        temp_vaule = full(At_PSD{i,1}(:,j));
        temp_vaule = temp_vaule./dd;
        temp_vaule = temp_vaule';
        I = find(temp_vaule);   % 找出不为0元素的位置
        count = length(I);
        temp_bara_subi = [temp_bara_subi, j*ones(1, count)];
        temp_bara_subj = [temp_bara_subj, i*ones(1, count)];
        temp_bara_val = [temp_bara_val, temp_vaule(I)];
        if isempty(I)   % 判断是否有非0元
            continue;
        end
        t_subk = [];
        t_subl = [];
        for ii = 1:bardim(i)
            t_subk = [t_subk, repmat(ii,1,ii)];
            t_subl = [t_subl, 1:ii];
        end
        temp_subk = [temp_subk,t_subk(I)];
        temp_subl = [temp_subl,t_subl(I)];
    end
    arr(i)=length(temp_bara_subi);
    bara_subi{i} = temp_bara_subi;
    bara_subj{i} = temp_bara_subj;
    bara_val{i} = temp_bara_val;
    bara_subk{i} = temp_subk;
    bara_subl{i} = temp_subl;
end

prob.bara.subi = horzcat(bara_subi{1:num_constraint});
prob.bara.subj = horzcat(bara_subj{1:num_constraint});
prob.bara.subk = horzcat(bara_subk{1:num_constraint});
prob.bara.subl = horzcat(bara_subl{1:num_constraint});
prob.bara.val = horzcat(bara_val{1:num_constraint});
prob.n=sum(prob.bardim);
prob.m=num_constraint;
end


