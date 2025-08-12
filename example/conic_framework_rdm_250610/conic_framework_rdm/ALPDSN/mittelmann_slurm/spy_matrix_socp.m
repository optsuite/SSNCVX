clear; clc;
root_dir = '../..';
test_dir = '..';
result_path = '/bicmr/home/optsuite/wt/projects/conic_framework_rdm/ALPDSN/mittelmann_slurm/result';
addpath([root_dir]);
startup(root_dir);
addpath("../../sdpt3fun/Mexfun/");
dir_data = addpath_data();
dataset = 'pass';

probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2a", ...
            "2013_firL2L1alph", "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firLinf",...
            "2013_wbNRL", "beam7", "beam30", "chainsing-50000-1", "chainsing-50000-2",...
            "chainsing-50000-3", "db-joint-soerensen", "db-plate-yield-line"];

for i = 1:18
    probname = probnames{i};
    fid_path = [result_path, '/', int2str(i), '/info.txt'];
    fid = fopen(fid_path, 'w');
    [model, model_original] = data2model(dataset,probname,dir_data);
    fprintf(fid, 'Detecting problem %s from mittelmann SOCP dataset ...\n', probname);
    fprintf(fid, '\nmodel_original:\n');
    print_std_model_info(fid, model_original);
    fprintf(fid, '\nmodel:\n');
    print_std_model_info(fid, model);
    fprintf(fid, '\n');

    At_original = model_original.At;
    K = model_original.K;
    Lchol.Rperm = speye(size(At_original{1}, 2));
    trans.sp_info = detect_sp_info(At_original, K, Lchol,2 , 0.1, fid);

    % Plot and save the first spy plot
    figure;
    spy(cell2mat(At_original));
    saveas(gcf, [result_path, '/', int2str(i), '/At.jpg']);
    close;

    % Plot and save the second spy plot
    figure;
    spy(trans.sp_info.At_row_sp_merged);
    saveas(gcf, [result_path, '/', int2str(i), '/At_sp.jpg']);
    close;

    % Plot and save the third spy plot
    figure;
    spy(trans.sp_info.AHAt_mat22_At);
    saveas(gcf, [result_path, '/', int2str(i), '/At_mat22.jpg']);
    close;

    % Plot and save the fourth spy plot
    figure;
    spy(trans.sp_info.lhs_structure);
    saveas(gcf, [result_path, '/', int2str(i), '/lhs.jpg']);
    close;
end
