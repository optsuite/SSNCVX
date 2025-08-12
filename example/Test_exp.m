addpath(genpath('../'));
clear
%% QP
problemtype = 'QP';
datadir = '../data/QP';
fname{1} = 'AUG2DCQP';
fname{2} = 'CONT-101';
fname{3} = 'CONT-201';
fname{4} = 'CVXQP3_L';
fname{5} = 'LISWET5';
fname{6} = 'LISWET11';
fname{7} = 'POWELL20';
fname{8} = 'Q25FV47';
fname{9} = 'QPILOTNO';
fname{10} = 'QSHIP12L';
fname{11} = 'STADAT1';
fname{12} = 'UBH1';
idxMaros = 1:12;
fname{13} = 'lipa50a';
fname{14} = 'sko56';
fname{15} = 'tai50a';
fname{16} = 'wil50';
fname{17} = 'esc64a';
idxQAP = 13:17;
fname{18} = 'be100.1';
fname{19} = 'bqp100-1';
fname{20} = 'gka5b';
fname{21} = 'Portfolio5';
fname{22} = 'Portfolio8';
seed = 2024;
rng(seed);
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



    [m,n] = size(AE);
    A =AE;
    b = A*ones(size(A,2),1);
    C = cell(1);
    C{1} = randn(size(A,2),1);
    x0 = zeros(m,1);
    opts = struct();
    opts.sigx4l = 0.5;
    opts.sigx4m = 0.5;
    opts.sigx4u = 0.5;
    At = {A'};

    [m ,n]=size(A);
    pblk{1} = struct;
    pblk{1}.type = 'l1';
    pblk{1}.size = n;
    pblk{1}.coefficient = 1;
    pblk{1}.coefficient2 = 2;



    f{1} = struct;
    f{1}.type = 'exp';
    f{1}.size = n;
    f{1}.coefficient = 1;

    [xopt, out] = SSNCVX([],pblk,[],f,[],[],[],[],[],[],[],opts);

     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     f2{1,1} = f{1};
     f2{2,1} = f{1};

    [xopt, out] = SSNCVX([],pblk2,[],f2,[],[],[],[],[],[],[],opts);


end
out.totaltime