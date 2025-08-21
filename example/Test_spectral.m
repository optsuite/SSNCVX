addpath(genpath('../'));
clear
%% Nuclear
ori_image = imread('peppers256.png');
ori_image = double(ori_image)/256;
colormap('gray')
noi_image = ori_image + 0.1*randn(size(ori_image));
imshow(noi_image);
1;
for i = 12% 1  3 11 12
    [n1,n2] = size(noi_image);
 



    

    %% opts setting
    opts.tol = 1e-5;
    opts.sigyl = 0.5;
    opts.sigym = 0.5;
    opts.sigyu = 0.5;

    x0 = zeros(n1,n2,1);

    %% pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'nuclear';
    pblk{1}.size = [n1,n2];
    pblk{1}.coefficient = 3;
    B.out_size = [n1, n2];

    %% f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = [2*n1,n2];
    f{1}.coefficient = 1;
    f{1}.shift = noi_image;


    opts.method = 'iterative';
    [xopt, out] = SSNCVX([],pblk,B,f,[],[],[],[],[],[],[],opts);
    %% Two block problem

     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};


     f2{1,1} = f{1};
     f2{2,1} = f{1};

    [xopt, out] = SSNCVX([],pblk2,B,f2,[],[],[],[],[],[],[],opts);


end
