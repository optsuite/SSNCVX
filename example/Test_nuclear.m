%% Test_nuclear: test the problem that has nuclear norm
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
addpath(genpath('../'));
clear
for i = 1:2
    %% Load image
    filename = sprintf('image%d.png', i);
    ori_image = imread(filename);
    
    %% Convert to grayscale (handling 1/2/3/4 channels)
    if size(ori_image, 3) == 1
        % Already grayscale (512x512)
        gray_image = ori_image;
    elseif size(ori_image, 3) == 2
        % Grayscale + Alpha (take only the first channel)
        gray_image = ori_image(:,:,1);
    elseif size(ori_image, 3) == 3
        % RGB -> convert to grayscale
        gray_image = rgb2gray(ori_image);
    elseif size(ori_image, 3) == 4
        % RGBA -> discard Alpha, then RGB to grayscale
        gray_image = rgb2gray(ori_image(:,:,1:3));
    else
        error('Unsupported image format: %d channels', size(ori_image, 3));
    end
    
    %% Normalize to [0,1]
    gray_image = double(gray_image) / 255;
    
    seed = 2024;
    rng(seed);
    noi_image = gray_image + 0.1 * randn(size(gray_image));
    total_elements = size(noi_image,1) * size(noi_image,2);
    ones_count = round(total_elements * 0.5);  
    vec = [ones(1, ones_count), zeros(1, total_elements - ones_count)];
    shuffled_vec = vec(randperm(total_elements));
    sampling_matrix = reshape(shuffled_vec, [size(noi_image,1), size(noi_image,2)]);
    [n1,n2] = size(noi_image);
    opts.tol = 1e-3;
    
    if exist('b','var')
    opts.b = b;
    end
    x0 = zeros(n1,n2,1);

    % opts setting
    opts.sigyl = 0.6;
    opts.sigym = 0.6;
    opts.sigyu = 0.8;
    opts.sigx4l = 0.51;
    opts.sigx4m = 0.51;
    opts.sigx4u = 0.51;
    opts.adaplambda = 1;

    % pblk setting
    pblk{1} = struct;
    pblk{1}.type = 'nuclear';
    pblk{1}.size = [n1,n2];
    pblk{1}.coefficient = 10;

    B.Bmap = @(u) sampling_matrix.*u;
    B.BTmap = @(y) sampling_matrix.*y;
    B.out_size = [n1, n2];
    % f setting
    f{1} = struct;
    f{1}.type = 'square';
    f{1}.size = [2*n1,n2];
    f{1}.coefficient = 1;
    f{1}.shift = noi_image;

    % solve
     [xopt, out] = SSNCVX(x0,pblk,B,f,[],[],[],[],[],[],[],opts);

     %% Two block problem
     pblk2{1,1} = pblk{1};
     pblk2{2,1} = pblk{1};
     f2{1,1} = f{1};
     f2{2,1} = f{1};

    [xopt, out] = SSNCVX([],pblk2,B,f2,[],[],[],[],[],[],[],opts);


end

