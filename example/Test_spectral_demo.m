%% Tutorial for Low-Rank Matrix Completion (LRMC)
% This script demonstrates how to solve a matrix completion problem using SSNCVX.
% The problem is defined as:
%
% min ||X||_*
% s.t. P_Omega(X) = P_Omega(M)
%
% where ||.||_* is the nuclear norm (sum of singular values),
% M is the original matrix, and Omega is the set of observed entries.
% The solver addresses the formulation:
% min 0.5 * || P_Omega(X) - P_Omega(M) ||_F^2 + lambda * ||X||_*

clear;clc;
addpath(genpath('../'));

%% Problem Data Construction
% Create a low-rank matrix
n1 = 10;
n2 = 10;
rank = 2;
rng(42);
M = randn(n1, rank) * randn(rank, n2);

% Create a mask for observed entries (50% observed)
p = 0.5;
Omega = rand(n1, n2) < p;
M_observed = M .* Omega;

lambda = 0.1; % Regularization parameter

%% pblk setting for the nuclear norm penalty
% p(X) = lambda * ||X||_*
pblk{1} = struct;
pblk{1}.type = 'nuclear';
pblk{1}.size = [n1, n2];
pblk{1}.coefficient = lambda;

%% Bmap and BTmap for the sampling operator
B.Bmap = @(X) X .* Omega;
B.BTmap = @(Y) Y .* Omega;
B.out_size = [n1, n2];

%% f setting for the data fitting term
% f(y) = 0.5 * ||y - M_observed||_F^2 where y = Bmap(X)
f{1} = struct;
f{1}.type = 'square';
f{1}.size = [n1, n2];
f{1}.coefficient = 0.5;
f{1}.shift = M_observed;

%% Initial guess and solver options
x0 = zeros(n1, n2);
opts = struct(); % Use default options

%% Call the solver
fprintf('Solving a small Matrix Completion problem...\n');
[xopt, out] = SSNCVX(x0, pblk, B, f, [], [], [], [], [], [], [], opts);

%% Display results
fprintf('Solver finished.\n');
fprintf('Total time: %f seconds\n', out.totaltime);
fprintf('Relative recovery error: %f\n', norm(xopt.var{1} - M, 'fro') / norm(M, 'fro'));

% Visualize the results
figure;
subplot(1, 3, 1);
imagesc(M);
title('Original Matrix');
axis square;

subplot(1, 3, 2);
imagesc(M_observed);
title('Observed Matrix');
axis square;

subplot(1, 3, 3);
imagesc(xopt.var{1});
title('Recovered Matrix');
axis square;