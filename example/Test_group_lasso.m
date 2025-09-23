%% Tutorial for Lasso Problem
% This script demonstrates how to solve a Lasso problem using SSNCVX.
% The problem is defined as:
%
% min 0.5*||Ax-b||_2^2 + lambda*||x||_1
%
% We will construct a small-scale problem to illustrate the usage.

clear;clc;
addpath(genpath('../'));

%% Problem Data Construction
% Generate a small random problem.
m = 10; % Number of measurements
n = 20; % Number of features

rng(42);
A = randn(m, n);
x_true = sprandn(n, 1, 0.2); % A sparse signal with 20% non-zero elements
b = A * x_true + 0.1 * randn(m, 1);

lambda = 0.1;

%% pblk setting for the l1-norm penalty
% p(x) = lambda*||x||_1
pblk{1} = struct;
pblk{1}.type = 'topk';
pblk{1}.topk = 10;
pblk{1}.size = n;
pblk{1}.coefficient = lambda;

%% f setting for the least-squares data fitting term
% f(y) = 0.5*||y-b||^2 where y = Ax
f{1} = struct;
f{1}.type = 'square';
f{1}.size = m; % The dimension of y is m
f{1}.coefficient = 0.5;
f{1}.shift = b;

%% Initial guess and solver options
x0 = zeros(n, 1);
Bt = A'; % The solver requires the transpose of A

opts = struct(); % Use default options

%% Call the solver
fprintf('Solving a small Lasso problem...\n');
opts.method = 'iterative';
[xopt, out] = SSNCVX(x0, pblk, Bt, f, [], [], [], [], [], [], [], opts);

%% Display results
fprintf('Solver finished.\n');
fprintf('Total time: %f seconds\n', out.totaltime);

% Plot results for comparison
figure;
stem(x_true, 'b-');
hold on;
stem(xopt.var{1}, 'r--');
legend('True Signal', 'Recovered Signal');
title('Lasso: True vs. Recovered Signal');
grid on;