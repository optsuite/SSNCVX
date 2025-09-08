%% Tutorial for Fused Lasso Problem
% This script demonstrates how to solve a Fused Lasso problem using SSNCVX.
% The problem is defined as:
%
% min 0.5*||Ax-b||_2^2 + lambda1*||x||_1 + lambda2*sum(|x_i - x_{i-1}|)
%
% We will construct a small-scale problem to illustrate the usage.

clear;clc;
addpath(genpath('../'));

%% Problem Data Construction
% Generate a small random problem.
m = 5;  % Number of measurements
n = 10; % Number of features

rng(42);
A = randn(m, n);
b = randn(m, 1);
x_true = zeros(n, 1);
x_true(3:6) = 1; % A sparse and piecewise-constant signal
b = A * x_true + 0.1 * randn(m, 1);

lambda1 = 0.1;
lambda2 = 0.1;

%% pblk setting for the Fused Lasso penalty
% p(x) = lambda1*||x||_1 + lambda2*sum(|x_i - x_{i-1}|)
pblk{1} = struct;
pblk{1}.type = 'fused';
pblk{1}.size = n;
pblk{1}.coefficient = lambda1;
pblk{1}.coefficient2 = lambda2;
[Bmap,BTmap] = FLBmap(n);
pblk{1}.Binput = struct();
pblk{1}.Binput.Bmap = Bmap;
pblk{1}.Binput.BTmap = BTmap;

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

opts = struct();
opts.method = 'direct'; % Use direct method rather than iterative

%% Call the solver
fprintf('Solving a small Fused Lasso problem...\n');
[xopt, out] = SSNCVX(x0, pblk, Bt, f, [], [], [], [], [], [], [], opts);

%% Display results
fprintf('Solver finished.\n');
fprintf('Total time: %f seconds\n', out.totaltime);

% Plot results for comparison
figure;
plot(x_true, 'b-', 'LineWidth', 2);
hold on;
plot(xopt.var{1}, 'r--', 'LineWidth', 2);
legend('True Signal', 'Recovered Signal');
title('Fused Lasso: True vs. Recovered Signal');
grid on;