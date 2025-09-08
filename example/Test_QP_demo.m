%% Tutorial for Quadratic Programming (QP)
% This script demonstrates how to solve a Quadratic Program using SSNCVX.
% The problem is defined as:
%
% min 0.5*x'*Q*x + c'*x
% s.t. Ax = b
%      lb <= x <= ub
%
% We will construct a small-scale problem to illustrate the usage.

clear;clc;
addpath(genpath('../'));

%% Problem Data Construction
n = 10; % Number of variables
m = 5;  % Number of equality constraints

% Generate a positive semi-definite matrix Q
rng(42);
Q_half = randn(n, n);
Q = Q_half' * Q_half; 
c = randn(n, 1);

A = randn(m, n);
x_true = rand(n, 1);
b = A * x_true;

lb_val = 0;
ub_val = Inf;

%% pblk setting for box constraints
% p(x) for lb <= x <= ub
pblk{1} = struct;
pblk{1}.type = 'box';
pblk{1}.size = n;
pblk{1}.l = lb_val * ones(n, 1);
pblk{1}.u = ub_val * ones(n, 1);

%% Q and C setting for the quadratic objective
% The objective is 0.5*x'*Q*x + c'*x
Q_in = Q;
C_in = {c};

%% At, lb, ub setting for equality constraints
% Ax = b
At_in = {A'};
lb_eq = b;
ub_eq = b;

%% Initial solver options
opts = struct(); % Use default options
opts.m = m;

%% Call the solver
fprintf('Solving a small QP problem...\n');
[xopt, out] = SSNCVX([], pblk, [], [], Q_in, C_in, [], [], At_in, lb_eq, ub_eq, opts);

%% Display results
fprintf('Solver finished.\n');
fprintf('Total time: %f seconds\n', out.totaltime);
fprintf('Constraint violation ||Ax-b||: %e\n', norm(A*xopt.var{1} - b));

% The optimal solution is in xopt.var{1}
fprintf('Optimal x (first 5 elements): \n');
disp(xopt.var{1}(1:5));