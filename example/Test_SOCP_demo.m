%% Tutorial for Second-Order Cone Programming (SOCP)
% This script demonstrates how to solve an SOCP problem using SSNCVX.
% The problem is defined as:
%
% min c'*x
% s.t. A*x = b
%      x in K
%
% where K is a product of second-order cones.
% K = K_1 x K_2 x ... x K_L
% K_j = {x_j = (t_j; v_j) | ||v_j||_2 <= t_j}
%
% We will solve a simple SOCP problem:
% min x1 + x2 + 2*x3
% s.t. x1 + 2*x2 + 3*x3 = 1
%      ||[x2; x3]|| <= x1

clear;clc;
addpath(genpath('../'));

%% Problem Data Construction
% min x1 + x2 + 2*x3  => c = [1; 1; 2]
C = {[1; 1; 2]};
n = 3;

% s.t. x1 + 2*x2 + 3*x3 = 1 => A = [1 2 3], b = 1
At = {[1; 2; 3]};
lb = 1;
ub = 1;
m = 1;

% s.t. ||[x2; x3]|| <= x1
% This defines a second-order cone constraint on x = [x1; x2; x3].
% The cone is of size 3.
pblk{1} = struct;
pblk{1}.type = 'q'; % 'q' for quadratic or second-order cone
pblk{1}.size = 3;
pblk{1}.coefficient = 1; % Not used in this context, but required

% The K structure is also needed by the solver options
K{1} = struct;
K{1}.type = 'q';
K{1}.size = 3;

%% Initial guess and solver options
opts = struct();
opts.method = 'direct';
opts.K = K;
opts.m = m;

%% Call the solver
fprintf('Solving a small SOCP problem...\n');
% For SOCP, the structure of the call is similar to QP.
% We pass empty placeholders for unused arguments.
[xopt, out] = SSNCVX([], pblk, [], [], [], C, [], [], At, lb, ub, opts);

%% Display results
fprintf('Solver finished.\n');
fprintf('Total time: %f seconds\n', out.totaltime);

% The optimal solution is in xopt.var{1}
sol = xopt.var{1};
fprintf('Optimal solution x = [%f, %f, %f]\n', sol(1), sol(2), sol(3));

% Verify constraints
fprintf('Constraint violation ||Ax-b||: %e\n', norm(At{1}'*sol - lb));
fprintf('SOCP constraint: x1 = %f, ||[x2;x3]|| = %f\n', sol(1), norm(sol(2:3)));
if sol(1) >= norm(sol(2:3)) - 1e-6
    fprintf('SOCP constraint is satisfied.\n');
else
    fprintf('SOCP constraint is VIOLATED.\n');
end