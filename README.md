<<<<<<< HEAD
## SSNCVX

- SSNCVX (**S**emi-**S**mooth **N**ewton method for **C**on**v**e**x** optimization) is both the name of **our software package** and the **underlying algorithmic framework**.
- It is designed to efficiently solve **convex composite optimization problems**, including those with nonsmooth terms and conic constraints, as well as multi-block structures.
- The **MATLAB** version is open source in this repo and specific details can be found in [this paper](#).
- The ​​**C**++​​ version, is specifically designed to handle conic programming problems and offers improved performance for large-scale instances.
- Visit [our website](#) for more information.

## Problem Formulation

SSNCVX is a general algorithmic framework for solving the following convex composite optimization problem:
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

$$
\begin{aligned}
\min_{\mathbf{x}} \quad & p(\mathbf{x}) + f(\mathcal{B}(\mathbf{x})) + \langle \mathbf{c}, \mathbf{x} \rangle + \frac{1}{2} \langle \mathbf{x}, \mathcal{Q}(\mathbf{x}) \rangle, \\
\text{s.t.} \quad & \mathbf{x} \in \mathcal{P}_1, ~~ \mathcal{A}(\mathbf{x}) \in \mathcal{P}_2,
\end{aligned}
$$

where $p(\mathbf{x})$ is a convex and nonsmooth function,  $\mathcal{A}: \mathcal{X} \rightarrow \mathbb{R}^m$, $\mathcal{B}: \mathcal{X} \rightarrow \mathbb{R}^l$ are linear operators, $f: \mathbb{R}^l \rightarrow \mathbb{R}$ is a convex function, $\mathbf{c} \in \mathcal{X}$, $\mathcal{Q}$ is a positive semidefinite matrix or operator, $\mathcal{P}_1 = \{\mathbf{x} \in \mathcal{X} \mid \texttt{l} \le \mathbf{x} \le \texttt{u}\}$ and $\mathcal{P}_2 = \{\mathbf{x}\in \mathbb{R}^m \mid \texttt{lb} \le \mathbf{x} \le \texttt{ub}\}$. 


Here we list some examples of the problem.

| Problem | Objective Function | Constraints | Function block |
|---------|--------------------|-------------|----------------|
|  | $\underbrace{\langle \mathbf{c},\mathbf{x} \rangle}_{\text{(I)}} + \underbrace{\frac{1}{2}\langle \mathbf{x},\mathcal{Q}(\mathbf{x}) \rangle}_{\text{(II)}} + \underbrace{f(\mathcal{B}(\mathbf{x}))}_{\text{(III)}} + p(\mathbf{x})$ | $\underbrace{\mathbf{x} \in \mathcal{P}_1}_{\text{(IV)}}$, $\underbrace{\mathcal{A}(\mathbf{x}) \in \mathcal{P}_2}_{\text{(V)}}.$ | (I) (II) (III) (IV) (V) |
| LP | $\langle \mathbf{c},\mathbf{x} \rangle $ | $\mathcal{A}(\mathbf{x}) = \mathbf{b},  \mathbf{x} \ge 0$ | (I)(V) |
| SOCP | $\langle \mathbf{c},\mathbf{x} \rangle$ | $\mathcal{A}(\mathbf{x}) = \mathbf{b},  \mathbf{x} \in \mathcal{Q}^n$ | (I)(V) |
| SDP | $\langle \mathbf{C},\mathbf{X} \rangle$ | $\mathcal{A}(\mathbf{X}) = \mathbf{b},  \mathbf{X} \succeq 0$ | (I)(V) |
| SDP with box constraints | $\langle \mathbf{C},\mathbf{X} \rangle $ | $\mathcal{A}(\mathbf{X}) = \mathbf{b},  \mathbf{x} \in \mathcal{P}_1, \mathbf{X} \succeq 0$ | (I)(IV)(V) |
| QP | $\langle \mathbf{x},\mathcal{Q}(\mathbf{x}) \rangle + \langle \mathbf{x},\mathbf{c}\rangle$ | $\texttt{l} \le \mathbf{x} \le \texttt{u}, \mathcal{A}(\mathbf{x}) = \mathbf{b}$ | (I)(II)(IV)(V) |
| QP with $\ell_1$ norm | $\langle \mathbf{x},\mathcal{Q}(\mathbf{x}) \rangle + \lambda \|\|\mathbf{x}\|\|_1$ | $\texttt{l} \le \mathbf{x} \le \texttt{u}, \mathcal{A}(\mathbf{x}) = \mathbf{b}$ | (I)(II)(III)(V) |
| Lasso | $\frac{1}{2}\|\|\mathcal{B}(\mathbf{x})-\mathbf{b}\|\|^2+ \lambda \|\|\mathbf{x}\|\|_1$ | - | (III) |
| Fused Lasso | $\frac{1}{2}\|\|\mathcal{B}(\mathbf{x})-\mathbf{b}\|\|^2 + \lambda_1 \|\|\mathbf{x}\|\|_1 + \lambda_2\|\|D\mathbf{x}\|\|_1$ | - | (III) |
| Group Lasso | $\frac{1}{2}\|\|\mathcal{B}(\mathbf{x})-\mathbf{b}\|\|^2+ \lambda \|\|\mathbf{x}\|\|_1$ |- | (III) |
| Top-k Lasso | $\frac{1}{2}\|\|\mathcal{B}(\mathbf{x})-\mathbf{b}\|\|^2 +  \lambda \sum_{i=1}^k \mathbf{x}_{[i]} $ |- | (III) |
| Low-rank matrix recovery | $\|\|\mathcal{B}(\mathbf{X}) - \mathbf{B}\|\|^2 + \lambda \|\|\mathbf{X}\|\|_*$ | - | (III) |
| Sparse covariance matrix estimation | $- \log(\text{det}(\mathbf{X})) + \text{tr}(\mathbf{XS}) + \lambda \|\|\mathbf{X}\|\|_*$ | - | (I)(III) |
| Sparse PCA | $- \langle \mathbf{L},\mathbf{x} \rangle + \lambda \|\|\mathbf{x}\|\|_1$ | $\text{tr}(\mathbf{x}) =1, \mathbf{x} \succeq 0$ | (III) |
| Basis pursuit | $\|\|\mathbf{x}\|\|_1$ | $ \mathcal{A}(\mathbf{x}) = \mathbf{b} $ | (V) |
| Robust PCA | $\|\|\mathbf{x}_1\|\|_* + \lambda \|\|\mathbf{x}_2\|\|_1$ | $ \mathbf{x}_1 + \mathbf{x}_2 = \mathbf{D} $ | (III)(V) |

## Installation

1.  **Prerequisites:**
    MATLAB R2020a or a later version. (For old versions of MATLAB, there are some issues with mex functions.)

2.  **Download the package:**
    Clone this repository from GitHub.

3.  **Compile mex functions (optional):**
    We have compiled mex functions for SSNCVX, but if there are any issues with the compiled mex functions, you can compile them yourself.
    
    Navigate to the `src\mexfun` directory within the package, delete files with suffix `.mexw64`, `.mexa64`, `.mexmaci64`, `.mexmaca64`. Then open `Installmex_ssm.m` in MATLAB and run it.

4.  **Verify the installation:**
    To verify that the installation is successful, you can run a simple example provided. Navigate to the `example` directory within the package and run the example script.

    If the script runs without any errors, the solver has been installed correctly.

5.  **Add datasets and set path (optional):**
    If you want to use scripts in `example` directory to solve problems mentioned in the paper, please modify `example\addpath_data.m` and add the path of the datasets to the `data_dir` variable.


## Examples

### LP

This tutorial will guide you through solving your first Linear Program (LP) using the SSNCVX solver. An LP is a fundamental type of optimization problem that involves minimizing a linear objective function subject to linear constraints.

--------------
**The Problem**

The standard form for a Linear Program that SSNCVX solves is:

```
minimize    c' * x
subject to  A * x = b
            l <= x <= u
```

Where `x` is the vector of optimization variables. For our first example, we will solve the following simple LP:

- **Objective:** `minimize -x1 - 2*x2`
- **Constraints:** 
  - `2*x1 + x2 <= 4`
  - `x1 + 2*x2 <= 4`
- **Bounds:** `x1 >= 0, x2 >= 0`

To match the solver's required format, we add slack variables (`s1`, `s2`) to convert the inequality constraints into equality constraints. The problem variables become `x = [x1; x2; s1; s2]`, and the problem is formulated as:

- **c:** `[-1; -2; 0; 0]`
- **A:** `[2, 1, 1, 0; 1, 2, 0, 1]`
- **b:** `[4; 4]`
- **l:** `[0; 0; 0; 0]`
- **u:** `[inf; inf; inf; inf]`

--------------
**MATLAB Implementation**

Here is the complete MATLAB script to solve this problem. You can copy and paste this into a new `.m` file and run it.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Define the Problem Data

% Our variables are x = [x1; x2; s1; s2]
n = 4;

% Objective vector
c = [-1; -2; 0; 0];

% Equality constraint matrix and vector
A = sparse([2, 1, 1, 0; ...
            1, 2, 0, 1]);
b = [4; 4];

% Lower and upper bounds
l = zeros(n, 1);
u = inf(n, 1);

%% 2. Set up the Solver Options

opts = struct();
opts.maxits = 1000;       % Maximum number of iterations
opts.stoptol = 1e-6;      % Stopping tolerance
opts.method = 'iterative';

%% 3. Define Problem Structure for SSNCVX

% Initial guess
x0 = zeros(n, 1);

% Objective function (f) - empty for a standard LP
f = {}; 

% Quadratic term (Q) - zero for an LP
Q = sparse(n, n);

% Problem block (pblk) - empty for standard constraints
pblk = {};

% Constraint matrices (At) and bounds (lb, ub) for Ax=b
At = {A'};
lb = b;
ub = b;

%% 4. Run the Solver

[xopt, out] = SSNCVX(x0, pblk, [], f, Q, {c}, {l}, {u}, At, lb, ub, opts);

%% 5. Display the Results

variables = xopt.var{1};

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Optimal objective value: %f\n', out.obj(end));
fprintf('Optimal solution x:\n');
fprintf('  - x1: %.4f\n', variables(1));
fprintf('  - x2: %.4f\n', variables(2));
```

### Lasso

Lasso (Least Absolute Shrinkage and Selection Operator) is a regression analysis method that performs both variable selection and regularization in order to enhance the prediction accuracy and interpretability of the statistical model it produces.

--------------
**The Problem**

The Lasso optimization problem is:

```
minimize (1/2) * ||Ax - b||₂² + λ * ||x||₁
```

Where:
- `A` is the design matrix.
- `b` is the response vector.
- `x` is the vector of coefficients to be determined.
- `λ` (lambda) is the regularization parameter, which controls the sparsity of the solution.

In this tutorial, we will solve a small Lasso problem with randomly generated data.

--------------
**MATLAB Implementation**

Here is the complete MATLAB script to solve this problem. You can copy and paste this into a new `.m` file and run it.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Generate Synthetic Data

m = 50;  % Number of samples
n = 200; % Number of features

% Generate a sparse ground truth signal
x_true = sprandn(n, 1, 0.1); 

% Generate a design matrix and response vector
A = randn(m, n);
b = A * x_true + 0.1 * randn(m, 1); % Add some noise

%% 2. Set up the Solver Options

lambda = 0.1; % Regularization parameter

opts = struct();
opts.maxits = 1000;
opts.stoptol = 1e-6;
opts.method = 'direct'; % Use a direct solver for the subproblem

%% 3. Define Problem Structure for SSNCVX

% Initial guess
x0 = zeros(n, 1);

% Objective function (f)
f = cell(1);
f{1} = struct('type', 'square', 'size', m, 'coefficient', 0.5, 'shift', -b);

% Regularizer (pblk)
pblk = cell(1);
pblk{1} = struct('type', 'l1', 'size', n, 'coefficient', lambda);

% Constraint matrix (B = A')
Bt = A';

%% 4. Run the Solver

[xopt, out] = SSNCVX(x0, pblk, Bt, f, [], [], [], [], [], [], [], opts);

%% 5. Display and Plot the Results

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Optimal objective value: %f\n', out.obj(end));

% Plot the true vs. recovered signal
figure;
stem(x_true, 'b', 'filled', 'DisplayName', 'True Signal');
hold on;
stem(xopt.var{1}, 'r', 'DisplayName', 'Recovered Signal');
legend;
title('Lasso: True vs. Recovered Signal');
xlabel('Feature Index');
ylabel('Coefficient Value');
grid on;
```

### Fused Lasso

The Fused Lasso is an extension of the Lasso that is particularly useful for problems where features have a natural ordering (e.g., time series or genomic data). It encourages sparsity in both the coefficients and their successive differences, leading to piecewise constant solutions.

--------------
**The Problem**

The Fused Lasso optimization problem is:

```
minimize (1/2) * ||Ax - b||₂² + λ₁ * ||x||₁ + λ₂ * ||Dx||₁
```

Where:
- `A` is the design matrix, `b` is the response vector, and `x` is the coefficient vector.
- `λ₁` is the regularization parameter for the L1 norm of `x`, encouraging sparsity.
- `λ₂` is the regularization parameter for the L1 norm of the differences `Dx`.
- `D` is a matrix that computes the differences between adjacent elements of `x`.

--------------
**MATLAB Implementation**

This example demonstrates how to solve a Fused Lasso problem. We'll generate a piecewise constant signal and try to recover it.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Generate Synthetic Data

n = 500; % Signal length

% Create a piecewise constant signal
x_true = zeros(n, 1);
x_true(101:150) = 2;
x_true(251:300) = -1.5;
x_true(401:420) = 1;

% Create a measurement matrix and noisy measurements
m = 200;
A = randn(m, n);
b = A * x_true + 0.1 * randn(m, 1);

%% 2. Set up the Solver Options

lambda1 = 0.1; % Sparsity regularization
lambda2 = 5;   % Fused (total variation) regularization

opts = struct();
opts.maxits = 2000;
opts.stoptol = 1e-5;
opts.method = 'direct';

%% 3. Define Problem Structure for SSNCVX

% Initial guess
x0 = zeros(n, 1);

% Objective function (f)
f = cell(1);
f{1} = struct('type', 'square', 'size', m, 'coefficient', 0.5, 'shift', -b);

% Regularizer (pblk)
pblk = cell(1);
pblk{1} = struct('type', 'fused', 'size', n, 'coefficient', lambda1, 'coefficient2', lambda2);

% Define the difference operator D for the fused term
[Bmap, BTmap] = FLBmap(n);
pblk{1}.Binput = struct('Bmap', Bmap, 'BTmap', BTmap);

% Constraint matrix (B = A')
Bt = A';

%% 4. Run the Solver

[xopt, out] = SSNCVX(x0, pblk, Bt, f, [], [], [], [], [], [], [], opts);

%% 5. Display and Plot the Results

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Optimal objective value: %f\n', out.obj(end));

% Plot the true vs. recovered signal
figure;
plot(x_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Signal');
hold on;
plot(xopt.var{1}, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Recovered Signal');
legend;
title('Fused Lasso: True vs. Recovered Signal');
xlabel('Feature Index');
ylabel('Coefficient Value');
grid on;
```

### QP

Quadratic Programming (QP) is a fundamental problem in mathematical optimization. It involves minimizing a quadratic objective function over a feasible region defined by linear equality and inequality constraints. QP has a wide range of applications, including portfolio optimization, support vector machines, and control theory.

--------------
**The Problem**

The standard form of a Quadratic Program that SSNCVX can solve is:

```
minimize    (1/2) * x' * Q * x + c' * x
subject to  A * x = b
            l <= x <= u
```

Where:
- `x` is the vector of optimization variables.
- `Q` is a symmetric matrix (the quadratic form matrix).
- `c` is a vector of linear coefficients.
- `A*x = b` represents the linear equality constraints.
- `l <= x <= u` represents the box constraints (lower and upper bounds) on the variables.

To illustrate how to solve a QP with this solver, let's consider a simple example. We aim to find the point `x` in 2D space that is closest to the origin, subject to the constraints that `x1 + x2 = 1` and `0 <= x1, x2 <= 0.7`.

This can be formulated as the following QP:

- **Objective:** `minimize (1/2) * (x1² + x2²)`, which corresponds to `Q = eye(2)` and `c = zeros(2,1)`.
- **Equality Constraint:** `x1 + x2 = 1`, which corresponds to `A = [1, 1]` and `b = 1`.
- **Box Constraints:** `l = [0; 0]` and `u = [0.7; 0.7]`.

--------------
**MATLAB Implementation**

Here is the complete MATLAB script to define and solve this problem using the SSNCVX solver.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Define the Problem Data

n = 2; % Number of variables
Q = speye(n);
c = zeros(n, 1);
A = sparse([1, 1]);
b = 1;
l = zeros(n, 1);
u = 0.7 * ones(n, 1);

%% 2. Set up the Solver Options

opts = struct();
opts.maxits = 1000;       % Maximum number of iterations
opts.stoptol = 1e-6;      % Stopping tolerance
opts.method = 'iterative';

%% 3. Define Problem Structure for SSNCVX

% Initial guess
x0 = zeros(n, 1);

% Objective function (f)
f = cell(1);
f{1} = struct('type', 'square', 'size', n, 'coefficient', 0.5);

% Constraints (pblk)
pblk = cell(1);
pblk{1} = struct('type', 'l2', 'size', n);

% Constraint matrices
At = {A'};

% Constraint bounds
lb = b;
ub = b;

%% 4. Run the Solver

[xopt, out] = SSNCVX(x0, pblk, [], f, Q, {c}, {l}, {u}, At, lb, ub, opts);

%% 5. Display the Results

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Optimal objective value: %f\n', out.obj(end));
fprintf('Optimal solution x:\n');
disp(xopt.var{1});
```

### SOCP

Second-Order Cone Programming (SOCP) is a class of convex optimization problems that includes Linear Programming (LP) and Quadratic Programming (QP) as special cases. SOCP problems involve minimizing a linear function over an intersection of affine sets and second-order cones.

--------------
**The Problem**

The standard form of an SOCP is:

```
minimize    f' * x
subject to  ||A_i * x - b_i||_2 <= c_i' * x + d_i,  for i = 1,...,m
            F * x = g
```

Where `x` is the optimization variable. The constraints `||A_i * x - b_i||_2 <= c_i' * x + d_i` are called second-order cone constraints.

The SSNCVX solver handles problems in a standard conic form:

```
minimize    <C, x>
subject to  A * x = b
            x in K
```

Where `K` is a product of simple cones (e.g., non-negative, second-order).

Let's solve a simple SOCP. We want to find the point `(x, y)` in the first quadrant of the unit circle that maximizes the sum `x + y`.

This can be formulated as:

- **Objective:** `maximize x + y` (or `minimize -x - y`)
- **Constraints:**
    - `sqrt(x² + y²) <= 1` (point is inside the unit circle)
    - `x >= 0, y >= 0` (point is in the first quadrant)

To fit this into the solver's format, we introduce a variable `t` and formulate the problem in 3D space `(t, x, y)`.

- **Variable:** `z = [t; x; y]`
- **Objective:** `minimize -x - y`, which is `[0, -1, -1]' * z`
- **Constraints:**
    - `sqrt(x² + y²) <= t` (This is the definition of the second-order cone)
    - `t = 1` (We are interested in the slice of the cone at `t=1`)
    - `x >= 0, y >= 0`

The solution should be `x = y = 1/sqrt(2)`. 

--------------
**MATLAB Implementation**

Here is the MATLAB script to define and solve this problem.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Define the Problem Data

% Objective function: minimize -x - y
C = {[0; -1; -1]};

% Linear constraint: t = 1
At = {[1; 0; 0]};
lb = 1;
ub = 1;

% Box constraints: x >= 0, y >= 0
% For t, we have no bounds (-inf, inf)
l = {[-inf; 0; 0]};
u = {[inf; inf; inf]};

%% 2. Set up the Solver Options

opts = struct();
opts.maxits = 100;
opts.stoptol = 1e-7;

%% 3. Define Problem Structure for SSNCVX

% Initial guess (not needed, can be empty)
x0 = [];

% Define the second-order cone constraint on our variable z = [t; x; y]
pblk = cell(1);
pblk{1} = struct('type', 'q', 'size', 3);

%% 4. Run the Solver

[xopt, out] = SSNCVX(x0, pblk, [], [], [], C, l, u, At, lb, ub, opts);

%% 5. Display the Results

solution = xopt.var{1};

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Optimal objective value: %f\n', out.obj(end));
fprintf('Optimal solution (t, x, y):\n');
disp(solution);
```

### SPCA

Sparse Principal Component Analysis (SPCA) is a modern variant of PCA that aims to find principal components with a small number of non-zero entries. This is highly desirable in high-dimensional settings (e.g., genomics, finance) where interpretability is key. By forcing the loadings to be sparse, we can identify which original features are most important for explaining the variance.

--------------
**The Problem**

One common formulation for finding a single sparse principal component is:

```
maximize    x' * A * x
subject to  ||x||_2 <= 1
            ||x||_1 <= k
```

Where:
- `A` is the sample covariance matrix of the data.
- `x` is the principal component (loading vector) we are looking for.
- `||x||_2 <= 1` is the standard PCA unit norm constraint.
- `||x||_1 <= k` is the sparsity-inducing L1-norm constraint, where `k` controls the sparsity.

This problem can be reformulated and solved as a semidefinite program (SDP), but for this solver, we can often tackle it more directly using its specialized structure.

--------------
**MATLAB Implementation**

In this tutorial, we will generate a synthetic dataset where the first principal component is sparse. Then, we will use the solver to recover it.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Generate Synthetic Data

p = 100; % Number of features
n = 80;  % Number of samples

% Create a sparse principal component
v1 = zeros(p, 1);
v1([1, 10, 20, 30, 40]) = 1;
v1 = v1 / norm(v1);

% Generate data with this component
Data = randn(n, p);
Data = Data - (Data * v1) * v1'; % Project out v1
Data = Data + (20 * randn(n, 1)) * v1'; % Add strong v1 component

% Compute the sample covariance matrix
A = cov(Data);

%% 2. Set up the Solver Options

lambda = 0.5; % Sparsity-controlling parameter

opts = struct();
opts.maxits = 500;
opts.stoptol = 1e-5;

%% 3. Define Problem Structure for SSNCVX

% We want to maximize x'Ax, which is equivalent to minimizing -x'Ax
% This is a quadratic objective.
Q = -A;

% Initial guess
x0 = randn(p, 1); x0 = x0/norm(x0);

% Define the constraints
pblk = cell(2,1);
% L2-norm ball constraint: ||x||_2 <= 1
pblk{1} = struct('type', 'l2', 'size', p, 'ub', 1);
% L1-norm ball constraint (sparsity): ||x||_1 <= k
% We use the 'l1' block for this.
pblk{2} = struct('type', 'l1', 'size', p, 'ub', sqrt(5) + 0.5); % k = sqrt(non-zero elements)

%% 4. Run the Solver

% No linear terms or constraints, so most arguments are empty.
[xopt, out] = SSNCVX({x0}, pblk, [], [], Q, [], [], [], [], [], [], opts);

solution = xopt.var{1};

%% 5. Display the Results

% Clean up small values for better visualization
solution(abs(solution) < 1e-4) = 0;

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Number of non-zero elements in recovered component: %d\n', nnz(solution));

% Plot the true and recovered components
figure;
stem(v1, 'b', 'DisplayName', 'True Component');
hold on;
stem(solution, 'r--', 'DisplayName', 'Recovered Component');
legend;
title('Sparse Principal Component Recovery');
```

### LRMC

Low-Rank Matrix Completion (LRMC) is the task of recovering a matrix from a small subset of its entries, under the assumption that the original matrix is low-rank. This has applications in recommendation systems (predicting user ratings), image inpainting, and system identification.

--------------
**The Problem**

The most common formulation for LRMC is the nuclear norm minimization problem:

```
minimize ||X||_*
subject to P_Ω(X) = P_Ω(M)
```

Where:
- `X` is the matrix to be recovered.
- `||X||_*` is the nuclear norm of `X` (the sum of its singular values), which is a convex proxy for the rank of `X`.
- `M` is the original (unknown) low-rank matrix.
- `Ω` is the set of indices of the observed entries.
- `P_Ω` is a projection operator that keeps the entries in `Ω` and sets others to zero.

We can rewrite this as a least-squares problem with nuclear norm regularization:

```
minimize (1/2) * ||P_Ω(X) - P_Ω(M)||_F² + λ * ||X||_*
```

--------------
**MATLAB Implementation**

In this tutorial, we will generate a low-rank matrix, sample a fraction of its entries, and then use the SSNCVX solver to recover the full matrix.

```matlab
% Add the solver to your MATLAB path
addpath(genpath('path/to/your/SSNCVX'));

%% 1. Generate Synthetic Low-Rank Data

n1 = 100; % Matrix rows
n2 = 100; % Matrix columns
rank = 5;  % Rank of the matrix

% Create a low-rank matrix
M = randn(n1, rank) * randn(rank, n2);

% Create a sampling mask (observe 50% of entries)
sampling_ratio = 0.5;
Omega = rand(n1, n2) < sampling_ratio;

% Observed data
M_observed = M .* Omega;

%% 2. Set up the Solver Options

lambda = 10; % Regularization parameter

opts = struct();
opts.maxits = 500;
opts.stoptol = 1e-4;
opts.method = 'iterative'; % Use an iterative solver for the subproblem

%% 3. Define Problem Structure for SSNCVX

% Initial guess (a matrix of zeros)
X0 = zeros(n1, n2);

% Regularizer (pblk) - Nuclear norm
pblk = cell(1);
pblk{1} = struct('type', 'nuclear', 'size', [n1, n2], 'coefficient', lambda);

% Objective function (f) - Least squares data fidelity
f = cell(1);
f{1} = struct('type', 'square', 'size', [n1, n2], 'coefficient', 0.5, 'shift', -M_observed);

% Sampling operator (B)
B.Bmap = @(u) u .* Omega;
B.BTmap = @(y) y .* Omega;
B.out_size = [n1, n2];

%% 4. Run the Solver

[Xopt, out] = SSNCVX(X0, pblk, B, f, [], [], [], [], [], [], [], opts);

recovered_matrix = Xopt.var{1};

%% 5. Display the Results

% Calculate the relative recovery error
recovery_error = norm(recovered_matrix - M, 'fro') / norm(M, 'fro');

fprintf('Solver terminated in %d iterations.\n', out.iter);
fprintf('Relative recovery error: %.4f\n', recovery_error);

% Visualize the matrices
figure;
subplot(1, 3, 1); imagesc(M); title('Original Matrix'); axis off;
subplot(1, 3, 2); imagesc(M_observed); title('Observed Entries'); axis off;
subplot(1, 3, 3); imagesc(recovered_matrix); title('Recovered Matrix'); axis off;
colormap('gray');
```

## References

This software is released under the [MIT License](https://opensource.org/licenses/MIT).

If you use SSNCVX in your research, please cite [this paper](#).

```bibtex
@article{deng2025augmented,
}
```

If you want to contact us for collaboration or other inquiries, please reach out via email at [to_be_determined@pku.edu.cn]


=======
# SSNCVX
>>>>>>> e1666d73c4c9dc5472e9b0bb7ce03fae4cb92407
