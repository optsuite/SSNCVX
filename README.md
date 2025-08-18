## SSNCVX

- SSNCVX (**S**emi-**S**mooth **N**ewton method for **C**on**v**e**x** optimization) is both the name of **our software package** and the **underlying algorithmic framework**.
- It is designed to efficiently solve **convex composite optimization problems**, including those with nonsmooth terms and conic constraints, as well as multi-block structures.
- The **MATLAB** version is open source in this repo and specific details can be found in [this paper](#).
- The ​​**C**++​​ version, SSNCVXpp, is specifically designed to handle conic programming problems and offers improved performance for large-scale instances. We have C++ API and Python API.
- Visit [our website](#) for more information.

## Problem Formulation

SSNCVX is a general algorithmic framework for solving the following convex composite optimization problem:

$$
\begin{aligned}
\min_{\mathbf{x}} \quad & p(\mathbf{x}) + f(\mathcal{B}(\mathbf{x})) + \langle \mathbf{c}, \mathbf{x} \rangle + \frac{1}{2} \langle \mathbf{x}, \mathcal{Q}(\mathbf{x}) \rangle, \\
\text{s.t.} \quad & \mathbf{x} \in \mathcal{P}_1, ~~ \mathcal{A}(\mathbf{x}) \in \mathcal{P}_2,
\end{aligned}
$$

where $p(\mathbf{x})$ is a convex and nonsmooth function,  $\mathcal{A}: \mathcal{X} \rightarrow \mathbb{R}^m$, $\mathcal{B}: \mathcal{X} \rightarrow \mathbb{R}^l$ are linear operators, $f: \mathbb{R}^l \rightarrow \mathbb{R}$ is a convex function, $\mathbf{c} \in \mathcal{X}$, $\mathcal{Q}$ is a positive semidefinite matrix or operator, $\mathcal{P}_1 = \{\mathbf{x} \in \mathcal{X} \mid \texttt{l} \le \mathbf{x} \le \texttt{u}\}$ and $\mathcal{P}_2 = \{\mathbf{x}\in \mathbb{R}^m \mid \texttt{lb} \le \mathbf{x} \le \texttt{ub}\}$. The choices of $p(\mathbf{x})$ provide flexibility to handle many types of problems. While the proposed model focuses on a single variable $\mathbf{x}$, it is indeed capable of solving the following more general problem with $N$ blocks of variables and shifting terms $\mathbf{b}_{1,i}$, $\mathbf{b}_{2,i}$ ($i =1,\cdots,N$):

$$
\begin{aligned}
\min_{\mathbf{x}_i} \quad & \sum_{i=1}^N p_i(\mathbf{x}_i - \mathbf{b}_{1,i}) + \sum_{i=1}^N f_i(\mathcal{B}_i(\mathbf{x}) - \mathbf{b}_{2,i}) + \sum_{i=1}^N \langle \mathbf{c}_i, \mathbf{x}_i \rangle + \sum_{i=1}^N \frac{1}{2} \langle \mathbf{x}_i, \mathcal{Q}_i(\mathbf{x}_i) \rangle, \\
\text{s.t.} \quad & \mathbf{x}_i \in \mathcal{P}_{1,i}, ~~ \sum_{i=1}^N \mathcal{A}_i (\mathbf{x}_i) \in \mathcal{P}_{2}, \quad i =1,\cdots,N,
\end{aligned}
$$

where $p_i,f_i,\mathbf{c}_i,\mathcal{Q}_i, \mathcal{P}_{1,i}$ and $\mathcal{P}_{2}$ satisfy the same assumptions as in the single-variable model. These models have widespread applications in engineering, image processing, and machine learning. We refer readers to Boyd et al. (2011), Anjos et al. (2011), and Wolkowicz et al. (2012) for concrete applications.

Here we list some examples of the problem.

| Problem | Objective Function | Constraints | Function block |
|---------|--------------------|-------------|----------------|
|  | $\underbrace{\langle \bm{c},\bm{x} \rangle}_{\text{(I)}} + \underbrace{\frac{1}{2}\langle \bm{x},\mathcal{Q}(\bm{x}) \rangle}_{\text{(II)}} + \underbrace{f(\mathcal{B}(\bm{x}))}_{\text{(III)}} + p(\bm{x})$ | $\underbrace{\bm{x} \in \mathcal{P}_1}_{\text{(IV)}}$, $\underbrace{\mathcal{A}(\bm{x}) \in \mathcal{P}_2}_{\text{(V)}}.$ | (I) (II) (III) (IV) (V) |
| LP | $\langle \bm{c},\bm{x} \rangle $ | $\mathcal{A}(\bm{x}) = \bm{b},  \bm{x} \ge 0$ | (I)(V) |
| SOCP | $\langle \bm{c},\bm{x} \rangle$ | $\mathcal{A}(\bm{x}) = \bm{b},  \bm{x} \in \mathcal{Q}^n$ | (I)(V) |
| SDP | $\langle \bm{C},\bm{X} \rangle$ | $\mathcal{A}(\bm{X}) = \bm{b},  \bm{X} \succeq 0$ | (I)(V) |
| SDP with box constraints | $\langle \bm{C},\bm{X} \rangle $ | $\mathcal{A}(\bm{X}) = \bm{b},  \bm{x} \in \mathcal{P}_1, \bm{X} \succeq 0$ | (I)(IV)(V) |
| QP | $\langle \bm{x},\mathcal{Q}(\bm{x}) \rangle + \langle \bm{x},\bm{c}\rangle$ | $\texttt{l} \le \bm{x} \le \texttt{u}, \mathcal{A}(\bm{x}) = \bm{b}$ | (I)(II)(IV)(V) |
| QP with $\ell_1$ norm | $\langle \bm{x},\mathcal{Q}(\bm{x}) \rangle + \lambda \|\|\bm{x}\|\|_1$ | $\texttt{l} \le \bm{x} \le \texttt{u}, \mathcal{A}(\bm{x}) = \bm{b}$ | (I)(II)(III)(V) |
| Lasso | $\frac{1}{2}\|\|\mathcal{B}(\bm{x})-\bm{b}\|\|^2+ \lambda \|\|\bm{x}\|\|_1$ | - | (III) |
| Fused Lasso | $\frac{1}{2}\|\|\mathcal{B}(\bm{x})-\bm{b}\|\|^2 + \lambda_1 \|\|\bm{x}\|\|_1 + \lambda_2\|\|D\bm{x}\|\|_1$ | - | (III) |
| Group Lasso | $\frac{1}{2}\|\|\mathcal{B}(\bm{x})-\bm{b}\|\|^2+ \lambda \|\|\bm{x}\|\|_1$ |- | (III) |
| Top-k Lasso | $\frac{1}{2}\|\|\mathcal{B}(\bm{x})-\bm{b}\|\|^2 +  \lambda \sum_{i=1}^k \bm{x}_{[i]} $ |- | (III) |
| Low-rank matrix recovery | $\|\|\mathcal{B}(\bm{X}) - \bm{B}\|\|^2 + \lambda \|\|\bm{X}\|\|_*$ | - | (III) |
| Sparse covariance matrix estimation | $- \log(\text{det}(\bm{X})) + \text{tr}(\bm{XS}) + \lambda \|\|\bm{X}\|\|_*$ | - | (I)(III) |
| Sparse PCA | $- \langle \bm{L},\bm{x} \rangle + \lambda \|\|\bm{x}\|\|_1$ | $\text{tr}(\bm{x}) =1, \bm{x} \succeq 0$ | (III) |
| Basis pursuit | $\|\|\bm{x}\|\|_1$ | $ \mathcal{A}(\bm{x}) = \bm{b} $ | (V) |
| Robust PCA | $\|\|\bm{x}_1\|\|_* + \lambda \|\|\bm{x}_2\|\|_1$ | $ \bm{x}_1 + \bm{x}_2 = \bm{D} $ | (III)(V) |

## Algorithms

* **Problem Formulation**: The algorithm starts by considering a multi-block convex optimization problem and formulating its dual.

$$
\begin{aligned}
    &\min_{\bm{y},\bm{z},\bm{s},\bm{r},\bm{v}}  \quad  \delta_{\mathcal{P}_2}^*(-\bm{y}) + f^*(\bm{-z}) +  p^*(-\bm{s}) + \frac{1}{2} \left< \mathcal{Q}\bm{v},\bm{v}\right> + \delta_{\mathcal{P}_1}^*(-\bm{r}), \\
    &\quad \text{s.t.} \quad  \mathcal{A}^*(\bm{y}) + \mathcal{B}^*\bm{z} + \bm{s} - \mathcal{Q}\bm{v} + \bm{r} = \bm{c}.
\end{aligned}
$$

* **Augmented Lagrangian**: Introduce slack variables and formulate the augmented Lagrangian function $\mathcal{L}_{\sigma}$ of the equivalent problem.

$$
\begin{aligned}
& \mathcal{L}_{\sigma}(\bm{y},\bm{s},\bm{z},\bm{r},\bm{v},\bm{o},\bm{q},\bm{t},\bm{x}_1,\bm{x}_2,\bm{x}_3,\bm{x}_4) = \delta^*_{\mathcal{P}_2}(-\bm{o}) + f^*(-\bm{q}) +p^*(-\bm{s}) - \left< \bm{b}_1,\bm{s}\right> + \frac{1}{2}\left< \mathcal{Q}(\bm{v}),\bm{v}\right>  \\
& \qquad + \delta_{\mathcal{P}_1}^*(-\bm{t})  + \frac{\sigma}{2}\left(\|\bm{o}-\bm{y} +  \frac{1}{\sigma}\bm{x}_1 \|_{\mathrm{F}}^2 + \|\bm{q}-\bm{z}+ \frac{1}{\sigma}\bm{x}_2 \|_{\mathrm{F}}^2 + \|\bm{t}-\bm{r}+ \frac{1}{\sigma}\bm{x}_3 \|^2 \right) \\
& \qquad + \frac{\sigma}{2}(\| \mathcal{A}^*(\bm{y}) + \mathcal{B}^*\bm{z} + \bm{s} - \mathcal{Q}\bm{v} + \bm{r} - \bm{c} + \frac{1}{\sigma}\bm{x}_4 \|_{\mathrm{F}}^2) - \frac{1}{2\sigma}\sum_{i=1}^4\|\bm{x}_i\|^2 .
\end{aligned}
$$

* **Saddle Point Problem**: By minimizing the augmented Lagrangian with respect to some variables, the problem is reformulated as a differentiable saddle point problem based on the Moreau envelope.
    $$
    \small
    \begin{aligned}
    \Phi_{\sigma}(\bm{w}) &  =  \underbrace{p^*(\text{prox}_{p^*/\sigma}(\bm{x}_4/\sigma -\mathcal{A}^*(\bm{y}) - \mathcal{B}^*\bm{z}  - \mathcal{Q}\bm{v} -\bm{r} + \bm{c} ) ) + \frac{1}{2\sigma}\|\text{prox}_{\sigma p} (\bm{x}_4 + \sigma(\mathcal{A}^*(\bm{y}) + \mathcal{B}^*\bm{z}  - \mathcal{Q}\bm{v} + \bm{r} - \bm{c})) \|^2}_{ \text{Moreau~envelope~} p^* }
    \\
    & \quad   +  \underbrace{\delta_{\mathcal{P}_1}^*(\text{prox}_{\delta^*_{\mathcal{P}_1} }(\bm{x}_3/\sigma - \bm{t} ) ) + \frac{1}{2\sigma}\|\Pi_{\mathcal{P}_1}(\bm{x}_3 -\sigma\bm{r}  ) \|^2}_{ \text{Moreau~envelope~} \delta^*_{\mathcal{P}_1} } +\underbrace{ \delta^*_{\mathcal{P}_2}(\text{prox}_{\delta^*_{\mathcal{P}_2}/\sigma}(\bm{x}_1/\sigma -\bm{y})) + \frac{1}{2\sigma}\|\Pi_{\mathcal{P}_2}(\bm{x}_1 - \sigma \bm{y}) \|^2}_{ \text{Moreau~envelope~} \delta^*_{\mathcal{P}_2}  }     \\
    & \quad +  \underbrace{f^*(\text{prox}_{f^*/\sigma}(\bm{x}_2/\sigma -\bm{z} )) + \frac{1}{2\sigma} \|\text{prox}_{\sigma f}(\bm{x}_2 -\sigma \bm{z}  )  \|^2}_{ \text{Moreau~envelope~} f^* } + \frac{1}{2}\left<\mathcal{Q}\bm{v}, \bm{v}\right> - \frac{1}{2\sigma}\sum_{i=1}^4 \|\bm{x}_i \|^2.
    \end{aligned}
    $$
    $$
        \min_{\bm{y},\bm{z},\bm{r},\bm{v}} \max_{\bm{x}_1,\bm{x}_2,\bm{x}_3,\bm{x}_4 } \Phi(\bm{y},\bm{z},\bm{r},\bm{v};\bm{x}_1,\bm{x}_2,\bm{x}_3,\bm{x}_4).
    $$

* **Nonlinear System**: The saddle point problem is equivalent to solving a system of nonlinear equations $F(\bm{w}) = 0$, where $F(\bm{w})$ is defined by the gradient of $\Phi$.
    $$
         F(\bm{w}) =
         \begin{pmatrix}
              \nabla_{\bm{y}} \Phi(\bm{w});\nabla_{\bm{z}} \Phi(\bm{w});\nabla_{\bm{r}} \Phi(\bm{w});\nabla_{\bm{v}} \Phi(\bm{w});
             - \nabla_{\bm{x}_1} \Phi(\bm{w});
             - \nabla_{\bm{x}_2} \Phi(\bm{w});
             - \nabla_{\bm{x}_3} \Phi(\bm{w});
             - \nabla_{\bm{x}_4} \Phi(\bm{w})
         \end{pmatrix} = 0.
    $$

* **Semismooth Newton Method**: An element $J^k$ from the generalized Jacobian $\hat{\partial} F(\bm{w}^k)$ is chosen to construct a Newton-like system. The search direction $\bm{d}^{k, i}$ is then computed by solving the following linear system:
    $$
    (J^k + \tau_{k,i} \mathcal{I}) \bm{d}^{k, i} = -  F(\bm{w}^k) + \bm{\varepsilon}^k,
    $$
    This direction is used to update the variables for the next iteration.
    $$
    \bar{\bm{w}}^{k,i}  = \bm{w}^k + \bm{d}^{k,i}.
    $$

## Installation

1.  **Prerequisites:**
    MATLAB R2020a or a later version. (For old versions of MATLAB, there are some issues with mex functions.)

2.  **Download the package:**
    Clone this repository from GitHub.

3.  **Compile mex functions (optional):**
    We have compiled mex functions for SSNCVX, but if there are any issues with the compiled mex functions, you can compile them yourself.
    
    Navigate to the `src\mexfun` directory within the package, delete files with suffix `.mexw64`, `.mexa64`, `.mexmaci64`, `.mexmaca64`. Then open `Installmex_ssm.m` in MATLAB and run it.

4.  **Verify the installation:**
    To verify that the installation is successful, you can run a simple example provided. Navigate to the `example` directory within the package and run the example script, `Test_simple_example.m`.

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

## Benchmarks

- In the following tables, `η` presents a metric of accuracy. Visit [our paper](#) to know more details about the benchmark results of SSNCVX against other solvers.
- For certain problems, parameters are tuned on a case-by-case basis to achieve optimal performance. We provide detailed [documentation](#) on how to configure these parameters for new problems.

--------------------------
**Lasso**

Note that for the datasets `pyrim`, `triazines`, `abalone`, `bodyfat`, `housing`, `mpg`, and `space_ga`, we expand their original features by using polynomial basis functions over
those features. For example, the last digit in `pyrim5` indicates that an order 5 polynomial is used to generate the basis functions. This naming convention is also used
in the rest of the expanded data sets.

| id | nnz | SSNCVX | | SSNAL | | SLEP | | ADMM | |
|---|---|---|---|---|---|---|---|---|---|
| | | η | time | η | time | η | time | η | time |
| uci_CT | 13 | 7.6e-7 | **0.64** | 4.4e-13 | 0.86 | 2.2e-2 | 35.95 | 7.7e-3 | 46.02 |
| log1p.E2006.train | 5 | 5.4e-7 | **17.3** | 1.5e-11 | 36.0 | 2.0e-2 | 1850.15 | 1.2e-1 | 3604.34 |
| E2006.test | 1 | 2.2e-11 | **0.17** | 4.3e-10 | 0.28 | 7.5e-12 | 1.11 | 7.9e-7 | 428.64 |
| log1p.E2006.test | 8 | 3.3e-8 | **2.83** | 2.5e-10 | 5.12 | 4.8e-2 | 447.56 | 1.2e-1 | 3603.64 |
| pyrim5 | 72 | 4.2e-16 | **1.82** | 5.7e-8 | 2.16 | 2.4e-2 | 106.09 | 1.5e-3 | 3600.52 |
| triazines4 | 519 | 2.6e-13 | **10.64** | 3.4e-9 | 11.23 | 8.3e-2 | 246.11 | 9.7e-3 | 3603.99 |
| abalone7 | 24 | 4.6e-11 | **0.75** | 1.8e-9 | 1.06 | 2.5e-3 | 34.57 | 3.7e-4 | 540.27 |
| bodyfat7 | 2 | 4.8e-13 | **0.79** | 1.4e-8 | 1.08 | 1.9e-6 | 28.10 | 8.4e-4 | 3609.63 |
| housing7 | 158 | 5.1e-13 | **1.83** | 6.3e-9 | 1.74 | 1.3e-2 | 46.60 | 1.1e-2 | 3601.26 |
| mpg7 | 47 | 4.4e-16 | **0.10** | 1.5e-8 | 0.14 | 7.4e-5 | 0.69 | 1.0e-6 | 63.41 |
| spacega9 | 14 | 4.7e-15 | **0.25** | 9.7e-9 | 1.01 | 1.9e-8 | 21.12 | 1.0e-6 | 294.52 |
| E2006.train | 1 | 3.9e-9 | **0.44** | 4.4e-10 | 0.87 | 1.4e-11 | 1.13 | 4.4e-5 | 1149.22 |

*The results on Lasso problem (λ = 10^{-3}‖B^T b‖∞).*

| id | nnz | SSNCVX | | SSNAL | | SLEP | | ADMM | |
|---|---|---|---|---|---|---|---|---|---|
| | | η | time | η | time | η | time | η | time |
| uci_CT | 44 | 2.6e-7 | **1.26** | 2.9e-12 | 1.75 | 1.8e-1 | 41.63 | 2.0e-3 | 49.88 |
| log1p.E2006.train | 599 | 3.0e-7 | **33.92** | 5.9e-11 | 68.83 | 3.3e-2 | 1835.32 | 1.2e-1 | 3608.17 |
| E2006.test | 1 | 2.6e-14 | **0.20** | 3.7e-9 | 0.29 | 2.4e-12 | 0.38 | 9.0e-7 | 268.11 |
| log1p.E2006.test | 1081 | 8.8e-9 | **13.72** | 2.7e-10 | 30.1 | 7.5e-2 | 455.56 | 1.6e-1 | 3606.60 |
| pyrim5 | 78 | 5.6e-16 | **2.01** | 5.0e-7 | 2.59 | 1.1e-2 | 108.93 | 3.1e-3 | 3601.09 |
| triazines4 | 260 | 9.5e-16 | **18.48** | 8.3e-8 | 34.44 | 9.2e-2 | 187.45 | 1.2e-2 | 3604.48 |
| abalone7 | 59 | 6.1e-12 | **1.63** | 1.2e-8 | 2.00 | 1.5e-2 | 43.91 | 1.0e-6 | 356.34 |
| bodyfat7 | 3 | 1.0e-16 | **1.14** | 9.7e-8 | 1.51 | 6.1e-4 | 41.98 | 1.3e-4 | 3601.89 |
| housing7 | 281 | 2.6e-11 | **2.51** | 1.2e-7 | 2.52 | 4.1e-2 | 52.60 | 3.6e-4 | 3601.09 |
| mpg7 | 128 | 1.8e-15 | **0.11** | 6.9e-8 | 0.18 | 5.8e-4 | 0.76 | 9.9e-7 | 11.67 |
| spacega9 | 38 | 3.1e-12 | **0.53** | 3.5e-7 | 0.72 | 9.0e-5 | 22.96 | 1.0e-6 | 53.23 |
| E2006.train | 1 | 5.6e-9 | **0.75** | 4.4e-9 | 0.88 | 1.0e-11 | 1.39 | 4.4e-5 | 1132.34 |

*The results of tested algorithms on Lasso problem (λ = 10^{-4}‖B^T b‖∞).*

--------------------------
**Fused Lasso**

| id | nnz(x) | nnz(Bx) | SSNCVX | | SSNAL | | SLEP | | ADMM | |
|---|---|---|---|---|---|---|---|---|---|---|
| | | | η | time | η | time | η | time | η | time |
| uci_CT | 8 | 1 | 6.3e-7 | **0.25** | 7.9e-7 | 0.42 | 1.8e-6 | 2.06 | 7.7e-3 | 41.75 |
| log1p.E2006.train | 31 | 2 | 2.8e-7 | **10.43** | 2.4e-7 | 14.02 | 1.2e-2 | 4889.15 | 1.2e-1 | 3623.18 |
| E2006.test | 1 | 1 | 1.5e-7 | **0.17** | 5.1e-7 | 0.33 | 4.8e-8 | 0.93 | 8.2e-7 | 1768.26 |
| log1p.E2006.test | 33 | 1 | 4.1e-7 | **2.60** | 8.1e-7 | 2.74 | 1.2e-2 | 1690.60 | 2.4e-2 | 3601.25 |
| pyrim5 | 1135 | 74 | 9.1e-7 | **2.34** | 4.5e-7 | 3.40 | 3.4e-2 | 238.43 | 2.4e-3 | 3601.20 |
| triazines4 | 2666 | 206 | 2.1e-7 | **10.24** | 9.8e-7 | 15.49 | 7.8e-2 | 585.70 | 2.8e-2 | 3601.89 |
| bodyfat7 | 63 | 8 | 3.0e-7 | **0.72** | 7.2e-9 | 1.35 | 9.9e-7 | 41.13 | 3.5e-3 | 3612.99 |
| abalone7 | 1 | 1 | 1.6e-7 | **0.83** | 5.3e-8 | 0.95 | 1.3e-3 | 32.51 | 6.4e-4 | 538.90 |
| housing7 | 205 | 47 | 7.6e-7 | **1.98** | 8.2e-7 | 2.73 | 5.0e-3 | 117.07 | 2.2e-2 | 3600.28 |
| mpg7 | 42 | 20 | 1.9e-7 | **0.08** | 1.8e-7 | 0.11 | 3.4e-6 | 3.19 | 6.3e-6 | 156.31 |
| spacega9 | 24 | 11 | 5.0e-8 | **0.27** | 1.2e-7 | 0.44 | 6.1e-8 | 5.32 | 9.9e-7 | 337.14 |
| E2006.train | 1 | 1 | 3.7e-7 | **0.42** | 4.0e-8 | 0.98 | 9.7e-12 | 0.39 | 4.3e-5 | 1196.42 |

*The results of tested algorithms on Fused Lasso problem (λ₁ = 10^{-3} ‖B*b‖∞ and λ₂ = 5λ₁).*

| id | nnz(x) | nnz(Bx) | SSNCVX | | SSNAL | | SLEP | | ADMM | |
|---|---|---|---|---|---|---|---|---|---|---|
| | | | η | time | η | time | η | time | η | time |
| uci_CT | 18 | 8 | 6.3e-7 | **0.40** | 8.9e-10 | 0.42 | 1.8e-6 | 2.06 | 7.7e-3 | 39.29 |
| log1p.E2006.train | 8 | 3 | 7.0e-7 | **8.37** | 1.5e-7 | 12.6 | 1.2e-2 | 4889.15 | 1.2e-1 | 3606.14 |
| E2006.test | 1 | 1 | 1.5e-7 | **0.17** | 2.9e-8 | 0.33 | 4.8e-8 | 0.93 | 7.7e-7 | 699.27 |
| log1p.E2006.test | 32 | 5 | 3.1e-9 | **3.07** | 1.2e-8 | 3.31 | 1.2e-2 | 1690.60 | 7.9e-2 | 3601.20 |
| pyrim5 | 327 | 97 | 9.1e-7 | **2.34** | 2.0e-7 | 3.06 | 3.4e-2 | 238.43 | 1.5e-3 | 3601.13 |
| triazines4 | 1244 | 286 | 8.2e-7 | **10.51** | 2.4e-7 | 12.63 | 7.8e-2 | 585.70 | 2.8e-2 | 3603.56 |
| bodyfat7 | 2 | 3 | 2.8e-8 | **0.81** | 4.7e-8 | 0.89 | 9.9e-7 | 41.13 | 2.7e-3 | 3606.85 |
| abalone7 | 26 | 15 | 3.7e-7 | **0.49** | 5.0e-9 | 1.17 | 1.3e-3 | 32.51 | 5.0e-4 | 545.23 |
| housing7 | 131 | 117 | 6.4e-7 | **1.46** | 3.9e-7 | 2.4 | 5.0e-3 | 117.07 | 2.0e-2 | 3603.08 |
| mpg7 | 32 | 39 | 6.7e-7 | **0.07** | 2.2e-7 | 0.15 | 3.4e-6 | 3.19 | 1.0e-6 | 77.58 |
| spacega9 | 14 | 13 | 8.7e-7 | **0.22** | 1.7e-7 | 0.44 | 6.1e-8 | 5.32 | 1.0e-6 | 333.39 |
| E2006.train | 1 | 1 | 4.2e-7 | **0.45** | 4.0e-7 | 1.12 | 9.7e-12 | 0.39 | 4.4e-5 | 1189.36 |

*The results of tested algorithms on Fused Lasso problem (λ₁ = 10^{-3} ‖B*b‖∞ and λ₂ = λ₁).*

--------------------------
**QP**

For synthetic data, `Q` and `c` are generated as follows:
```MATLAB
p = 0.01*n;
F = sprandn(n, p, 0.1); D = sparse(diag(sqrt(p)*rand(n,1)));
Q = cov(F') + D;
c = randn(n,1);
```

| problem | SSNCVX | | HIGHS | |
|---|---|---|---|---|
| | η | time | η | time |
| Aug2D | 2.9e-11 | **0.25** | - | - |
| Aug2DC | 7.5e-13 | **0.20** | - | - |
| Aug2DCQP | 7.5e-13 | **0.18** | - | - |
| Aug2DQP | 1.7e-16 | **0.31** | - | - |
| BOYD1 | 2.0e-7 | **47.80** | - | - |
| BOYD2 | 2.3e-9 | **0.29** | 4.3e-6 | 3667.91 |
| CONT-100 | 7.0e-12 | **1.23** | 7.8e-4 | 122.04 |
| CONT-101 | 0.0e+0 | **0.07** | 4.5e-3 | 3600.04 |
| CONT-200 | 4.3e-8 | **3.96** | 3.2e-3 | 3600.09 |
| CONT-201 | 0.0e+0 | **0.16** | - | - |
| CONT-300 | 0.0e+0 | **0.24** | 0.0e+0 | 4011.90 |
| DTOC-3 | 8.9e-18 | **0.39** | - | - |
| LISWET1 | 2.3e-18 | **0.15** | 6.8e-6 | 0.70 |
| UBH1 | 4.8e-9 | **0.28** | - | - |
| random512_1 | 7.2e-11 | **0.36** | 2.1e-7 | 1.10 |
| random512_2 | 7.4e-13 | **0.40** | 2.5e-7 | 1.11 |
| random1024_1 | 2.2e-9 | **1.41** | 4.0e-7 | 2.32 |
| random1024_2 | 2.7e-8 | **0.81** | 2.5e-7 | 2.32 |
| random2048_1 | 1.7e-7 | **3.40** | 2.5e-7 | 3.96 |
| random2048_2 | 2.6e-10 | **2.92** | 1.4e-7 | 4.06 |

*Computational results of tested algorithms on portfolio optimization.*

--------------------------
**SOCP**

Note that the MATLAB solvers (`SSNCVX` and `SDPT3`) solve the preprocessed datasets with preprocessing time excluded. This preprocessing, which typically requires several seconds, significantly reduced solution times for some instances (e.g., `firL2a`), making these solvers appear faster for such problems. However, as geometric means are calculated with a 10-second shift, the exclusion has a negligible impact on the overall results.

| id | SSNCVX | | SDPT3 | | ECOS | | MOSEK | |
|---|---|---|---|---|---|---|---|---|
| | η | time | η | time | η | time | η | time |
| beam7 | - | - | - | - | 1.0e-7 | 206.0 | 6.0e-4 | 19.7 |
| beam30 | - | - | - | - | 3.0e-7 | 2464.7 | 3.0e-6 | 96.5 |
| chainsing-50000-1 | 1.5e-7 | 5.8 | 6.9e-7 | 5.5 | - | - | 1.6e-6 | 3.8 |
| chainsing-50000-2 | 7.3e-7 | 14.4 | 7.0e-7 | 9.5 | - | - | 1.0e-7 | 4.1 |
| chainsing-50000-3 | 5.0e-9 | 15.7 | 1.4e-7 | 19.4 | - | - | 1.0e-8 | 2.0 |
| db-joint-soerensen | - | - | - | - | - | - | 2.0e-8 | 36.3 |
| db-plate-yield-line | 8.5e-7 | 597.2 | 8.7e-7 | 217.6 | - | - | 5.0e-7 | 6.2 |
| dsNRL | 1.0e-6 | 859.2 | 8.9e-7 | 567.8 | - | - | 8.2e-10 | 67.1 |
| firL1 | 5.3e-11 | 101.6 | 7.8e-7 | 582.0 | 3.0e-8 | 1305.2 | 3.1e-9 | 20.5 |
| firL1Linfalph | 8.4e-7 | 509.6 | 7.5e-7 | 916.2 | 3.0e-8 | 2846.6 | 4.0e-9 | 91.8 |
| firL1Linfeps | 7.0e-7 | 86.4 | 8.2e-7 | 179.1 | 2.0e-9 | 2530.8 | 3.0e-8 | 27.5 |
| firL2a | 1.4e-8 | 0.4 | 6.1e-7 | 0.1 | 2.0e-9 | 944.6 | 2.0e-13 | 4.4 |
| firL2L1alph | 1.1e-7 | 37.4 | 7.3e-7 | 131.7 | 3.0e-9 | 201.5 | 2.2e-10 | 5.8 |
| firL2L1eps | 2.0e-9 | 159.5 | 6.2e-7 | 586.0 | 2.0e-8 | 796.6 | 3.5e-9 | 17.2 |
| firL2Linfalph | 7.9e-7 | 89.1 | 7.9e-7 | 799.9 | - | - | 9.0e-9 | 41.7 |
| firL2Linfeps | 5.2e-7 | 72.4 | 8.0e-9 | 251.2 | 5.0e-10 | 687.1 | 1.0e-8 | 29.9 |
| firLinf | 1.4e-7 | 280.2 | 7.1e-7 | 576.7 | 5.0e-9 | 3478.7 | 1.0e-8 | 123.6 |
| wbNRL | 8.7e-7 | 20.1 | 5.9e-7 | 151.2 | 5.0e-9 | 1332.6 | 2.4e-9 | 11.8 |
| geomean | - | 155.0 | - | 267.8 | - | 1731.4 | - | 22.7 |

*The results on Hans Mittelmann's SOCP benchmark.*

--------------------------
**SPCA**

For random examples, `L`, the quadratic coefficient matrix, is generated by: `L = (1/||u||_2) * u * u^T + V * V^T`, where `u = [1, 1/2, ..., 1/n]` and each entry of `V ∈ R^{n×n}` is randomly uniformly chosen from `[0,1]`.

| problem | SSNCVX |  | superSCS |  |
|---|---|---|---|---|
|  | η | time | η_K | time |
| 20news | 2.0e-12 | 0.8 | 1.0e-6 | 9.6 |
| bibtex | 1.2e-11 | 76.6 | 2.7e-1 | 3626.4 |
| colon_cancer | 5.5e-12 | 45.9 | 4.9e-1 | 3647.9 |
| delicious | 2.6e-12 | 2.9 | 2.5e-3 | 2813.5 |
| dna | 1.2e-13 | 0.3 | 1.0e-6 | 29.2 |
| gisette | 2.5e-12 | 1190.0 | 7.0e-1 | 3703.5 |
| madelon | 5.9e-15 | 16.7 | 4.4e-5 | 3343.6 |
| mnist | 4.0e-17 | 15.7 | 1.0e-6 | 195.4 |
| protein | 3.5e-11 | 3.7 | 8.7e-3 | 2334.1 |
| random1024_1 | 9.3e-18 | 2.8 | 3.2e-2 | 3603.3 |
| random1024_2 | 4.4e-18 | 2.7 | 1.9e-3 | 3604.8 |
| random1024_3 | 1.3e-17 | 2.8 | 1.4e-3 | 3608.3 |
| random2048_1 | 7.8e-18 | 3.3 | 2.3e-1 | 3605.5 |
| random2048_2 | 5.1e-18 | 3.5 | 5.9e-2 | 3607.0 |
| random2048_3 | 1.5e-18 | 2.3 | 1.5e-2 | 3608.2 |
| random4096_1 | 8.2e-18 | 73.4 | N/A | 3655.4 |
| random4096_2 | 3.5e-18 | 73.1 | 1.2e-2 | 3638.0 |
| random4096_3 | 6.7e-19 | 72.4 | 9.6e-3 | 3645.0 |
| random512_1 | 4.3e-18 | 0.6 | 1.0e-6 | 252.0 |
| random512_2 | 1.1e-17 | 0.6 | 8.1e-3 | 2938.5 |
| random512_3 | 5.7e-18 | 0.6 | 8.2e-3 | 2802.0 |
| usps | 2.4e-13 | 1.1 | 1.0e-6 | 229.8 |

*Computational results of SSNCVX and superSCS on SPCA.*

--------------------------
**LRMC**

| Problem | SSNCVX | | ADMM | | PG | | APG | |
|---|---|---|---|---|---|---|---|---|
| | η | Time | η | Time | η | Time | η | Time |
| Image1 | 1.5e-9 | **20.1** | 9.9e-9 | 84.5 | 9.6e-9 | 122.4 | 9.9e-9 | 55.8 |
| Image2 | 4.3e-9 | **22.1** | 1.0e-8 | 84.0 | 9.8e-9 | 120.9 | 9.6e-9 | 54.5 |
| Image3 | 5.3e-9 | **23.2** | 9.9e-9 | 82.8 | 9.6e-9 | 119.5 | 9.3e-9 | 53.9 |
| Image4 | 3.3e-9 | **25.3** | 9.7e-9 | 84.1 | 9.8e-9 | 121.1 | 9.8e-9 | 54.6 |
| Image5 | 7.4e-9 | **20.3** | 9.5e-9 | 83.7 | 9.7e-9 | 120.4 | 9.9e-9 | 54.4 |
| Image6 | 1.9e-9 | **20.9** | 1.0e-8 | 83.5 | 9.8e-9 | 120.4 | 9.7e-9 | 54.3 |
| Image7 | 1.6e-9 | **20.2** | 9.9e-9 | 82.2 | 9.9e-9 | 118.3 | 9.7e-9 | 53.1 |
| Image8 | 2.3e-9 | **20.8** | 9.8e-9 | 83.0 | 9.7e-9 | 120.0 | 9.7e-9 | 53.9 |

*Comparison of tested algorithms on the LRMC problem.*

## References

This software is released under the [MIT License](https://opensource.org/licenses/MIT).

If you use SSNCVX in your research, please cite [this paper](#).

```bibtex
@article{deng2025augmented,
}
```

If you want to contact us for collaboration or other inquiries, please reach out via email at [to_be_determined@pku.edu.cn]


