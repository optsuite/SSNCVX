%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-11 21:45:36
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%%  In many algorithms for conic programming, the Jacobian of some operator (e.g., projection, smoothing projection) have a special structure.  We discuss the cases of different cones.
% Denote the Jacobian at point Z as J(Z).
% semidenite cone: 
%   Z \in S^n, J(Z) is a operator from S^n to S^n with
%   J(Z)[B] = shift * B + U * ( Sigma .* (U' * B * U)) * U', B \in S^n
%   where Z = U * Lambda * U' is the spectral decomposition of Z.
%   and Sigma is a matrix of size [n, n] with structure
%   Sigma = [Sigma11, Sigma12;
%           Sigma12', 0];
% quadratic cone:
%   Z \in R^n, J(Z) is a matrix of size [n, n] with
%   J(Z) = shift * I + Sigma1 * u1 * u1' + Sigma2 * u2 * u2';
%   where Z = lambda1 * u1 + lambda2 * u2 is the spectral decomposition of Z. 
%   shift, Sigma1, Sigma2 are scalars, u1, u2 are vectors of length n.
%   For quadratic cone, u1, u2 have explicit expressions as follows:
%   u1 = 1 / 2 * [1; Z(2:n) / norm(Z(2:n))];
%   u2 = 1 / 2 * [1; - Z(2:n) / norm(Z(2:n))];
% nonnegative orthant / Euclead space
%   Z \in R^n, J(Z) is a digonal matrix of size [n, n] with
%   J(Z) = diag(Sigma)
%  where Sigma is a vector of length n
%
% Hence we consider representing this structure with a class ProjJac.
% It has fields: dd, Dsch1, Dsch2, P1, P2, posidx. Each of them is a cell.
% Their meaning is as follows:
% Assume cone.size = [n1, n2, ..., nk], where n1 + n2 + ... + nk = n
% semidenite cone ( cone.type = 's' ): 
%    dd: eigenval of Z, i.e., dd = diag(Lambda), size(dd) = [n, 1]
%    Dsch1: Sigma11, sometimes Sigma11 is a scalar multiple of identity matrix, then Dsch1 is a scalar or a vector of length k
%    Dsch2: Sigma12, size(Dsch2) = [length(posidx), n - length(posidx)]
%    posidx: index of positive eigenvalues of Z
%    P1: eigenvec of positive eigenvalues of Z, i.e., P1 = U(:, posidx)
%    P2: eigenvec of nonpoitive eigenvalues of Z, i.e., P2 = U(:, ~posidx)
%    shift: shift, scalar
% quadratic cone ( cone.type = 'q' ):
%    dd: eigenvalues of Z, i.e., dd(i, :) = [zi(1) + norm(z(2:ni)), zi(1) - norm(z(2:ni))], size(dd) = [k, 2]
%    Dsch1: Sigma1, size(Dsch1) = [k, 1]
%    Dsch2: Sigma2, size(Dsch2) = [k, 1]
%    posidx: index of positive eigenvalues of Z, i.e., posidx = find(dd > 0)
%    P1: concatenation of u1, size(P1) = [n, 1]
%    P2: concatenation of u2, size(P2) = [n, 1]
%    shift: shift, size(shift) = [k, 1]
% nonnegative orthant / Euclead space ( cone.type = 'l' or 'u' ):
%    dd: []
%    Dsch1: []
%    Dsch2: Sigma, size(Dsch2) = [n, 1]
%    posidx: index of positive elements of Z
%    P1: []
%    P2: []
%    shift: []

%%  In many algorithms for conic programming, 
% the lhs of newton system is constrcuted as 
%    sigma * A * J * A' + epsilon * I 
% hence we add two fields sig and epsilon


% Besides, we add a field precond to indicate whether the preconditioning is needed.

% Now we implement the ProjJac as a struct rather than a class
% because users can add custom fields to it more conveniently.


function obj = ProjJac(num_cone)
    obj = struct();
    obj.dd = cell(num_cone, 1);
    obj.Dsch1 = cell(num_cone, 1);
    obj.Dsch2 = cell(num_cone, 1);
    obj.P1 = cell(num_cone, 1);
    obj.P2 = cell(num_cone, 1);
    obj.posidx = cell(num_cone, 1);
    obj.shift = cell(num_cone, 1);
    obj.precond = 0;
end

