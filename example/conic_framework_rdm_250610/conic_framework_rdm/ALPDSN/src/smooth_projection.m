%%************************************************************************
%% smooth_project: compute smoothing projection onto the cone of positive
%%           semidefinite matrices.
%%
%% [Xp,par] = smooth_projection(blk,X,par, smoothfun, smooth_par);
%%
%%************************************************************************

function [Xp, par, d2Xp] = smooth_projection(blk, X, smoothfun, d1smoothfun, d2smoothfun, smooth_par)

tol = 1e-15;  %% must be small as it will affect gap
addtol = 1e-6;
Xp = cell(size(X));
d2Xp = cell(size(X));
par.P1 = cell(size(X));     par.P2 = cell(size(X));
par.Dsch2 = cell(size(X)); par.dd = cell(size(X));
par.posidx = cell(size(X));
%%
for p = 1:size(blk,1)
   pblk = blk(p,:);
   if strcmp(pblk{1},'s')
      blktmp = pblk{1,2};
      if (length(blktmp) > 1) || (blktmp < 100) % only one subblock
         n = sum(blktmp);
         if smooth_par > 1e-16
            [Xptmp,P,Dsch2,dd,posidx, Dsch11, d2Xptmp] = smooth_projectfun(X{p},tol, smoothfun, d1smoothfun, d2smoothfun, smooth_par);
         else
            [Xptmp,P,Dsch2,dd,posidx] = projectfun(X{p},tol);
            Dsch11 = [];
            d2Xptmp = [];
         end
         %%***** perturbation *****
         Dsch2 = max(addtol,Dsch2);
         Dsch11 = max(addtol,Dsch11);
         if ~isempty(posidx)
            P1 = P(:,posidx);
            P2 = P(:,setdiff(1:n,posidx));
         else
            P1 = [];
            P2 = P;
         end
      else   % more than one subblocks or one very large subblock
         n = sum(blktmp);
         if smooth_par > 1e-16
            [Xptmp,P,Dsch2,dd,posidx, Dsch11, d2Xptmp] = smooth_projectfun_sparse(X{p},tol, smoothfun, d1smoothfun, d2smoothfun, smooth_par);
         else
            [Xptmp,P,Dsch2,dd,posidx] = projectfun_sparse(X{p},tol);
            Dsch11 = [];
            d2Xptmp = [];
         end
         %%***** perturbation *****
         Dsch2 = max(addtol,Dsch2);
         Dsch11 = max(addtol,Dsch11);
         if ~isempty(posidx)
            P1 = P(:,posidx);
            P2 = P(:,setdiff(1:n,posidx));
         else
            P1 = [];
            P2 = P;
         end
      end
   elseif strcmp(pblk{1},'l')
      P1 = []; P2 = [];
      n  = sum(pblk{1,2});
      Xptmp = zeros(n,1);
      d2Xptmp = zeros(n,1);
      %%***** perturbation *****
      %Dsch2 = addtol*ones(n,1);
      Dsch2 = zeros(n,1);
      Dsch11 = [];
      posidx = find(X{p} > tol);
      if ~isempty(posidx)
         Xptmp(posidx)  = abs(X{p}(posidx));
         d2Xptmp(posidx) = d2smoothfun(X{p}(posidx), smooth_par);
         Dsch2(posidx) = ones(length(posidx),1);
      end
      dd = X{p};
   elseif strcmp(pblk{1},'u')
      P1 = []; P2 = [];
      Xptmp = X{p};
      d2Xptmp = d2smoothfun(X{p}, smooth_par);
      %%***** perturbation *****
      Dsch2 = [];
      Dsch11 = [];
      posidx = [];
      dd = [];
   end
   Xp{p}     = Xptmp;
   d2Xp{p}   = d2Xptmp;
   par.P1{p} = P1;
   par.P2{p} = P2;
   par.dd{p} = dd;
   par.posidx{p} = posidx;
   par.Dsch2{p} = Dsch2;
   par.Dsch11{p} = Dsch11;
end
end
%%***************************************************************************

function [Xp,V,Dsch2,d,posidx, Dsch11, d2Xp] = smooth_projectfun(X, tol, smoothfun, d1smoothfun, d2smoothfun, smooth_par)
% [Xp,V,Dsch2,d,posidx, Dsch11] = smooth_projectfun(X, tol) % euiqvalent to projectfun(X,tol)

% smoothfun: a mapping g(lambda, smooth_par) from R to R or element mapping from R^n to R^n
%              its function value is nonnegative and nondecreasing in lambda
% d1smoothfun: the derivative of smoothfun w.r.t. lambda
% d2smoothfun: the derivative of smoothfun w.r.t. mu
% smooth_par: nonnegative scalar


assert(smooth_par >= 0);
if  nargin == 2 || smooth_par < 1e-16
   [Xp,V,Dsch2,d,posidx] = projectfun(X,tol);
   if nargout == 6
      Dsch11 = [];
   end
   return;
end

n = length(X);
exist_mexeig = exist('mexeig','file');
X(abs(X)<1e-14) = 0;
X = 0.5*(X+X');
if (exist_mexeig==3)
   %[V,D] = mexeig(full(X));
   [V,D] = eig(full(X));
else
   [V,D] = eig(full(X));
end
d = diag(D);
[d,idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:,idx);
posidx = find(d > tol);
if isempty(posidx)
   Xp = sparse(n,n);
   Dsch2 = [];
   Dsch11 = [];
   d2Xp = sparse(n,n);
else
   r = length(posidx); s = n-r;
   negidx = [r+1:n];
   dp = abs(d(posidx));
   dn = abs(d(negidx));
   hdp = smoothfun(dp, smooth_par);  % the smooth function is denoted by h in the paper hence the name hdp
   d2hdp = d2smoothfun(dp, smooth_par); % the derivative of h with respect to mu
   Vtmp = V(:,posidx)*diag(sqrt(hdp));
   Xp = Vtmp*Vtmp';
   Xp = 0.5*(Xp+Xp');
   Vtmp2 = V(:,posidx)*diag(sqrt(d2hdp));
   d2Xp = Vtmp2*Vtmp2';
   d2Xp = 0.5*(d2Xp+d2Xp');
   
   Dsch2 = (hdp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
   d1hdp = zeros(size(dp));
   idxeq = find_equal_values(dp);
   d1hdp(idxeq) = d1smoothfun(dp(idxeq), smooth_par); % the derivative of hdp, only need to compute at the equal values
   Dsch11 = compute_Dsch11(dp, hdp, d1hdp);
end
end


function [Xp,V,Dsch2,d,posidx] = projectfun(X,tol)
n = length(X);
exist_mexeig = exist('mexeig','file');
X(abs(X)<1e-14) = 0;
X = 0.5*(X+X');
if (exist_mexeig==3)
   %[V,D] = mexeig(full(X));
   [V,D] = eig(full(X));
else
   [V,D] = eig(full(X));
end
d = diag(D);
[d,idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:,idx);
posidx = find(d > tol);
if isempty(posidx)
   Xp = sparse(n,n);
   Dsch2 = [];
elseif (length(posidx) == n)
   Xp = X;
   Dsch2 = [];
else
   r = length(posidx); s = n-r;
   negidx = [r+1:n];
   dp = abs(d(posidx));
   dn = abs(d(negidx));
   Vtmp = V(:,posidx)*diag(sqrt(dp));
   Xp = Vtmp*Vtmp';
   Xp = 0.5*(Xp+Xp');
   Dsch2 = (dp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
end
end

function [Xp,V,Dsch2,d,posidx, Dsch11, d2Xp] = smooth_projectfun_sparse(X, tol, smoothfun, d1smoothfun, d2smoothfun, smooth_par)

% smoothfun: a mapping g(lambda, smooth_par) from R to R or element mapping from R^n to R^n
%              its function value is nonnegative and nondecreasing in lambda
% dsmoothfun: the derivative of smoothfun
% smooth_par: nonnegative scalar


n = length(X);
exist_mexeig = exist('mexeig','file');
X(abs(X)<1e-14) = 0;

X = 0.5*(X+X');
pm = symamd(X);
Xperm = X(pm,pm);
[t,~] = etree(Xperm);
idx0 = find(t == 0); % split the matrix into blocks
len0 = length(idx0);
Dsub = zeros(n,1);
d = zeros(n,1);
Vsub = cell(n,1);
offset = 1;
for it = 1:len0
   idx = offset:idx0(it);
   if (exist_mexeig==3)
      %[Vsub0,Dsub0] = mexeig(full(Xperm(idx,idx)));
      [Vsub0,Dsub0] = eig(full(Xperm(idx,idx)));
   else
      [Vsub0,Dsub0] = eig(full(Xperm(idx,idx)));
   end
   Dsub0 = diag(Dsub0);
   Vsub{it} = Vsub0;
   Dsub(idx) = Dsub0;
   offset = idx0(it)+1;
end
V = blkdiag(Vsub{:});
V(pm,pm) = V;
d(pm) = Dsub;

[d,idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:,idx);
posidx = find(d > tol);
if isempty(posidx)
   Xp = sparse(n,n);
   Dsch2 = [];
   Dsch11 = [];
   d2Xp = sparse(n,n);
else
   r = length(posidx); s = n-r;
   negidx = [r+1:n];
   dp = abs(d(posidx));
   dn = abs(d(negidx));
   hdp = smoothfun(dp, smooth_par);  % the smooth function is denoted by h in the paper hence the name hdp
   d2hdp = d2smoothfun(dp, smooth_par); % the derivative of h with respect to mu
   Vtmp = V(:,posidx)*diag(sqrt(hdp));
   Xp = Vtmp*Vtmp';
   Xp = 0.5*(Xp+Xp');
   Vtmp2 = V(:,posidx)*diag(sqrt(d2hdp));
   d2Xp = Vtmp2*Vtmp2';
   d2Xp = 0.5*(d2Xp+d2Xp');
   
   Dsch2 = (hdp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
   d1hdp = zeros(size(dp));
   idxeq = find_equal_values(dp);
   d1hdp(idxeq) = d1smoothfun(dp(idxeq), smooth_par); % the derivative of hdp, only need to compute at the equal values
   Dsch11 = compute_Dsch11(dp, hdp, d1hdp);
end
end
%%************************************************************************


function [Xp,V,Dsch2,d,posidx] = projectfun_sparse(X,tol)

n = length(X);
exist_mexeig = exist('mexeig','file');
X(abs(X)<1e-14) = 0;

X = 0.5*(X+X');
pm = symamd(X);
Xperm = X(pm,pm);
[t,~] = etree(Xperm);
idx0 = find(t == 0);
len0 = length(idx0);
Dsub = zeros(n,1);
d = zeros(n,1);
Vsub = cell(n,1);
offset = 1;
for it = 1:len0
   idx = offset:idx0(it);
   if (exist_mexeig==3)
      %[Vsub0,Dsub0] = mexeig(full(Xperm(idx,idx)));
      [Vsub0,Dsub0] = eig(full(Xperm(idx,idx)));
   else
      [Vsub0,Dsub0] = eig(full(Xperm(idx,idx)));
   end
   Dsub0 = diag(Dsub0);
   Vsub{it} = Vsub0;
   Dsub(idx) = Dsub0;
   offset = idx0(it)+1;
end
V = blkdiag(Vsub{:});
V(pm,pm) = V;
d(pm) = Dsub;

[d,idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:,idx);
posidx = find(d > tol);
if isempty(posidx)
   Xp = sparse(n,n);
   Dsch2 = [];
elseif (length(posidx) == n)
   Xp = X;
   Dsch2 = [];
else
   r = length(posidx); s = n-r;
   negidx = [r+1:n];
   dp = abs(d(posidx));
   dn = abs(d(negidx));
   Vtmp = V(:,posidx)*diag(sqrt(dp));
   Xp = Vtmp*Vtmp';
   Xp = 0.5*(Xp+Xp');
   Dsch2 = (dp*ones(1,s))./(dp*ones(1,s) + ones(r,1)*dn');
end
end
%%************************************************************************
