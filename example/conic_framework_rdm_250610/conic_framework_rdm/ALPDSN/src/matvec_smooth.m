%%***************************************NEWTopts.sigPow*********************************
%% matvec: matrix-vector multiply.
%% matrix B = sig* (A* hatDPi*At) = sig * A(PxP) hatDsch (PtxPt)At
%%
%% SDPNAL:
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh
%%************************************************************************

function By = matvec(blk,At,par,y,ecoe,AL)
% ecoe is of no use here
if (nargin < 5); AL = []; end;
if isempty(AL); existAL = 0; else; existAL = 1; end

if ~isfield(par,'P1t'); par.P1t = ops_sdpnal(par.P1,'transpose'); end
if ~isfield(par,'P2t'); par.P2t = ops_sdpnal(par.P2,'transpose'); end
if ~isfield(par,'Dsch21'); par.Dsch21 = ops_sdpnal(par.Dsch2,'transpose'); end

N = length(y);
if (norm(y) == 0); By = zeros(N,1); return; end
%%
yorg = y;
if (existAL)
   if strcmp(AL.matfct_options,'chol')
      y(AL.perm) = AL.R \ y;
   elseif strcmp(AL.matfct_options,'spcholmatlab')
      y(AL.perm) = mexbwsolve_sdpnal(AL.Rt,y);
   end
end
%%
By = zeros(N,1);
for p = 1:size(blk,1)
   pblk = blk(p,:);
   n = sum(pblk{2});
   if strcmp(pblk{1},'s')
      rr = size(par.P1{p},2);
      if rr > 0
         Aty = Atyfun_sdpnal(pblk,At{p},y);   % S
         tmp0 = par.P1{p}'*Aty;  % U = Q_{\alpha}^T S
         tmp1 = (par.Dsch11{p} .* (tmp0*par.P1{p}) ) * par.P1{p}';
         tmp2 = (par.Dsch2{p} .* (tmp0*par.P2{p})) * par.P2t{p} ;
         tmp3 = par.P1{p}*(0.5*tmp1 + tmp2);
         By = By + par.sig*AXfun_sdpnal(pblk,At{p}, tmp3+tmp3' );
      end
   elseif strcmp(pblk{1},'l')
      if (~isempty(par.Dsch2{p}))
         tmp = par.Dsch2{p}.*(At{p}*y);
         By = By + par.sig*(tmp'*At{p})';
      end
   elseif strcmp(pblk{1},'u')
      tmp = At{p}*y;
      By = By + par.sig*(tmp'*At{p})';
   end
end
if (existAL)
   if strcmp(AL.matfct_options,'chol')
      By = AL.Rt \ By(AL.perm);
   elseif strcmp(AL.matfct_options,'spcholmatlab')
      By = mexfwsolve_sdpnal(AL.R,By(AL.perm));
   end
end
%
%    if (par.use_proximal);
%       By = By + par.H2.*(y/par.sighat);
%    else
%       sighat = max([1e4,10*par.sig]);
%       By = By + y/sighat;
%    end

if isfield(par,'epsilon')
   By = By + par.epsilon*yorg;
end
%%************************************************************************
