


function By = matvec_y2(K,At,par,y,AL)

if (nargin < 5); AL = []; end;
if isempty(AL); existAL = 0; else; existAL = 1; end


N = length(y);
if (norm(y) == 0); By = zeros(N,1); return; end
%%
yorg = y;
if ~AL.isidentity && (existAL)
    if strcmp(AL.matfct_options,'chol')
        y(AL.perm) = AL.R \ y;
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        y(AL.perm) = mexbwsolve(AL.Rt,y);
    end
end
%%
By = zeros(N,1);
for p = 1:length(K)
    cone = K{p};

    n = sum(cone.size);
    if strcmp(cone.type,'s')
        rr = size(par.P1{p},2);
        %          Aty = Atyfun_sdpnal(pblk,At{p},y);
        Aty = Atyfun(cone,At{p},y);
        if (rr > 0 && rr < n)
            if (rr <= n/2)
                tmp0 = par.P1{p}'*Aty;
                tmp1 = (tmp0*par.P1{p})*par.P1{p}';
                tmp2 = par.Dsch12{p}.*(tmp0*par.P2{p});
                tmp2 = tmp2*par.P2t{p};
                tmp3 = par.P1{p}*(0.5*par.Dsch1{p}*tmp1 + tmp2);
                %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},tmp3+tmp3');
                By = By + par.sig*AXfun(cone,At{p},tmp3+tmp3');
            else
                tmp0 = par.P2{p}'*Aty;
                tmp1 = (tmp0*par.P2{p})*par.P2t{p};
                tmp2 = (par.Dsch1{p}-par.Dsch12{p}').*(tmp0*par.P1{p});
                tmp2 = tmp2*par.P1t{p};
                tmp3 = par.P2{p}*(0.5*par.Dsch1{p}*tmp1 + tmp2);
                %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},par.Dsch1{p}*Aty-tmp3-tmp3');
                By = By + par.sig*AXfun(cone,At{p},par.Dsch1{p}*Aty-tmp3-tmp3');
            end
        elseif (rr == n)
            %             By = By + par.sig*AXfun_sdpnal(pblk,At{p},Aty);
            By = By + par.Dsch1{p}*par.sig*AXfun(cone,At{p},Aty);
        end
    elseif strcmp(cone.type,'l')
        if (~isempty(par.Dsch12{p}))
            tmp = par.Dsch12{p}.*(At{p}*y);
            By = By + par.sig*(tmp'*At{p})';
        end
    elseif strcmp(cone.type,'u')
        tmp = At{p}*y;
        By = By + par.sig*(tmp'*At{p})';
    end
end
if ~AL.isidentity && (existAL)
    if strcmp(AL.matfct_options,'chol')
        By = AL.Rt \ By(AL.perm);
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        By = mexfwsolve(AL.R,By(AL.perm));
    end
end
%
%    if (par.use_proximal);
%       By = By + par.H2.*(y/par.sighat);
%    else
%       sighat = max([1e4,10*par.sig]);
%       By = By + y/sighat;
%    end

%    if isfield(par,'epsilon')
By = By + par.epsilon*yorg;

%    end
%%************************************************************************
