function [By1,By2] = matvec_yz(K, At, par, y1, y2, AL)

if (nargin < 5); AL = []; end
if isempty(AL); existAL = 0; else; existAL = 1; end
N = length(y1);
% if (norm(y1) == 0); By1 = zeros(N,1); return; end
%%
yorg1 = y1;
yorg2 = y2;
if ~AL.isidentity  && (existAL)
    if strcmp(AL.matfct_options,'chol')
        y1(AL.perm) = AL.R \ y1;
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        y1(AL.perm) = mexbwsolve(AL.Rt,y1);
    end
end

%% perform A * (par.sig * par.Dsch + par.epsilon * identity) * A' * y
By1 = zeros(N,1);
By2 = 0; %modified
for p = 1: length(K)
    cone = K{p};
    n = sum(cone.size);
    if strcmp(cone.type,'s')
        rr = size(par.P1{p}, 2);
        Aty = Atyfun(cone, At{p}, y1)+y2;
        if (rr > 0 && rr < n)
            if (rr <= n/2)
                tmp0 = par.P1{p}'*Aty;
               tmp1 = (tmp0*par.P1{p})*par.P1{p}';         
               tmp2 = par.Dsch2{p}.*(tmp0*par.P2{p});
               tmp2 = tmp2*par.P2t{p}; 
               tmp3 = par.P1{p}*(0.5*par.Dsch1{p}*tmp1 + tmp2);
                By2 =  tmp3+tmp3';
                By1 = By1 + par.sig*AXfun(cone, At{p}, By2);
            else
               tmp0 = par.P2t{p}*Aty;
               tmp1 = (tmp0*par.P2{p})*par.P2t{p};         
               tmp2 = (par.Dsch1{p}-par.Dsch2{p}').*(tmp0*par.P1{p});
               tmp2 = tmp2*par.P1t{p}; 
               tmp3 = par.P2{p}*(0.5*par.Dsch1{p}*tmp1 + tmp2);
                By2 =  par.Dsch1{p}*Aty-tmp3-tmp3';
                By1 = By1 + par.sig*AXfun(cone, At{p}, By2);
            end
        elseif (rr == n)
            By1 = By1 + par.sig*AXfun(cone, At{p}, Aty);
        end
    elseif strcmp(cone.type, 'q')
        if (~isempty(par.Dsch2{p}))
            Aty = At{p} * y1;
            P1blkdiag = blk_spdiag(par.P1{p}, cone.size);
            P2blkdiag = blk_spdiag(par.P2{p}, cone.size);
            tmp = repelem(par.shift{p}, cone.size, 1) .* Aty + ...
                P1blkdiag * spdiag(par.Dsch1{p}) * (P1blkdiag' * Aty)+ ...
                P2blkdiag * spdiag(par.Dsch2{p}) * (P2blkdiag' * Aty);
            By1 = By1 + par.sig*(tmp'*At{p})';
        end
    elseif strcmp(cone.type, 'l')
        if (~isempty(par.Dsch2{p}))
            tmp = par.Dsch2{p}.*(At{p}*y1);
            By1 = By1 + par.sig*(tmp'*At{p})';
        end
    elseif strcmp(cone.type,'u')
        tmp = At{p}*y1;
        By1 = By1 + par.sig*(tmp'*At{p})';
    end
end
if (existAL) && ~AL.isidentity 

    if strcmp(AL.matfct_options,'chol')
        By1 = AL.Rt \ By1(AL.perm);
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        By1 = mexfwsolve(AL.R,By1(AL.perm));
    end
end




By1 = By1 + par.epsilon1*yorg1;

By2 = By2 +  yorg2 .* (par.Dh4 + par.epsilon2) ;

%%************************************************************************
