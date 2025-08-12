


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
sumlen = 0;
By = zeros(N,1);
for p = 1:length(K)
    cone = K{p};
    % sumlen = sumlen + length(cone.size);
    if strcmp(cone.type,'s')
        cumsumlen = [0 cumsum(cone.size)];
        cumsumsquare = (cone.size.^2 + cone.size)/2;
        cumsumsquare = [0 cumsum(cumsumsquare )];
        %          Aty = Atyfun_sdpnal(pblk,At{p},y);
        Aty = Atyfun(cone,At{p},y);
        for  jj = 1:length(cone.size)
            n = cone.size(jj);
            rr = size(par.P1{sumlen+jj},2);
            if (rr > 0 && rr < n)
                if (rr <= n/2)
                    tmp0 = par.P1{sumlen+jj}'*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                    tmp1 = (tmp0*par.P1{sumlen+jj})*par.P1{sumlen+jj}';
                    tmp2 = par.Dsch12{sumlen+jj}.*(tmp0*par.P2{sumlen+jj});
                    tmp2 = tmp2*par.P2t{sumlen+jj};
                    tmp3 = par.P1{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
                    tmpcone = cone;
                    tmpcone.size = cone.size(jj);
                    %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},tmp3+tmp3');
                    By = By + par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),tmp3+tmp3');
                else
                    tmp0 = par.P2t{sumlen+jj}*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                    tmp1 = (tmp0*par.P2{sumlen+jj})*par.P2t{sumlen+jj};
                    tmp2 = (par.Dsch1{sumlen+jj}-par.Dsch12{sumlen+jj}').*(tmp0*par.P1{sumlen+jj});
                    tmp2 = tmp2*par.P1t{sumlen+jj};
                    tmp3 = par.P2{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
                    tmpcone = cone;
                    tmpcone.size = cone.size(jj);
                    %                By = By + par.sig*AXfun_sdpnal(pblk,At{p},par.Dsch1{p}*Aty-tmp3-tmp3');
                    By = By + par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),par.Dsch1{sumlen+jj}*Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))-tmp3-tmp3');
                end
            elseif (rr == n)
                tmpcone = cone;
                tmpcone.size = cone.size(jj);
                %             By = By + par.sig*AXfun_sdpnal(pblk,At{p},Aty);
                By = By + par.Dsch1{sumlen+jj}*par.sig*AXfun(tmpcone,At{p}(cumsumsquare(jj)+1:cumsumsquare(jj+1),:),Aty(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1)));
            end
        end
    elseif strcmp(cone.type,'q')
        if ~isfield(par, 'socp_formula') || par.socp_formula == 1 || par.socp_formula == 2
            Aty = At{p}*y;
            tmp1 = par.Q1{p}'*Aty.*par.Dsch1{sumlen+1};
            tmp1 = repelem(tmp1,cone.size,1).*par.P1{sumlen+1};
            tmp2 = par.Q2{p}'*Aty.*par.Dsch2{sumlen+1};
            tmp2 = repelem(tmp2,cone.size,1).*par.P2{sumlen+1};
            tmp3 = repelem(par.shift{p},cone.size,1).*Aty;
            tmp = tmp1+tmp2+tmp3;
            By = By + par.sig*(tmp'*At{p})';
        elseif par.socp_formula == 3
            Aty = At{p}*y;
            tmp1 = par.Q3{p}'*Aty.*par.Dsch3{p};
            tmp1 = repelem(tmp1,cone.size,1).*par.P3{p};
            tmp3 = par.shift3{p}.*Aty;
            tmp = tmp1+tmp3;
            By = By + par.sig*(tmp'*At{p})';
        end
    elseif strcmp(cone.type,'l') || strcmp(cone.type,'b') || strcmp(cone.type,'u') || strcmp(cone.type,'b2l')
        if (~isempty(par.Dsch12{sumlen+1}))
            tmp = par.Dsch12{sumlen+1}.*(At{p}*y);
            By = By + par.sig*(tmp'*At{p})';
        end
    elseif strcmp(cone.type,'u')
        tmp = At{p}*y;
        By = By + par.sig*(tmp'*At{p})';
    end
    sumlen = sumlen + strcmp(cone.type, 's') * length(cone.size) + (1 - strcmp(cone.type, 's')) * 1;
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
