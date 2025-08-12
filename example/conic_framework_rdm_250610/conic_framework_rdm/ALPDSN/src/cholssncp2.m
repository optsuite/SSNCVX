function [dy,relres,flag] = cholssncp2(K,At,rhs1,par,L,pcgTol,AL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m  =  length(rhs1);
Ntmp = [];
y = rhs1;
% AL.R = speye(m,m);
% AL.Rt = AL.R;
%%

%%
m  =  length(rhs1);
Ntmp = [];
y = rhs1;
% AL.R = speye(m,m);
% AL.Rt = AL.R;
schur = speye(m,m)*par.epsilon*AL.Rt*AL.R;
% schur = zeros(m,m);
if isempty(AL); existAL = 0; else; existAL = 1; end
if ~AL.isidentity && (existAL)
    if strcmp(AL.matfct_options,'chol')% R'*R = A*A'
        rhs1(AL.perm) = AL.Rt * rhs1;
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        rhs1(AL.perm) = mexbwsolve(inv(AL.R),rhs1);
    end
end
%%
for p = 1: length(K)
    cone = K{p};
    
    if strcmp(cone.type,'s')
      if ~issparse(par.Q1{p}) ; par.Q1{p} = sparse(par.Q1{p}); end
      if ~issparse(par.Q2{p}) ; par.Q1{p} = sparse(par.Q2{p}); end
      pblk = {K{p}.type,K{p}.size};
%       pblk2 = {par.cone.type,par.cone.size }
      tmp1 = mexskron(pblk,par.Q1{p}',par.Q1{p}');
      tmp2 = mexskron(pblk,par.Q1{p}',par.Q2{p}');
%       tmp3 = mexskron(pblk,par.Q2{p}',par.Q1{p}');

%       tmp = tmp1 + 2*tmp2 ;
      schurtmp1 = tmp1*At{p,1}.*sqrt(svecd(pblk, par.Sigma1{p} ));
      schurtmp2 = 2*tmp2*At{p,1}.*sqrt(svecd(pblk, par.Sigma2{p} + par.Sigma2{p}'));
      schurtmp = schurtmp1'*schurtmp1 + schurtmp2'*schurtmp2;     
%       tmp1 = mexskron(pblk,par.Q1{p}',par.Q1{p}');
%       tmp2 = mexskron(pblk,par.Q1{p}',par.Q2{p}');
% %       tmp3 = mexskron(pblk,par.Q2{p}',par.Q1{p}');
% 
%       tmp = tmp1 + 2*tmp2 ;
%       schurtmp = tmp*At{p,1}.*sqrt(svecd(pblk,par.Sigma1{p} + par.Sigma2{p} + par.Sigma2{p}'));
%       schurtmp = schurtmp'*schurtmp;     

%       tmp1 = mexskron(pblk,par.Q1{p}',par.Q1{p}');
%       tmp2 = mexskron(pblk,par.Q1{p}',par.Q2{p}');
%       tmp3 = mexskron(pblk,par.Q2{p}',par.Q1{p}');
% 
% %       tmp = tmp1 + tmp2 + tmp3 ;
%       schurtmp1 = tmp1*At{p,1}.*sqrt(svecd(pblk,par.Sigma1{p}));
%       schurtmp2 = tmp2*At{p,1}.*sqrt(svecd(pblk, par.Sigma2{p} ));
%       schurtmp3 = tmp3*At{p,1}.*sqrt(svecd(pblk, par.Sigma2{p}'));
%       schurtmp = schurtmp1 + schurtmp2 + schurtmp3;
%       schurtmp = schurtmp'*schurtmp; 



      %% schurtmp = 0.5*(schurtmp + schurtmp');
%       if (norm(par.permA(p,:)-[1:m]) > 0)
%           Perm = spconvert([(1:m)', par.permA(p,:)', ones(m,1)]);
%           schur = schur + Perm'*schurtmp*Perm;
%       else
        schur = schur + schurtmp;
    elseif strcmp(cone.type,'l')
        schurtmp = par.Dsch12{p}.*(At{p});
        schur = schur + par.sig*(schurtmp'*At{p})';
    end

end
          if ~issparse(schur); schur = sparse(schur); end;
          L.matfct_options = 'spchol'; 
          [L.R,indef,L.perm] = chol(schur,'vector'); 
          L.Rt  = L.R';
          L.matdim = length(schur); 
          diagR = full(diag(L.R)).^2; 
          coeff.mat11 = schur;
          coeff.mat22 = [];
          coeff.mat12 = [];
          [dy,resnrm,solve_ok] = symqmr(coeff,rhs1,L,[],[],10);
          1;
% schur2 = zeros(m,m);
% for ii = 1:m
%     aa = sparse(m,1);
%     aa(ii) = 1;
%     schur2(:,ii) = matvec_y2mit(K, At, par, aa, AL );
% end
% 1;


%%
% symm = 0;
% for p = 1: length(K)
%     cone = K{p};
%     if strcmp(cone.type,'s')
%     n  = par.K{p}.size;
%     for jj = 1:length(n)
%     Sigma{p} = zeros(n(jj),n(jj));
%     r = size(par.Dsch12{p},1);
%     if r > 0
%     Sigmatmp = [par.Dsch1{p}*ones(r,r) par.Dsch12{p}];
%     Sigma{p}(1:r,:) = Sigmatmp;
%     Sigma{p}(:,1:r) = Sigmatmp';
%     else
%        Sigma{p} = 0;
%     end
%     Q{p} = [par.P1{p} par.P2{p}];
%     Ntmp = [];
%     blk = {K{p}.type,K{p}.size};
%     pblk = blk;
%    
%     for k =1:m
%         isspAk = 0;
%         Ak = mexsmat(blk,At{p},isspAk,1,k);
%         Ak = Ak;
%         idx = 1:k;
% %         tmp1 = Prod3(pblk,par.P1{p}',Ak,par.P1{p},symm);
% %         tmp2 = Prod3(pblk,par.P1{p}',Ak,par.P2{p},symm);
% %         tmp1 = tmp1.*sqrt(par.Dsch1{p});
% %         tmp2 = tmp2.*sqrt(par.Dsch2{p});
% %         tmp = [tmp1,tmp2;tmp2',sparse(size(tmp2,2),size(tmp2,2)) ];
%         tmp = Prod3(pblk,Q{p}',Ak,Q{p},symm);
%         tmp = tmp.*sqrt(Sigma{p});
%         tmp = svec(pblk,tmp);
%         Ntmp = [Ntmp, tmp];
% %         tmp2 = schur(idx,k) + Ntmp'*tmp;
%         tmp2 = schur(idx,k) + mexinprod(blk,Ntmp,tmp,k,p);
%         schur(idx,k) = tmp2;
%         schur(k,idx) = tmp2';
%     end
%     1;
%     end
%     elseif strcmp(cone.type,'l')
%         if (~isempty(par.Dsch12{p}))
%             schurtmp = At{p}'.*(par.Dsch12{p}.*At{p});
%             schur = schur + schurtmp;
%         end
%     end
% 
% end

% R = chol(schur);

% dy = schur \ rhs1;
% relres2 = norm(schur*dy - rhs1,'fro');
flag = 1;

if ~AL.isidentity && (existAL)
    if strcmp(AL.matfct_options,'chol')
        dy = AL.R * dy(AL.perm);
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        dy = mexfwsolve(AL.R,dy(AL.perm));
    end
end

% relres = norm(y - matvec_y2mit(K, At, par, dy, AL ) )/ (1 + norm(rhs1));
relres =  resnrm;

1;
end

%%*******************************************************************
%% schurmat_sblk: compute Schur complement matrix corresponding to 
%%                SDP blocks. 
%%
%% symm = 0, HKM
%%      = 1, NT
%%*****************************************************************
%% SDPT3: version 4.0
%% Copyright (c) 1997 by
%% Kim-Chuan Toh, Michael J. Todd, Reha H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************

%    function schur = schurmat_sblk(blk,At,par,schur,p,X,Y) 
% 
%    global  nnzschur nzlistschur
% 
%    iter = par.iter; 
%    smallblkdim = par.smallblkdim; 
% 
%    if isempty(smallblkdim); smallblkdim = 50; end
%    if (nargin == 7); symm = 0; else; symm = 1; Y = X; end; 
%    m = length(schur);    
%    pblk = blk(p,:); 
%    if (iter == 1)
%       nnzschur(size(blk,1),1) = m*m; 
%       nzlistschur = cell(size(blk,1),1); 
%    end
% %%
%    if (max(pblk{2}) > smallblkdim) | (length(pblk{2}) <= 10)
%       %%
%       %% compute schur for matrices that are very sparse. 
%       %%
%       m1 = size(At{p,1},2); 
%       if issparse(schur); schur = full(schur); end;  
%       J = 0;     
%       if (J > 0)
%          if issparse(X{p}) & ~issparse(Y{p}); X{p} = full(X{p}); end
%          if ~issparse(X{p}) & issparse(Y{p}); Y{p} = full(Y{p}); end
%          if (iter <= 3)     
%             [schur,nnzschur(p),nzlisttmp] = mexschur(pblk,At{p,1},par.nzlistA{p,1},...
%             par.nzlistA{p,2},par.permA(p,:),Y{p},X{p},J,symm,schur); 
%             if (nnzschur(p) == mexnnz(nzlisttmp)) 
%                nzlistschur{p} = nzlisttmp;
%             else
%                nzlistschur{p} = []; 
%             end
%          else
%             if isempty(nzlistschur{p})
%                schur = mexschur(pblk,At{p,1},par.nzlistA{p,1},...
%                par.nzlistA{p,2},par.permA(p,:),Y{p},X{p},J,symm,schur);
%             else
%                schur = mexschur(pblk,At{p,1},par.nzlistA{p,1},...
%                par.nzlistA{p,2},par.permA(p,:),Y{p},X{p},J,symm,schur,nzlistschur{p});
%             end 
%          end
%       end  
%       %%
%       %% compute schur for matrices that are not so sparse or dense.
%       %% 
%       if (m1 < m) %% for low rank constraints
%          ss = [0, cumsum(pblk{3})]; 
%          len = sum(pblk{3});
%          dd = At{p,3};
%          DD = spconvert([dd(:,2:4); len,len,0]);
%          XVD = X{p}*At{p,2}*DD; 
%          YVD = Y{p}*At{p,2}*DD;
%       end
%       L = max(find(par.nzlistAsum{p,1} < inf)) -1;  
%       if (J < L)
%          len = par.nzlistAsum{p,1}(J+1); list = par.nzlistAsum{p,2}(1:len,:); 
%       end      
%       if (m1 > 0)
%          for k = J+1:m 
%             if (k<=m1) 
%                isspAk = par.isspA(p,k);
%                Ak = mexsmat(blk,At,isspAk,p,k);
%                if (k <= L) 
%                   idx1 = par.nzlistAsum{p,1}(k)+1; idx2 = par.nzlistAsum{p,1}(k+1);
%                   list = [list; par.nzlistAsum{p,2}(idx1:idx2,:)]; 
%                   list = sortrows(list,[2 1]); 
%                   tmp = Prod3(pblk,X{p},Ak,Y{p},symm,list); 
%                else
%                   tmp = Prod3(pblk,X{p},Ak,Y{p},symm);
%                end
%             else %%--- for low rank constraints
%                idx = [ss(k-m1)+1 :ss(k-m1+1)]; 
%                tmp = XVD(:,idx)* (Y{p}*At{p,2}(:,idx))';
%             end
%             if (~symm)
%                tmp = 0.5*(mexsvec(pblk,tmp) + mexsvec(pblk,tmp'));
%             else
%                tmp = mexsvec(pblk,tmp);                
%             end 
%             permk = par.permA(p,k);   
%             idx  = par.permA(p,1:min(k,m1));             
%             tmp2 = schur(idx,permk) + mexinprod(blk,At,tmp,min(k,m1),p);
%             schur(idx,permk) = tmp2; 
%             schur(permk,idx) = tmp2';
%          end         
%       end
%       if (m1 < m) %% for low rank constraints
%          m2 = m - m1;
%          XVtmp = XVD'*At{p,2};
%          YVtmp = At{p,2}'*YVD;
%          for k = 1:m2
%             idx0 = [ss(k)+1 : ss(k+1)]; 
%             tmp = XVtmp(:,idx0) .* YVtmp(:,idx0);  
%             tmp = tmp*ones(length(idx0),1);             
%             tmp3 = schur(m1+[1:m2],m1+k) + mexqops(pblk{3},tmp,ones(length(tmp),1),1); 
%             schur(m1+[1:m2],m1+k) = tmp3; 
%          end
%       end
%    else  
%       %%--- for SDP block where each sub-block is small dimensional
%       if issparse(X{p}) & ~issparse(Y{p}); Y{p} = sparse(Y{p}); end
%       if ~issparse(X{p}) & issparse(Y{p}); X{p} = sparse(X{p}); end
%       tmp = mexskron(pblk,X{p},Y{p});
%       schurtmp = At{p,1}'*tmp*At{p,1};      
%       %% schurtmp = 0.5*(schurtmp + schurtmp');
%       if (norm(par.permA(p,:)-[1:m]) > 0)
%          Perm = spconvert([(1:m)', par.permA(p,:)', ones(m,1)]); 
%          schur = schur + Perm'*schurtmp*Perm;
%       else
%          schur = schur + schurtmp;
%       end    
%    end
%%*******************************************************************


