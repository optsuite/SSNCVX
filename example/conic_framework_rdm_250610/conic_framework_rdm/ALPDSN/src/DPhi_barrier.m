%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-09 16:59:05
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%%*************************************************************************
%% compute P*(Dsch.*(Pt*H*P))*Pt + shift * H

%% SDPNAL: 
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh 
%%************************************************
%%*************************************************************************  

function Y = DPhi_barrier(K, par, H)
 
   Y = MatCell(length(K));
   for p = 1: length(K)
      cone = K{p};
      n = sum(cone.size); 
      if strcmp(cone.type,'s')
         % 
         % Y = P * (Sigma .* (P'* H * P)) * P'
         % where Sigma = [Dsch1{p}*I, Dsch2{p}; 
         %                Dsch2{p}', 0];
         %  note that Dsch1{p} is a scalar here
         Y{p} = sparse(n,n); 
        tmp0 = par.P1{p}'*H{p};
        %                 tmp1 = (tmp0*par.P1{p})*par.P1{p}';
        tmp2 = par.Dsch2{p}.*(tmp0*par.P1{p});
        tmp2 = tmp2*par.P1t{p};
        Y{p} = par.P1{p}* tmp2;


%          Y{p} = Y{p} + par.shift{p} * H{p};
      elseif strcmp(cone.type,'q')
         % For one subblock, Y = Dsch11 * I + Dsch2(1) * P1 * P1' + Dsch2(2) * P2 * P2' 
         % For multiblock 
         % Y = diag(repelem(Dsch11, cone.size, 1) ) 
         %   + blk_spdiag(P1, cone.size) * diag(Dsch2(:, 1)) * blk_spdiag(P1, cone.size)' 
         %   + blk_spdiag(P2, cone.size) * diag(Dsch2(:, 2)) * blk_spdiag(P2, cone.size)'
         Y{p} = sparse(n,1); 
         if (~isempty(par.Dsch2{p})) 
            P1blkdiag = blk_spdiag(par.P1{p}, cone.size);
            P2blkdiag = blk_spdiag(par.P2{p}, cone.size);
            Y{p} = repelem(par.shift{p}, cone.size, 1) .* H{p} + ...
                  P1blkdiag * spdiag(par.Dsch1{p}) * (P1blkdiag' * H{p})+ ...
                  P2blkdiag * spdiag(par.Dsch2{p}) * (P2blkdiag' * H{p});
         end
      elseif strcmp(cone.type,'l')
         Y{p} = sparse(n,1); 
         if (~isempty(par.Dsch2{p}))
            Y{p} = par.Dsch2{p}.*H{p}; 
         end
      elseif strcmp(cone.type,'u')
         Y{p} = H{p}; 
      end
   end   
