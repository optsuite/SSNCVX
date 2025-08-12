%  Tao Wei: Edit 'q' cone on 2024-09-19
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

function Y = DPhimit(K, par, H)

Y = MatCell(length(K));
sumlen = 0;
for p = 1: length(K)
    cone = K{p};
    
    if strcmp(cone.type,'s')
        %       sumlen = sumlen + length(cone.size);
        cumsumlen = [0 cumsum(cone.size)];
        cumsumlens = [0 cumsum(cone.size.^2)];
        n = sum(cone.size);
        n2 = sum(cone.size.^2);
        %
        % Y = P * (Sigma .* (P'* H * P)) * P'
        % where Sigma = [Dsch1{p}*I, Dsch2{p};
        %                Dsch2{p}', 0];
        %  note that Dsch1{p} is a scalar here
        Y{p} = sparse(n,n);
        row_index = zeros(n2,1);
        col_index = zeros(n2,1);
        value_index = zeros(n2,1);

        for  jj = 1:length(cone.size)
            nn = cone.size(jj);
            rr = size(par.P1{sumlen+jj},2);
            cj1 = cumsumlen(jj)+1;
            cj2 = cumsumlen(jj+1);
            cjs1 = cumsumlens(jj) + 1;
            cjs2 = cumsumlens(jj+1);
            [row, col] = meshgrid(cj1:cj2, cj1:cj2) ;
            row_index(cjs1:cjs2 ) = row(:);
            col_index(cjs1:cjs2) = col(:);  
            if (rr > 0 && rr < nn)
                if (rr <= nn/2)
                    tmp0 = par.P1{sumlen+jj}'*H{p}(cj1:cj2,cj1:cj2) ;                  % tmp0: U = Q_{\alpha}^\top H
                    tmp1 = (tmp0*par.P1{sumlen+jj})*par.P1{sumlen+jj}';      % par.P1{p}: Q_{\alpha}
                    tmp2 = par.Dsch2{sumlen+jj}.*(tmp0*par.P2{sumlen+jj});  % par.Dsch2{p}: \nu_{\alpha \bar{\alpha}}; par.P2{p}: Q_{\bar{\alpha}}
                    tmp2 = tmp2*par.P2{sumlen+jj}';                  % temp2: \nu_{\alpha \bar{\alpha}} \circ (U Q_{\bar{\alpha}}) Q_{\bar{\alpha}}^\top
                    tmp3 = par.P1{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
%                     Y{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))  = tmp3+tmp3';                 
                    tmp3 = tmp3+tmp3';     
                    value_index(cjs1:cjs2) =  tmp3(:);
                else
                    tmp0 = par.P2{sumlen+jj}'*H{p}(cj1:cj2,cj1:cj2) ;
                    tmp1 = (tmp0*par.P2{sumlen+jj})*par.P2{sumlen+jj}';
                    tmp2 = (par.Dsch1{sumlen+jj}-par.Dsch2{sumlen+jj}').*(tmp0*par.P1{sumlen+jj});
                    tmp2 = tmp2*par.P1{sumlen+jj}';
                    tmp3 = par.P2{sumlen+jj}*(0.5*par.Dsch1{sumlen+jj}*tmp1 + tmp2);
%                     Y{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1)) = par.Dsch1{sumlen+jj}*H{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))-tmp3-tmp3';
                    tmp3 = par.Dsch1{sumlen+jj}*H{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))-tmp3-tmp3';
%                     tmp3 = tmp3+tmp3';     
                    value_index(cjs1:cjs2 ) = tmp3(:);
                end
            elseif (rr == nn)
%                 Y{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1))  = par.Dsch1{sumlen+jj}*H{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                tmp3 = par.Dsch1{sumlen+jj}*H{p}(cumsumlen(jj)+1:cumsumlen(jj+1),cumsumlen(jj)+1:cumsumlen(jj+1));
                value_index(cumsumlens(jj)+1:cumsumlens(jj+1) ) = tmp3(:);
            end
        end
        Y{p} = sparse(row_index,col_index,value_index,n,n);
        1;
    elseif strcmp(cone.type,'q')
        if ~isfield(par, 'socp_formula') || par.socp_formula == 1 || par.socp_formula == 2
            n = sum(cone.size);
            Y{p} = sparse(n,1);
            tmp1 = par.Q1{p}'*H{p}.*par.Dsch1{sumlen+1};
            tmp1 = repelem(tmp1,cone.size,1).*par.P1{sumlen+1};
            tmp2 = par.Q2{p}'*H{p}.*par.Dsch2{sumlen+1};
            tmp2 = repelem(tmp2,cone.size,1).*par.P2{sumlen+1};
            tmp3 = repelem(par.shift{p},cone.size,1).*H{p};
            Y{p} = tmp1+tmp2+tmp3;
        elseif par.socp_formula == 3
            n = sum(cone.size);
            Y{p} = sparse(n,1);
            if ~isfield(par, 'Q3')
                ncols = length(cone.size);
                nrows = sum(cone.size);
                Prow_index = 1:nrows;
                Pcol_index = repelem(1:ncols, cone.size);
                Pvalue1 = par.P3{p}(:);

                Q3 = sparse(Prow_index, Pcol_index, Pvalue1, nrows, ncols);
            else
                Q3 = par.Q3{p};
            end
            tmp1 = Q3'*H{p}.*par.Dsch3{sumlen+1};
            tmp1 = repelem(tmp1,cone.size,1).*par.P3{sumlen+1};
            tmp3 = par.shift3{p}.*H{p};
            Y{p} = tmp1+tmp3;
        end
    elseif strcmp(cone.type,'l') || strcmp(cone.type,'b') || strcmp(cone.type,'u') || strcmp(cone.type,'b2l')
        n = sum(cone.size);
        Y{p} = sparse(n,1);
        if (~isempty(par.Dsch2{sumlen+1}))
            Y{p} = par.Dsch2{sumlen+1}.*H{p};
        end
    elseif strcmp(cone.type,'u')
        Y{p} = H{p};
    end
    sumlen = sumlen + strcmp(cone.type, 's') * length(cone.size) + (1 - strcmp(cone.type, 's')) * 1;
end

