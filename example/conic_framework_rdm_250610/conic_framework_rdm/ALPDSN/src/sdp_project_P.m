

function Xne = sdp_project_P(K,X,trans)

Xne = MatCell(length(X));

%%
for p = 1: length(K)
    cone = K{p} ;
    if strcmp(cone.type,'s')
        Xptmp = X{p};
%         Xptmp(Xptmp<trans.xL{p}) = 0;
%         Xptmp(Xptmp>trans.xU{p}) = 0;
        Xptmp = max(Xptmp,trans.xL{p});
%                end
%                 if trans.isU
        Xptmp = min(Xptmp,trans.xU{p});
    elseif strcmp(cone.type,'l')
        n  = sum(cone.size);
        Xptmp = zeros(n,1);
        %%***** perturbation *****
        %Dsch12 = addtol*ones(n,1);
        %           posidx1 = find(X{p} > tol);
        %           posidx2 = find(X{p} < tol);
        %           if ~isempty(posidx)
        %              Xptmp(posidx)  = abs(X{p}(posidx));
        %           end
        Xptmp = X{p};
        Xptmp(Xptmp<trans.xL{p}) = 0;
        Xptmp(Xptmp>trans.xU{p}) = 0;
    elseif strcmp(cone.type,'u')
        Xptmp = X{p};
        %%***** perturbation *****
    end
    Xne{p}  = Xptmp;
end
end
%%***************************************************************************

