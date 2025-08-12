function [value] = pdual_value(pblk, X, params)
plength = length(pblk);
value = 0;
eps = 1e-10;
for i =1:plength
    if  strcmp(pblk{i}.type,'l1')
        % if max(X{i}) <= pblk{i}.coefficient + eps
        %     tmpvalue = 0;
        % else
        %     tmpvalue = inf;
        % end
        % value = value + tmpvalue;
        if isfield(pblk{i},'shift')
            value = value + dot_ssn(pblk{i}.shift,X{i});
        end
    elseif strcmp(pblk{i}.type,'linfty')
        % if norm(X{i},1) <= pblk{i}.coefficient + eps
        %     tmpvalue = 0;
        % else
        %     tmpvalue = inf;
        % end
        % value = value + tmpvalue;
    elseif strcmp(pblk{i}.type,'linfcon')
        value = value + pblk{i}.coefficient*abs(X{i});
    elseif strcmp(pblk{i}.type,'box')
        xP = max(X{i},0);
        xN = min(X{i},0);
        value = value + dot_ssn_inf(xN,pblk{i}.l) + dot_ssn_inf(xP,pblk{i}.u) ;
    elseif strcmp(pblk{i}.type,'square')
        if isfield(pblk{i},'shift')
            value = value + dot_ssn(pblk{i}.shift,X{i}) + 0.5^2/pblk{i}.coefficient*norm(X{i},'fro')^2;
        else
            value = value + 0.5^2/pblk{i}.coefficient*norm(full(X{i}),'fro')^2;
        end
    elseif strcmp(pblk{i}.type,'l2')
        % if norm(X{i}) <= pblk{i}.coefficient + eps
        tmpvalue = 0;
        % else
        %     tmpvalue = inf;
        % end
        value = value + tmpvalue;
        if isfield(pblk{i},'shift')
            value = value + dot_ssn(pblk{i}.shift,X{i});
        end
    elseif strcmp(pblk{i}.type,'exp')
        Xitmp = X{i}/pblk{i}.coefficient;
        value = value + pblk{i}.coefficient*sum(Xitmp.*log(Xitmp) - Xitmp,'all' );
    elseif strcmp(pblk{i}.type,'log')
        Xitmp = -X{i}/pblk{i}.coefficient;
        value = value - pblk{i}.coefficient*( size(Xitmp,1)+ sum(log(Xitmp),'all' ) - size(Xitmp,1)*log(pblk{i}.coefficient) );
        % elseif strcmp(pblk{i}.type,'logsumexp')
        %     Xitmp = X{i}/pblk{i}.coefficient;
        %     value = value + sum(X{i}.*log(Xitmp) - pblk{i}.coefficient*log(pblk{i}.coefficient),'all' );
    elseif strcmp(pblk{i}.type,'linftycon')
        Xitmp = X{i};
        if isfield(pblk{i},'shift')
            value = value  + dot_ssn(pblk{i}.shift,X{i}) + pblk{i}.coefficient*norm(Xitmp,1);
        else
            value = value + pblk{i}.coefficient*norm(Xitmp,1);
        end
    elseif strcmp(pblk{i}.type,'l1topk')
    elseif strcmp(pblk{i}.type,'max')
    elseif strcmp(pblk{i}.type,'nuclear')
    elseif strcmp(pblk{i}.type,'l2con')
        Xitmp = X{i};
        if isfield(pblk{i},'shift')
            value = value  + dot_ssn(pblk{i}.shift,X{i}) + pblk{i}.coefficient*norm(Xitmp,2);
        else
            value = value + pblk{i}.coefficient*norm(Xitmp,2);
        end
    elseif strcmp(pblk{i}.type,'l1con')
        Xitmp = X{i};
        if isfield(pblk{i},'shift')
            value = value  + dot_ssn(pblk{i}.shift,X{i}) + pblk{i}.coefficient*norm(Xitmp,'inf');
        else
            value = value + pblk{i}.coefficient*norm(Xitmp,'inf');
        end
    elseif strcmp(pblk{i}.type,'logdet')
        Xitmp = X{i};
        value = value + pblk{i}.coefficient*(-size(Xitmp,1) - log(det(-Xitmp/pblk{i}.coefficient)));
    elseif strcmp(pblk{i}.type,'huber')
        if isfield(pblk{i},'shift')
            value = value + dot_ssn(pblk{i}.shift,X{i}) + 0.5^2/pblk{i}.coefficient*norm(X{i},'fro')^2;
        else
            % if norm(X{i},'inf') <= pblk{i}.coefficient2 + eps
            value = value + 0.5/pblk{i}.coefficient*norm(full(X{i}),'fro')^2;
            % else
            %     value = inf;
            % end
        end

    end
end
end