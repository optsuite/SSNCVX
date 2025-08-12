function [value] = pvalue(pblk, X, trans)
plength = length(pblk);
value = 0;
if ~iscell(X)
    tmpX = X;
    clear X
    X = {tmpX};
end
for i =1:plength
    if  strcmp(pblk{i}.type,'l1')

        if isfield(pblk{i},'shift' )
        value = value + pblk{i}.coefficient*norm(reshape(X{i} - pblk{i}.shift,[],1),1);
        else
        value = value + pblk{i}.coefficient*norm(reshape(X{i},[],1),1); 
        end
    elseif strcmp(pblk{i}.type,'square')
        if isfield(pblk{i},'shift' )
        value = value + pblk{i}.coefficient*norm( X{i} - pblk{i}.shift,'fro')^2;
        else
        value = value + pblk{i}.coefficient*norm( X{i},'fro')^2;
        end
    elseif strcmp(pblk{i}.type,'linfty')
        value = value + pblk{i}.coefficient*norm(X{i},'inf');
    elseif strcmp(pblk{i}.type,'fused')
        tmpX = reshape(X{i},[],1);
        value = value + pblk{i}.coefficient*norm(tmpX,1) + pblk{i}.coefficient2*(norm(tmpX(1:end-1) - tmpX(2:end),1));
    elseif strcmp(pblk{i}.type,'box') % We leave the constraint function as empty since the iteration point may not be feasilble
    elseif strcmp(pblk{i}.type,'s')
    elseif strcmp(pblk{i}.type,'l')
    elseif strcmp(pblk{i}.type,'q')
    elseif strcmp(pblk{i}.type,'l1con')
    elseif strcmp(pblk{i}.type,'l2con')
    elseif strcmp(pblk{i}.type,'linftycon')
    elseif strcmp(pblk{i}.type,'l1l2')
        for j = 1:size(X{i},2)
        value = value + pblk{i}.coefficient*norm(X{i}(:,j),2);
        end
    elseif strcmp(pblk{i}.type,'topk')
        descent = sort(X{i},'descend');
        value = value + pblk{i}.coefficient*sum(descent(1:pblk{i}.topk));
    elseif strcmp(pblk{i}.type,'l1topk')
        descent = sort(abs(X{i}),'descend');
        value = value + pblk{i}.coefficient*sum(descent(1:pblk{i}.topk));
    elseif strcmp(pblk{i}.type,'max')
        value = value + pblk{i}.coefficient*max(X{i});
    elseif strcmp(pblk{i}.type,'l2')
        if isfield(pblk{i},'shift' )
        value = value + pblk{i}.coefficient*norm(X{i} - pblk{i}.shift,2);
        else
        value = value + pblk{i}.coefficient*norm(X{i},2);
        end
    elseif strcmp(pblk{i}.type,'exp')
        value = value + pblk{i}.coefficient*sum(exp(X{i}),'all');
    elseif strcmp(pblk{i}.type,'log')
        value = value - pblk{i}.coefficient*sum(log(X{i}),'all');
    % elseif strcmp(pblk{i}.type,'logsumexp')
    %     value = value + pblk{i}.coefficient*log(sum(exp(X{i}),'all'));
    elseif strcmp(pblk{i}.type,'nuclear')
        S = svd(X{i});
        value = value + pblk{i}.coefficient*sum(S,'all');
    elseif strcmp(pblk{i}.type,'l2l2')
        value = value + pblk{i}.coefficient*norm(X{i},2);
    elseif strcmp(pblk{i}.type,'logdet')
        value = value - pblk{i}.coefficient*log(det(X{i}));
    elseif strcmp(pblk{i}.type,'huber')
        value = value + pblk{i}.coefficient*( 0.5*sum((X{i}.*(abs(X{i})<pblk{i}.coefficient2)).^2,'all') + sum(pblk{i}.coefficient2*(abs(X{i}) - pblk{i}.coefficient2/2).*(abs(X{i})>pblk{i}.coefficient2),'all'));
    else %user-defined
        value = value + pblk{i}.coefficient*pblk{i}.pobj(X{i});
    end
end
end