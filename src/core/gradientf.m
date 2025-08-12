function  value = gradientf(f, X, trans)
plength = length(f);
for i =1:plength
    if  strcmp(f{i}.type,'square')
        if isfield(f{i},'shift')
        value{i} = 2*f{i}.coefficient*(X{i} ) - f{i}.shift;
        else
        value{i} = 2*f{i}.coefficient*(X{i});
        end
    elseif strcmp(f{i}.type,'exp')
        value{i} = -f{i}.coefficient*log(max(-X{i},1e-7));
    elseif strcmp(f{i}.type,'logdet')
        value{i} = f{i}.coefficient*inv(X{i});
    elseif strcmp(f{i}.type,'log')
        value{i} = f{i}.coefficient*1./X{i};
    end
end
end