%% gradientf: compuate the gradient of the dual function of f
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
function  value = dual_gradientf(f, X, trans)
plength = length(f);
for i =1:plength
    if  strcmp(f{i}.type,'square')
        if isfield(f{i},'shift')
        value{i,1} = -0.5/f{i}.coefficient*(X{i} ) + f{i}.shift;
        else
        value{i,1} = -0.5/f{i}.coefficient*(X{i});
        end
    elseif strcmp(f{i}.type,'exp')
        Xtmp = X{i}/f{i}.coefficient;
        if isfield(f{i},'shift')
            value{i,1} = log(max(-Xtmp,1e-7)) + f{i}.shift; 
        else
        value{i,1} = log(max(-Xtmp,1e-7));
        end
    elseif strcmp(f{i}.type,'log')
        Xtmp = X{i}/f{i}.coefficient;
        if isfield(f{i},'shift')
            value{i,1} = f{i}.coefficient./Xtmp + f{i}.shift; 
        else
        value{i,1} = f{i}.coefficient./Xtmp;
        end
    elseif strcmp(f{i}.type,'logdet')
        if iscell(X)
            Xtmp = X{i};
        else
            Xtmp = X;
        end
         if isfield(f{i},'shift')
             value{i,1} = inv(Xtmp/f{i}.coefficient) + f{i}.shift; 
         else
        value{i,1} = inv(Xtmp/f{i}.coefficient);
         end
    end
end
end
