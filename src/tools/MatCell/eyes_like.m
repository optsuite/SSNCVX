function Y = eyes_like(X)
    Y = cell(size(X));
    for i = 1:length(X)
        if size(X{i},2) ==1
        Y{i} = ones((size(X{i})));
        else
        Y{i} = eye((size(X{i})));
        end
    end
end