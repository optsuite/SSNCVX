function out = spdiag(input)
    if length(input) == 0
        out = sparse(0,0);
    else
        out = spdiags(input(:),0,length(input),length(input));
    end
    
end