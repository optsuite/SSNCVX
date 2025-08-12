function q = fwsolve(L,r)
    if (L.isidentity)
        q = r;
    else
        if strcmp(L.matfct_options,'chol')
            q = mextriang(L.R,r(L.perm),2);
        elseif strcmp(L.matfct_options,'spcholmatlab')
            q = mexfwsolve(L.R,r(L.perm,1));
        end
    end
end
    