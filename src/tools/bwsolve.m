function q = bwsolve(L,r)
    if (L.isidentity)
        q = r;
    else
        if strcmp(L.matfct_options,'chol')
            q(L.perm,1) = mextriang(L.R, r, 1);
        elseif strcmp(L.matfct_options,'spcholmatlab')
            q(L.perm,1) = mexbwsolve(L.Rt, r);
        end
    end
end