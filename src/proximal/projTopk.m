function [y,rr] = projTopk(x,lambda,k)
% [tmp, rr] = proj_inf(x,lambda);
nrmx = norm(x,1);
if nrmx > k
    f = @(mu) sum(max( min(x/lambda - mu,1),0 )) - k;
    options = optimset('Display', 'off');
    mu = fsolve(f,0,options);
    tmp = max(min(x/lambda - mu,1) ,0);
    y = x - lambda*tmp;
    rr = x/lambda > mu;
    rr = 1 - rr;
else
    y = x;
    rr = ones(length(x),1);
end