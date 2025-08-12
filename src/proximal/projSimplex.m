function [y,rr] = projSimplex(x,lambda)
% [tmp, rr] = proj_inf(x,lambda);
% nrmx = norm(x,1);
% if nrmx > lambda
    f = @(mu) sum(max(x/lambda - mu*ones(size(x)),0)) - 1;
    options = optimset('Display', 'off');
    mu = fsolve(f,0,options);
    tmp = max(x/lambda - mu ,0);
    y = x - lambda*tmp;
    rr = x/lambda > mu;
    rr = 1 - rr;
% else
%     y = x;
%     rr = ones(length(x),1);
% end
% rr = ~rr;