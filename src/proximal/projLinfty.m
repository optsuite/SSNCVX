function [y,rr] = projLinfty(x,lambda)
% [tmp, rr] = proj_inf(x,lambda);
nrmx = norm(x,1);
if nrmx > lambda
    f = @(mu) sum(abs(sign(x).*max(abs(x/lambda) - mu*ones(size(x)),0))) - 1;
    options = optimset('Display', 'off');
    mu = fsolve(f,0,options);
    tmp = sign(x).*max(abs(x/lambda) - mu ,0);
    y = x - lambda*tmp;
    rr = abs(x/lambda) > mu;
    rr = 1 - rr;
else
    y = zeros(size(x));
    rr = zeros(size(x));
end
% rr = ~rr;