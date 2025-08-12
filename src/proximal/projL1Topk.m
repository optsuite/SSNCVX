function [y,rr] = projL1Topk(x,lambda,k)
% [tmp, rr] = proj_inf(x,lambda);
nrmx = norm(x,1);
if nrmx > k
    f = @(mu) sum(abs(sign(x).*max( min(abs(x/lambda) - mu,1),0 ))) - k;
    options = optimset('Display', 'off');
    mu = fsolve(f,0,options);
    tmp = sign(x).*min(max(abs(x/lambda) - mu,0) ,1);
    y = x - lambda*tmp;
    rr = abs(x/lambda) > mu & tmp > -1 & tmp < 1;
    rr = 1 - rr;
else
    y = min(max(x,-1),1);
    rr = y == x;
end 