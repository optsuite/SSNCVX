

function  [z1,z2,resnrm,solve_ok] = psqmrplus2(matvecfname,K,At,b1,b2,par,L,tol,maxit,AL)

N1 = length(b1);
y0 = zeros(N1,1);
zb0 = zeros(size(b2{1}));

solve_ok = 1;
miniter = 1;
printlevel = 0;
if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end
if isfield(par,'printlevel'); printlevel = par.printlevel; end

%%
z1 = y0;
z2 = zb0;
r1 = b1;
r2 = b2{1} ;
err = norm(r1) ; resnrm(1) = err; minres = err;
%%
q1 = r1;
q2 = r2 ;
% tau_old = sqrt(q1(:)'*q1(:)+q2(:)'*q2(:));
tau_old = norm(q1)+sqrt(sum(q2.*q2,"all"));
rho_old  = r1'*q1+r2(:)'*q2(:);
theta_old = 0;
d1 = zeros(N1,1);
d2 = zeros(size(b2{1}));
res1 = r1; res2 = r2;
Ad1 = zeros(N1,1);
Ad2 = zeros(size(b2{1}));
%%
%% main loop
%%
tiny = 1e-30;
for iter = 1:maxit
    [Aq1,Aq2] = feval(matvecfname,K,At,par,q1,q2,AL);
    sigma = q1'*Aq1 + sum(q2(:).*Aq2(:),"all");
    if (abs(sigma) < tiny)
        solve_ok = 2;
        if (printlevel); fprintf('s1'); end
        break;
    else
        alpha = rho_old/sigma;
        r1 = r1 - alpha*Aq1;
        r2 = -alpha*Aq2 + r2;
    end
    u1 = r1;
    u2 = r2;
    %%
           theta = (norm(u1)+sqrt(sum(u2.*u2,"all")))/tau_old;
%     theta = sqrt((sum(u1.*u1,"all")+sum(u2.*u2,"all")))/tau_old;
    c = 1/sqrt(1+theta^2);
    tau = tau_old*theta*c;
    gam = (c^2*theta_old^2);
    eta = (c^2*alpha);
    d1 = gam*d1 + eta*q1;

    d2 = gam*d2+eta*q2;
    z1 = z1 + d1;
    z2 = z2 + d2;

    %%----- stopping conditions ----
    Ad1 = gam*Ad1 + eta*Aq1;
    res1 = res1 - Ad1;


    Ad2 = gam*Ad2+eta*Aq2;
    res2 = res2-Ad2;
%     err = sqrt(res1(:)'*res1(:)+sum(res2.*res2,"all"));
           err = norm(res1) + sqrt(sum(res2.*res2,"all"));
    resnrm(iter+1) = err;
    if (err < minres); minres = err; end
    if (err < tol) && (iter > miniter) && (b1'*z1 + sum(b2{1}.*z2,'all') > 0)
        break;
    end

    if (abs(rho_old) < tiny)
        solve_ok = 2;
        fprintf('s2');
        break;
    else
        %           rho  = r1'*u1 + r2(:)'*u2(:);
        rho  = r1'*r1 + sum(r2.*r2,"all");
        beta = rho/rho_old;
        q1 = r1 + beta*q1;
        q2 = r2+beta*q2;
    end
    rho_old = rho;
    tau_old = tau;
    theta_old = theta;
end
if (iter == maxit); solve_ok = -2; end
if (solve_ok ~= -1)
    if (printlevel); fprintf(' '); end
end

