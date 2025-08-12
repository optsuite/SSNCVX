%%**********************************************************************
%%**********************************************************************

%   (AP'*AP + mu I) d = rhs
    function [xi,resnrm,solve_ok,solver] = Classic_Lasso_linsys_solver1(opts,rhs,par)
    
    m = length(rhs);
    pp = ~par.rr;
    
    mu = par.mu;
    Ayes = isfield(opts,'A');
    if Ayes
        AP = - par.DD * par.AP;
    end
    solver = 'd_pcg';
    dn = 10000;
    sp = sum(pp);
%%
    if (m <= dn) && Ayes  
       if (m <= 1000) %1000
          solver = 'd_direct';
       elseif sp <= max(0.01*par.n,dn)
          solver = 'd_direct';
       end
    end
    if sp <= 0.7*m && Ayes && sp <=dn 
       solver = 'p_direct';
    end
    if (m > 5e3 && sp >= 1.5e3) || (m>2000 && sp > 1.5e3) || (m > 100 && sp > 1e4)
       solver = 'd_pcg';
%        solver = 'p_direct';
    end
    
%%

    if strcmp(solver,'d_pcg')
       if Ayes
          %AP = Ainput.A(pp,:); AP = DD*AP;
          if false
             tmp = sum(AP.*AP,2);
             par.precond = 1;
             par.invdiagM = 1./(1 + par.sigma*tmp);
          end
          [xi,~,resnrm,solve_ok] = ...
          psqmry('matvec_ClassicLasso1',AP,rhs,par); 
       else
          [xi,~,resnrm,solve_ok] = ...
          psqmry('matvec_ClassicLasso_Amap1',opts,rhs,par); 
       end
    elseif strcmp(solver,'d_direct')
       %AP = Ainput.A(pp,:); AP = DD*AP;
       sigAPAt = (AP'*AP);
       if m <= 1500      
          M = par.mu * eye(m) + sigAPAt;
          xi = M\rhs;
       else
          M = par.mu * speye(m,m) + sigAPAt;  
          L = mychol(M,m);
          xi = mylinsysolve(L,rhs);
       end
       resnrm = 0; solve_ok = 1;
    elseif strcmp(solver,'p_direct')
       %AP = Ainput.A(pp,:); AP = DD*AP;
%        APT = AP';
       rhstmp = AP*rhs/mu; 
       PAtAP = AP*AP';
       if sp <= 1500
          M = eye(sp) + PAtAP/mu;
          tmp = M\rhstmp;
       else
          M = speye(sp,sp) + PAtAP;
          L = mychol(M,sp);
          tmp = mylinsysolve(L,rhstmp);
       end
       resnrm = 0; solve_ok = 1;
       xi = rhs/mu - AP'*tmp/mu;
    end
%%**********************************************************************
