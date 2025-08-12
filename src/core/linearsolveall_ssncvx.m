
function  [out,resnrm,solve_ok] = linearsolveall_ssncvx(matvecfname,rhs,par,L,tol,maxit,params)
sub_maxit = 100;


if strcmp(params.method,'iterative')
    if params.nblock == 1
        y0 =  params.y0;
        zb0 = params.z0{1};
        r0 =  params.r0{1};
        v0 =  params.v0{1};
        AL = params.Lchol;
        solve_ok = 1;
        miniter = 1;
        printlevel = 0;
        if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end
        if isfield(par,'printlevel'); printlevel = par.printlevel; end

        %%
        z1 = y0;
        z2 = zb0;
        z3 = r0;
        z4 = v0;
        r1 = rhs.rhsy;
        r2 = rhs.rhsz{1};
        r3 = rhs.rhsr{1};
        r4 = rhs.rhsv{1};

        %%
        q.q1 = rhs.rhsy;
        q.q2 = rhs.rhsz ;
        q.q3 = rhs.rhsr;
        q.q4 = rhs.rhsv;
        tau_old = sqrt(norm(q.q1)^2 + norm(q.q2)^2 + norm(q.q3)^2 + norm(q.q4)^2 );
        rho_old  = dot_ssn(r1,q.q1) + dot_ssn(r2,q.q2) + dot_ssn(r3,q.q3) + dot_ssn(r4,q.q4);
        theta_old = 0;
        d1 = y0;
        d2 = zb0;
        d3 = r0;
        d4 = v0;
        res1 = r1; res2 = r2; res3 = r3; res4 = r4;

        err = norm(r1) ; resnrm(1) = err; minres = err;
        Ad1 = y0;
        Ad2 = zb0;
        Ad3 = r0;
        Ad4 = v0;

        tiny = 1e-30;

        for iter = 1:maxit
            Aq = feval(matvecfname,params.At,par,q,params);

            sigma = dot_ssn(q.q1,Aq.Ax1) + dot_ssn(q.q2,Aq.Ax2) + dot_ssn(q.q3,Aq.Ax3) + dot_ssn(q.q4,Aq.Ax4);
            if (abs(sigma) < tiny)
                solve_ok = 2;
                if (printlevel); fprintf('s1'); end
                break;
            else
                alpha = rho_old/sigma;
                r1 = r1 - alpha*Aq.Ax1;
                r2 = -alpha*Aq.Ax2 + r2;
                r3 = -alpha*Aq.Ax3 + r3;
                r4 = -alpha*Aq.Ax4 + r4;
            end
            u1 = r1;
            u2 = r2;
            u3 = r3;
            u4 = r4;
            %%
            % theta = sqrt((norm(u1)^2 + norm(u2)^2  + norm(u3)^2 + norm(u4)^2 ))/tau_old;
            theta = sqrt(u1(:)'*u1(:) + u2(:)'*u2(:) + u3(:)'*u3(:) + u4(:)'*u4(:))/tau_old;
            c = 1/sqrt(1+theta^2);
            tau = tau_old*theta*c;
            gam = (c^2*theta_old^2);
            eta = (c^2*alpha);
            d1 = gam*d1 + eta*q.q1;
            d2 = gam*d2 + eta*q.q2{1};
            d3 = gam*d3 + eta*q.q3{1};
            d4 = gam*d4 + eta*q.q4{1};


            z1 = z1 + d1;
            z2 = z2 + d2;
            z3 = z3 + d3;
            z4 = z4 + d4;

            %%----- stopping conditions ----
            Ad1 = gam*Ad1 + eta*Aq.Ax1;
            Ad2 = gam*Ad2 + eta*Aq.Ax2;
            Ad3 = gam*Ad3 + eta*Aq.Ax3;
            Ad4 = gam*Ad4 + eta*Aq.Ax4;

            res1 = res1 - Ad1;
            res2 = res2 - Ad2;
            res3 = res3 - Ad3;
            res4 = res4 - Ad4;
            err = sqrt(res1(:)'*res1(:) + res2(:)'*res2(:) + res3(:)'*res3(:) + res4(:)'*res4(:) );
            resnrm(iter+1) = err;
            if (err < minres); 
                minres = err;
            end
            if (err < tol) && (iter > miniter)
                break;
            end

            if (abs(rho_old) < tiny)
                solve_ok = 2;
                fprintf('s2');
                break;
            else
                rho = r1(:)'*r1(:) + r2(:)'*r2(:) + r3(:)'*r3(:) + r4(:)'*r4(:);
                beta = rho/rho_old;
                q.q1 = r1 + beta*q.q1;
                q.q2{1} = r2 + beta*q.q2{1};
                q.q3{1} = r3 + beta*q.q3{1};
                q.q4{1} = r4 + beta*q.q4{1};
            end
            rho_old = rho;
            tau_old = tau;
            theta_old = theta;
        end
        z.q1 = z1;
        z.q2 = z2;
        z.q3 = z3;
        z.q4 = z4;

        out.dy = z1;
        out.dz{1} = z2;
        out.dr{1} = z3;
        out.dv{1} = z4;
        if (iter == maxit); solve_ok = -2; end
        if (solve_ok ~= -1)
            if (printlevel); fprintf(' '); end
        end
    else
        y0 =  params.y0;
        zb0 = params.z0;
        r0 =  params.r0;
        v0 =  params.v0;
        AL = params.Lchol;
        solve_ok = 1;
        miniter = 1;
        printlevel = 0;
        if isfield(par,'minitpsqmr'); miniter = par.minitpsqmr; end
        if isfield(par,'printlevel'); printlevel = par.printlevel; end

        %%
        z1 = y0;
        z2 = zb0;
        z3 = r0;
        z4 = v0;
        r1 = rhs.rhsy;
        r2 = rhs.rhsz;
        r3 = rhs.rhsr;
        r4 = rhs.rhsv;

        %%
        q.q1 = rhs.rhsy;
        q.q2 = rhs.rhsz ;
        q.q3 = rhs.rhsr;
        q.q4 = rhs.rhsv;
        tau_old = sqrt(norm(q.q1)^2 + norm(q.q2)^2 + norm(q.q3)^2 + norm(q.q4)^2 );
        rho_old  = dot_ssn(r1,q.q1) + dot_ssn(r2,q.q2) + dot_ssn(r3,q.q3) + dot_ssn(r4,q.q4);
        theta_old = 0;
        d1 = y0;
        d2 = zb0;
        d3 = r0;
        d4 = v0;
        res1 = r1; res2 = r2; res3 = r3; res4 = r4;

        err = norm(r1) ; resnrm(1) = err; minres = err;
        Ad1 = y0;
        Ad2 = zb0;
        Ad3 = r0;
        Ad4 = v0;
        %%
        %% main loop
        %%
        tiny = 1e-30;

        for iter = 1:maxit
            Aq = feval(matvecfname,params.At,par,q,params);
            sigma = dot_ssn(q.q1,Aq.Ax1) + dot_ssn(q.q2,Aq.Ax2) + dot_ssn(q.q3,Aq.Ax3) + dot_ssn(q.q4,Aq.Ax4);
            if (abs(sigma) < tiny)
                solve_ok = 2;
                if (printlevel); fprintf('s1'); end
                break;
            else
                alpha = rho_old/sigma;
                r1 = r1 - alpha*Aq.Ax1;
                r2 = -alpha*Aq.Ax2 + r2;
                r3 = -alpha*Aq.Ax3 + r3;
                r4 = -alpha*Aq.Ax4 + r4;
            end
            u1 = r1;
            u2 = r2;
            u3 = r3;
            u4 = r4;
            %%
            theta = sqrt((norm(u1)^2 + norm(u2)^2  + norm(u3)^2 + norm(u4)^2 ))/tau_old;
            % theta
            c = 1/sqrt(1+theta^2);
            tau = tau_old*theta*c;
            gam = (c^2*theta_old^2);
            eta = (c^2*alpha);
            d1 = gam*d1 + eta*q.q1;
            d2 = gam*d2 + eta*q.q2;
            d3 = gam*d3 + eta*q.q3;
            d4 = gam*d4 + eta*q.q4;


            z1 = z1 + d1;
            z2 = z2 + d2;
            z3 = z3 + d3;
            z4 = z4 + d4;

            %%----- stopping conditions ----
            Ad1 = gam*Ad1 + eta*Aq.Ax1;
            Ad2 = gam*Ad2 + eta*Aq.Ax2;
            Ad3 = gam*Ad3 + eta*Aq.Ax3;
            Ad4 = gam*Ad4 + eta*Aq.Ax4;

            res1 = res1 - Ad1;
            res2 = res2 - Ad2;
            res3 = res3 - Ad3;
            res4 = res4 - Ad4;
            err = sqrt(norm(res1)^2 + norm(res2)^2 + norm(res3)^2 + norm(res4)^2 );
            % err = sqrt(norm(res1)^2 + norm(res3)^2);
            resnrm(iter+1) = err;
            if (err < minres); minres = err; end
            if (err < tol) && (iter > miniter)
                break;
            end

            if (abs(rho_old) < tiny)
                solve_ok = 2;
                fprintf('s2');
                break;
            else
                rho = norm(r1)^2 + norm(r2)^2 + norm(r3)^2 + norm(r4)^2;
                beta = rho/rho_old;
                q.q1 = r1 + beta*q.q1;
                q.q2 = r2 + beta*q.q2;
                q.q3 = r3 + beta*q.q3;
                q.q4 = r4 + beta*q.q4;
            end
            rho_old = rho;
            tau_old = tau;
            theta_old = theta;
        end
        z.q1 = z1;
        z.q2 = z2;
        z.q3 = z3;
        z.q4 = z4;

        out.dy = z1;
        out.dz = z2;
        out.dr = z3;
        out.dv = z4;
        if (iter == maxit); solve_ok = -2; end
        if (solve_ok ~= -1)
            if (printlevel); fprintf(' '); end
        end

    end
elseif strcmp(params.method,'direct')
    use_socp_solver = 1;
    for p = 1:length(params.pblk)
        if ~(strcmp(params.pblk{p}.type, 'l') || strcmp(params.pblk{p}.type, 'q') || strcmp(params.pblk{p}.type, 'u') || strcmp(params.pblk{p}.type, 'b2l'))
            use_socp_solver = 0;
            break;
        end
    end
    if use_socp_solver
        opts.K = params.pblk;
        opts.At = params.At;
        rhsz = rhs.rhsy;
        optcg.epsilon = par.epsilon1;
        optcg.sig = par.sigma;
        for p = 1:length(params.pblk)
            if strcmp(params.pblk{p}.type, 'l') || strcmp(params.pblk{p}.type, 'u') || strcmp(params.pblk{p}.type, 'b2l')
                optcg.Dsch12{p} = par.D11{p};
            elseif strcmp(params.pblk{p}.type, 'q')
                % optcg.Dsch1{p} = par.qD1;
                % optcg.Dsch2{p} = par.qD2;
                optcg.Dsch1{p} = par.D11{p};
                optcg.Dsch2{p} = par.D12{p};
                optcg.Q1{p} = par.D4{p}.Q1;
                optcg.Q2{p} = par.D4{p}.Q2;
                optcg.P1{p} = par.D4{p}.P1;
                optcg.P2{p} = par.D4{p}.P2;
                optcg.shift{p} = par.qshift;
                % elseif strcmp(params.pblk{p}.type, 'u')
                %     optcg.Dsch12{p} = par.uD;
                % elseif strcmp(params.pblk{p}.type, 'b2l')
                %     optcg.Dsch12{p} = par.b2lD;
            end
        end

        [dy,resnrm,solve_ok] = socp_ldl_solver(opts,rhsz,optcg);

        out.dy = dy;
        for p = 1:length(params.pblk)
            out.dz{p,1} = [];
            out.dr{p,1} = [];
            out.dv{p,1} = [];
        end
    elseif params.nblock == 1 && (strcmp(params.pblk{1}.type,'l1') || strcmp(params.pblk{1}.type,'l1con') || strcmp(params.pblk{1}.type,'linftycon'))
        optcg.rr = ~params.D.D4{1};
        if isfield(params.D,'D2')
            optcg.rr =  ~ logical(params.D.D2{1} + params.D.D4{1});
        end

        optcg.mu = 1+par.epsilon2;
        optcg.DD = sqrt( par.tDD);
        if isfield(params.D,'D2')
            optcg.DD = sqrt( par.tDD2 );
        end
        optcg.n = length(optcg.rr);
        optcg.Ayes = params.Ayes;
        optcg.sub_maxit = sub_maxit;
        rhsz = rhs.rhsz{1};

        opts.Ayes = params.Ayes;
        opts.A = params.B;
        tmp = params.D.D4{1};
        if isfield(params.D,'D2')
            AP = opts.A(logical(params.D.D4{1} + params.D.D2{1}  ),:);
        else
            AP = opts.A(:,tmp)';
        end
        tolpsqmr = min(5e-3, 0.05*norm(rhsz));

        const2 = 1;
        % if times > 1 && prim_ratio > 0.5
        %     const2 = 0.5*const2;
        % end
        % if (dual_ratio > 1.1); const2 = 0.5*const2; end
        tolpsqmr = const2*tolpsqmr;
        optcg.tol = tolpsqmr;
        optcg.AP = AP;
        [dz,resnrm,solve_ok,solver] = Classic_Lasso_linsys_solver1(opts,rhsz,optcg);

        out.dy = {[]};
        out.dz = {dz};
        out.dr = {[]};
        out.dv = {[]};
    elseif params.nblock == 1 && strcmp(params.pblk{1}.type,'fused')
        optcg.AP = par.AP;
        optcg.Ph = par.Ph;
        optcg.PU = par.PU;
        optcg.PU1 = par.PU1;
        optcg.V1 = par.V1;
        optcg.mu1 = 1+par.epsilon2;
        optcg.PU2 = par.PU2;
        optcg.numblk1 = par.numblk1;
        optcg.lenP = par.lenP;
        optcg.info_u = par.info_u;
        rhsz = rhs.rhsz{1};
        opts.Ayes = params.Ayes;
        opts.A = params.B;
        [dz,resnrm,solve_ok,solver] =  Fused_Lasso_linsys_solver1(opts,rhsz,optcg);;

        out.dy = {[]};
        out.dz = {dz};
        out.dr = {[]};
        out.dv = {[]};
    end
end



