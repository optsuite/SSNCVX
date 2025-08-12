function [dy, relres, flag] = solve_direct(K, At, rhs1, iHWy, Lchol, trans)
    %% use direct method to solve linear system
    % lhs* dy = rhs1
    % lhs = iHWy.sig * Rt^{-1} * S * A * iHWy.Dsch * At * Sinv * R^{-1} + epsilon * I
    % Sinv * y is equivalent to y(perm) = y
    % S * y is equivalent to y = y(perm)
    % lhs * dy = rhs is equivalent to 
    % (iHWy.sig  * A * iHWy.Dsch * At  + epsilon * Sinv * Rt * R * S) * (Sinv * R^{-1} * dy) = Sinv * Rt * rhs
    % compute RtR = Sinv * Rt * R * S 
    % if ~isfield(Lchol, 'RtR')
    %     Lchol.RtR = Lchol.Rt * Lchol.R;
    %     Lchol.RtR(Lchol.perm, Lchol.perm) = Lchol.RtR;
    % end

    direct_solver = trans.direct_solver;

    % compute lhs = iHWy.sig * A * iHWy.Dsch * At  + epsilon * Sinv * Rt * R * S
    m = size(At{1}, 2);
    par = iHWy; 
    AHAt = struct();
    AHAt.mat11 = sparse(m, m);

    sumlen = 0;

    for p =1: length(K)
        cone = K{p};
        if strcmp(cone.type, 's')
            if ~isempty(par.Sigma1{p})
                pblk = {K{p}.type,K{p}.size};
                tmp1 = mexskron(pblk,par.Q1{p}',par.Q1{p}');
                tmp2 = mexskron(pblk,par.Q1{p}',par.Q2{p}');
                schurtmp1 = tmp1*At{p,1}.*sqrt(svecd(pblk, par.Sigma1{p} ));
                schurtmp2 = 2*tmp2*At{p,1}.*sqrt(svecd(pblk, par.Sigma2{p} + par.Sigma2{p}'));
                AHAt.mat11 = AHAt.mat11 + schurtmp1'*schurtmp1 + schurtmp2'*schurtmp2;     
            end
        elseif strcmp(cone.type, 'q')
            if (~isempty(par.Dsch2{sumlen+1}))
                temp1 = At{p}' * par.Q1{p};
                temp1 = temp1 .* par.Dsch1{sumlen+1}' * temp1';
                temp2 = At{p}' * par.Q2{p};
                temp2 = temp2 .* par.Dsch2{sumlen+1}' * temp2';
                temp3 = At{p}' .* repelem(par.shift{p}, cone.size, 1)' * At{p};
                AHAt.mat11 = AHAt.mat11 + temp3 + temp1 + temp2;
             end    
        elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'b') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
            if (~isempty(par.Dsch12{p}))
                AHAt.mat11 = AHAt.mat11 + At{p}' .* par.Dsch12{sumlen+1}' * At{p};
            end
        elseif strcmp(cone.type, 'u')
        
        end
        sumlen = sumlen + strcmp(cone.type, 's') * length(cone.size) + (1 - strcmp(cone.type, 's')) * 1;
    end

    AHAt.mat11 = par.sig * AHAt.mat11;
    AHAt.mat11 = AHAt.mat11 + par.epsilon * Lchol.RtR;
    % AHAt.mat11 = AHAt.mat11 + par.epsilon * speye(m);
    
    % compute rhs1_temp = Sinv * Rt * rhs1
    rhs1_temp = zeros(size(rhs1));
    rhs1_temp(Lchol.perm) = Lchol.Rt * rhs1;
    %rhs1_temp = rhs1;

    % solve lhs * ( Sinv * R^{-1} * dy) = rhs
    % dy = AHAt.mat \ rhs1_temp;

    if strcmp(direct_solver, 'backslash')
        dy = AHAt.mat11 \ rhs1_temp;
    elseif strcmp(direct_solver, 'ldl')
        [L_ldl, D_ldl, P_ldl] = ldl(AHAt.mat11, 'vector'); % 'vector' for permutation vector
        bp_ldl = rhs1_temp(P_ldl, :);
        y_ldl = L_ldl \ bp_ldl;
        z_ldl = D_ldl \ y_ldl;
        x_permuted = L_ldl' \ z_ldl;
        dy = zeros(size(rhs1_temp, 1), size(rhs1_temp, 2));
        dy(P_ldl, :) = x_permuted;
    elseif strcmp(direct_solver, 'pardiso')
        dy = solve_pardiso(AHAt.mat11, rhs1_temp, trans.pardiso_info);
    end

    % recover dy 
    dy = Lchol.R * dy(Lchol.perm);

    % check residual
    relres = norm(rhs1 - matvec_y2mit(K, At, iHWy, dy, Lchol) ) / (1 + norm(rhs1));
    
    % fprintf("resdue of computing dy: %e\n", relres);
    flag = 1;
end
