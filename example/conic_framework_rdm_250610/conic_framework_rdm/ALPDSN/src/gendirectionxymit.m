
  
function [dx,dy,trans]= gendirectionxymit(Ftz,K,At,trans,cgopts, method, model_original)
    
b_cone_index = trans.b_cone_index;

m = size(At{1}, 2);
Amap  = trans.Amap;
ATmap = trans.ATmap;
sig = trans.sigma; %注意这里与SSNSDP的不同之处
tau1 = trans.NEWT.tau1;
tau2 = trans.NEWT.tau2;
CG_maxit = cgopts.CG_maxit;
pcgTol = cgopts.CG_tol;
Lchol = trans.Lchol;
Fx = Ftz.FX;
Fy = Ftz.FY;
rhsx = -Fx;
rhsy = -Fy;
%Ftz.par: \Sigma

% D = Ftz.par;

% Dtau2 = ProjJac_ops(D, K, 'affine', - 1 / sig, 1 / sig + tau2); % D^{\tau_2}

% invDtau2 = ProjJac_ops(Dtau2, K, 'affine_inv', 1 , 0); % (D^{\tau_2})^{-1}


tmptau = ( 1 + sig*tau2 );


% check Dtau2 * invDtau2 = I

% r = rhsx ;
% r1 = DPhi(K, Dtau2, r);
% r2 = - 1 / sig * (DPhi(K, D, r) + tmptau * r);
% norm(r1 - r2)
% norm(DPhi(K, invDtau2, r1))
% norm(DPhi(K, invDtau2, r1) - r)

% D * (D^{\tau_2})^{-1} = - sig * I -  sig * (1+ sig*tau_2) * (D + (1 + sig*tau_2) * I)^{-1}
% tmpinvD = ProjJac_ops(D, K, 'affine_inv', 1, tmptau); % tmpinvD = (D + (1 + sig*tau_2) * I)^{-1}
% DinvDtau2 = ProjJac_ops(tmpinvD, K, 'affine', sig * tmptau, - sig);


% iHW = DinvDtau2;

% original code 

% for p =1: length(K)
%     cone = K{p};
%     if strcmp(cone.type, 's')
%         iHW.Dsch2{p} = sig * D.Dsch2{p} ./ (1 + sig * tau2 - D.Dsch2{p}); %iHW: \tilde{\Sigma} 
%         iHW.Dsch2t{p} = iHW.Dsch2{p}';
%         iHW.Dsch1{p} = 1 / tau2 ;
%     elseif strcmp(cone.type, 'q')
%         iHW.Dsch1{p} = zeros(size(D.Dsch1{p}));
%         iHW.Dsch2{p} = zeros(size(D.Dsch2{p}));
%         iHW.shift{p} = zeros(size(D.shift{p}));
%         idx1 = (D.dd{p}(:, 2) > 0);   % both eigenvalues are positive
%         idx2 = (D.dd{p}(:, 1) >= 0 & D.dd{p}(:, 2) <= 0); % one eigenvalue is positive
%         idx3 = (D.dd{p}(:, 1) < 0); % both eigenvalues are negative
%         iHW.Dsch1{p}(idx1) = 0;
%         iHW.Dsch2{p}(idx1) = 0;
%         iHW.shift{p}(idx1) = 1 / tau2;
%         iHW.Dsch1{p}(idx2) = - sig * tmptau ./ (D.shift{p}(idx2) + tmptau) .* (D.Dsch1{p}(idx2) ./ (D.Dsch1{p}(idx2) + D.shift{p}(idx2) + tmptau) );
%         iHW.Dsch2{p}(idx2) = - sig * tmptau ./ (D.shift{p}(idx2) + tmptau) .* (D.Dsch2{p}(idx2) ./ (D.Dsch2{p}(idx2) + D.shift{p}(idx2) + tmptau) );
%         iHW.shift{p}(idx2) = - sig * D.shift{p}(idx2) ./ ( D.shift{p}(idx2) + tmptau) ;
%         iHW.Dsch1{p}(idx3) = 0;
%         iHW.Dsch2{p}(idx3) = 0;
%         iHW.shift{p}(idx3) = 0;
%     elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
%         iHW.Dsch2{p} = sig * D.Dsch2{p} ./ (1 + sig * tau2 - D.Dsch2{p});
%         iHW.Dsch2t{p} = iHW.Dsch2{p}';
%     end
% end

iHW = Ftz.par;

% for k = 1:length(iHW.Dsch12)
%     iHW.Dsch2{k} = (sig*iHW.Dsch2{k})./(1+tau2*sig-iHW.Dsch2{k});
%     iHW.Dsch2t{k} = iHW.Dsch2{k}';
%     iHW.Dsch1{k} = 1/tau2;
% end

sumlen = 0;

for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 's')
        for jj = 1:length(cone.size)
            iHW.Dsch2{sumlen+jj} = sig * iHW.Dsch2{sumlen+jj} ./ (1 + sig * tau2 - iHW.Dsch2{sumlen+jj});
            iHW.Dsch2t{sumlen+jj} = iHW.Dsch2{sumlen+jj}';
            iHW.Dsch1{sumlen+jj} = 1 / tau2;
        end
        sumlen = sumlen + length(cone.size);
    elseif strcmp(cone.type, 'q')
        if trans.socp_formula == 1 || trans.socp_formula == 2
            iHW.Dsch1{sumlen+1} = - iHW.Dsch1{sumlen+1} / sig;
            iHW.Dsch2{sumlen+1} = - iHW.Dsch2{sumlen+1} / sig;
            iHW.shift{p} = - iHW.shift{p} / sig + 1 / sig + tau2;
            % inverse
            iHW.Dsch1{sumlen+1} = - iHW.Dsch1{sumlen+1} ./ (iHW.shift{p} .* (iHW.shift{p} + iHW.Dsch1{sumlen+1}));
            iHW.Dsch2{sumlen+1} = - iHW.Dsch2{sumlen+1} ./ (iHW.shift{p} .* (iHW.shift{p} + iHW.Dsch2{sumlen+1}));
            iHW.shift{p} = 1 ./ iHW.shift{p};
            % multiply by D (i.e. Ftz.par)
            temp_Dsch1 = Ftz.par.Dsch1{sumlen+1} .* iHW.Dsch1{sumlen+1} + Ftz.par.shift{p} .* iHW.Dsch1{sumlen+1} + Ftz.par.Dsch1{sumlen+1} .* iHW.shift{p};
            temp_Dsch2 = Ftz.par.Dsch2{sumlen+1} .* iHW.Dsch2{sumlen+1} + Ftz.par.shift{p} .* iHW.Dsch2{sumlen+1} + Ftz.par.Dsch2{sumlen+1} .* iHW.shift{p};
            iHW.Dsch1{sumlen+1} = temp_Dsch1;
            iHW.Dsch2{sumlen+1} = temp_Dsch2;
            iHW.shift{p} = Ftz.par.shift{p} .* iHW.shift{p};
        end

        if trans.socp_formula == 2 || trans.socp_formula == 3
            cumsum_index = [0 cumsum(cone.size)];
            tip_index = 1 + cumsum_index(1:end-1);
            body_example_index = 1 + tip_index;
            Lambda = Ftz.par.shift3{p}; % won't be modified, can use reference in c++
            Lambda_tip = Lambda(tip_index); % won't be modified, can use reference in c++
            Lambda_body = Lambda(body_example_index); % won't be modified, can use reference in c++
            Lambda1 = tmptau - Lambda;
            Lambda1_tip = Lambda1(tip_index); % won't be modified, can use reference in c++
            Lambda1_body = Lambda1(body_example_index); % won't be modified, can use reference in c++

            vec_a = Ftz.par.Dsch3{p}; % won't be modified, can use reference in c++
            vec_gamma = zeros(length(cone.size), 1);
            for jj = 1:length(cone.size)
                u = Ftz.par.P3{p}(cumsum_index(jj)+1:cumsum_index(jj+1)); % won't be modified, can use reference in c++
                vec_gamma(jj) = u(1)^2 / Lambda1_tip(jj) + u(2:end)' * u(2:end) / Lambda1_body(jj);
            end
            vec_c = vec_a ./ (1 - vec_a .* vec_gamma);
        end

        sumlen = sumlen + 1;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b') || strcmp(cone.type, 'b2l')
        iHW.Dsch2{sumlen+1} = (sig*iHW.Dsch2{sumlen+1})./(1+tau2*sig-iHW.Dsch2{sumlen+1});
        iHW.Dsch2t{sumlen+1} = iHW.Dsch2{sumlen+1}';
        iHW.Dsch1{sumlen+1} = 1/tau2;
        sumlen = sumlen + 1;
    end
end

if trans.socp_formula == 3
    N24rhsx = DPhimit_socp_ddinv(K, iHW, rhsx, sig, Lambda1, vec_c);
else
    N24rhsx = DPhimit(K, iHW, rhsx);
end
N24rhsx = Amap(N24rhsx);             %N_2(N_4)^{-1} Fx 
rhs1 = rhsy - N24rhsx;

% iHWy = D * Dtau2^{-1} * D + sig * D = sig * tmptau * I - sig * tmptau^2 * (D + tmptau * I)^{-1};
% invDtmptau = ProjJac_ops(D, K, 'affine_inv', 1, tmptau); % (D + tmptau * I)^{-1}
% iHWy = ProjJac_ops(invDtmptau, K, 'affine', -sig*tmptau^2, sig*tmptau);
% 
% iHWy.sig = 1;
% iHWy.epsilon = tau1;

iHWy = Ftz.par;
iHWy.sig = 1;

% for k = 1:length(iHW.Dsch12)
%     iHWy.Dsch12{k} = (iHWy.Dsch2{k}.*(1+sig*tau2)*sig)./(1+sig*tau2-iHWy.Dsch2{k});
%     iHWy.Dsch12t{k} = iHWy.Dsch2{k}';
%     iHWy.Dsch1{k} = (1+sig*tau2)/tau2;
% end

sumlen = 0;

for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 's')
        for jj = 1:length(cone.size)
            iHWy.Dsch12{sumlen+jj} = (iHWy.Dsch2{sumlen+jj}.*(1+sig*tau2)*sig)./(1+sig*tau2-iHWy.Dsch2{sumlen+jj});
            iHWy.Dsch12t{sumlen+jj} = iHWy.Dsch2{sumlen+jj}';
            iHWy.Dsch1{sumlen+jj} = (1+sig*tau2)/tau2;
        end
        sumlen = sumlen + length(cone.size);
    elseif strcmp(cone.type, 'q')
        if trans.socp_formula == 1 || trans.socp_formula == 2
            iHWy.Dsch1{sumlen+1} = Ftz.par.Dsch1{sumlen+1} .* iHW.Dsch1{sumlen+1} + Ftz.par.shift{p} .* iHW.Dsch1{sumlen+1} + Ftz.par.Dsch1{sumlen+1} .* (iHW.shift{p} + sig);
            iHWy.Dsch2{sumlen+1} = Ftz.par.Dsch2{sumlen+1} .* iHW.Dsch2{sumlen+1} + Ftz.par.shift{p} .* iHW.Dsch2{sumlen+1} + Ftz.par.Dsch2{sumlen+1} .* (iHW.shift{p} + sig);
            iHWy.shift{p} = Ftz.par.shift{p} .* (iHW.shift{p} + sig);
        end

        if trans.socp_formula == 2 || trans.socp_formula == 3
            vec_c0 = Lambda_tip ./ Lambda1_tip;
            vec_c1 = Lambda_body ./ Lambda1_body; % c0 and c1 have no relation with c
            vec_a_plus_acgamma = vec_a .* (1 + vec_c .* vec_gamma);
            vec_a_plus_a2gamma_plus_a2cgamma2 = vec_a .* (1 + vec_gamma .* vec_a_plus_acgamma);
            % b_0 = cc_0^2 + 2a(1+c\gamma)c_0 + a + a^2\gamma + a^2c\gamma^2, b_1 = cc_0c_1 + a(1+c\gamma)(c_0 + c_1) + a + a^2\gamma + a^2c\gamma^2, b_2 = cc_1^2 + 2a(1+c\gamma)c_1 + a + a^2\gamma + a^2c\gamma^2
            vec_b0 = vec_c .* vec_c0.^2 + 2*vec_a_plus_acgamma .* vec_c0 + vec_a_plus_a2gamma_plus_a2cgamma2;
            vec_b1 = vec_c .* vec_c0 .* vec_c1 + vec_a_plus_acgamma .* (vec_c0 + vec_c1) + vec_a_plus_a2gamma_plus_a2cgamma2;
            vec_b2 = vec_c .*vec_c1.^2 + 2*vec_a_plus_acgamma .* vec_c1 + vec_a_plus_a2gamma_plus_a2cgamma2;
            vec_u0 = Ftz.par.P3{p}(tip_index); % won't be modified, can use reference in c++
            iHWy.shift3{p} = Lambda.^2 ./ Lambda1 + Lambda;
            iHWy.shift3{p}(tip_index) = iHWy.shift3{p}(tip_index) + (vec_b0 - vec_b1.^2 ./ vec_b2) .* vec_u0.^2;
            iHWy.shift3{p} = iHWy.shift3{p} * sig;
            iHWy.Dsch3{p} = sig * vec_b2;
            iHWy.P3{p} = Ftz.par.P3{p};
            iHWy.P3{p}(tip_index) = iHWy.P3{p}(tip_index) .* vec_b1 ./ vec_b2;
            ncols = length(cone.size);
            nrows = sum(cone.size);
            Prow_index = 1:nrows;
            Pcol_index = repelem(1:ncols, cone.size);
            Pvalue3 = iHWy.P3{p}(:);
            iHWy.Q3{p} = sparse(Prow_index, Pcol_index, Pvalue3, nrows, ncols);
        end

        sumlen = sumlen + 1;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b') || strcmp(cone.type, 'b2l')
        iHWy.Dsch12{sumlen+1} = (iHWy.Dsch2{sumlen+1}.*(1+sig*tau2)*sig)./(1+sig*tau2-iHWy.Dsch2{sumlen+1});
        iHWy.Dsch12t{sumlen+1} = iHWy.Dsch2{sumlen+1}';
        iHWy.Dsch1{sumlen+1} = (1+sig*tau2)/tau2;
        sumlen = sumlen + 1;
    end
end

if isfield(iHWy, 'Sigma1')
    for k = 1:length(iHWy.Sigma1)
    %     [i, j, s] = find(iHWy.Sigma{k});
    %     s = (s .* (1 + sig * tau2) * sig) ./ (1 + sig * tau2 - s);
    %     iHWy.Sigma{k} = sparse(i, j, s, size(iHWy.Sigma{k}, 1), size(iHWy.Sigma{k}, 2));   
    
        [i, j, s] = find(iHWy.Sigma1{k});
        s = (s .* (1 + sig * tau2) * sig) ./ (1 + sig * tau2 - s);
        iHWy.Sigma1{k} = sparse(i, j, s, size(iHWy.Sigma1{k}, 1), size(iHWy.Sigma1{k}, 2));
        [i, j, s] = find(iHWy.Sigma2{k});
        s = (s .* (1 + sig * tau2) * sig) ./ (1 + sig * tau2 - s);
        iHWy.Sigma2{k} = sparse(i, j, s, size(iHWy.Sigma2{k}, 1), size(iHWy.Sigma2{k}, 2));
    end
end
iHWy.epsilon = tau1;
iHWy.K = K;

L = struct();

if strcmp(method, 'direct')
    At_original = model_original.At;
    if ~isfield(trans.Lchol_original, 'RtR')
        trans.Lchol_original.RtR = trans.Lchol_original.Rt * trans.Lchol_original.R;
        trans.Lchol_original.RtR(trans.Lchol_original.perm, trans.Lchol_original.perm) = trans.Lchol_original.RtR;
        trans.Lchol_original.Rperm = sparse(size(trans.Lchol_original.R, 1), size(trans.Lchol_original.R, 2));
        trans.Lchol_original.Rperm(:, trans.Lchol_original.perm) = trans.Lchol_original.R;
    end
    if ~isfield(trans, 'sp_info')
        trans.sp_info = detect_sp_info(At_original, K, trans.Lchol_original, trans.system_opt, trans.socp_formula, trans.sys2_sparse_strategy, 0.1);
        trans.standard_nnz = nnz(trans.sp_info.lhs_structure);
    end
    if strcmp(trans.direct_solver, 'pardiso') && ~isfield(trans, 'pardiso_info')
        dll_path = '../pardiso_solver';
        addpath(dll_path);
        % trans.pardiso_info = pardisoinit(-2, 0);
        % trans.pardiso_info = pardisoreorder(tril(trans.sp_info.lhs_structure), trans.pardiso_info, false);
        trans.pardiso_info = 1;
    end
end

if (b_cone_index == 0)
    Lchol_original = trans.Lchol_original;
    [dy,relres,flag] = solve_AHAt(method, trans, @matvec_y2mit, K, At, rhs1, iHWy, L, pcgTol, CG_maxit, Lchol_original);
else
    Lchol_original = trans.Lchol_original;
    t_constant = full(At{b_cone_index}(end, end));

    At_b = model_original.At{b_cone_index};
    u = rhs1;
    b_size = K{b_cone_index}.size / 2;
    u1 = u(1:end - b_size);
    u2 = u(end - b_size + 1:end);

    iHWy_temp = iHWy;
    d_full = iHWy_temp.Dsch12{b_cone_index};
    d_b_hat = d_full(1:b_size);
    d_2_hat = d_full(b_size + 1:end);
    d_sum = d_b_hat + d_2_hat + tau1 / (t_constant^2) * ones(b_size, 1);
    u2_temp = d_b_hat ./ d_sum .* u2 .* (1 / t_constant);
    r1 = u1 - fwsolve(Lchol_original, At_b' * u2_temp);
    % r1 = u1 - At_b' * u2_temp;

    %% calculate v1
    iHWy_temp.Dsch12{b_cone_index} = d_b_hat - d_b_hat ./ d_sum .* d_b_hat;
    iHWy_temp.Dsch12t{b_cone_index} = iHWy_temp.Dsch12{b_cone_index}';

    [v1, relres, flag] = solve_AHAt(method, trans, @matvec_y2mit,model_original.K,model_original.At,r1,iHWy_temp,L,pcgTol,CG_maxit,Lchol_original);

    r2 = u2 - t_constant * d_b_hat .* At_b * bwsolve(Lchol_original, v1);
    % r2 = u2 - t_constant * d_b_hat .* (At_b * v1);
    v2 = r2 ./ d_sum .* (1 / t_constant^2);
    dy = [v1; v2];
end

if length(relres) > 1
    trans.cgres = relres;
else
    trans.cgres = [relres, relres];
end
%%


trans.flag = flag;
trans.cgiter = length(relres);

rhs2 = ATmap(dy);
rhs2 = DPhimit(K, Ftz.par, rhs2);

rhs2 = rhs2 + rhsx ;
rhs2 = rhs2 / tmptau;
% iHWx = ProjJac_ops(D, K, 'affine_inv', - 1 / (sig * tmptau), - 1 / sig); % - (tmptau * D^{\tau_2})^{-1}

% dx= DPhimit(K, iHWx, rhs2) ;

iHWx = Ftz.par;
% for k = 1:length(iHW.Dsch12)
%     iHWx.Dsch2{k} = (iHWx.Dsch2{k}*sig)./(1+sig*tau2-iHWx.Dsch2{k});
%     iHWx.Dsch2t{k} = iHWx.Dsch2{k}';
%     iHWx.Dsch1{k} = 1/tau2;
% end

sumlen = 0;

for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 's')
        for jj = 1:length(cone.size)
            iHWx.Dsch2{sumlen+jj} = (iHWx.Dsch2{sumlen+jj}*sig)./(1+sig*tau2-iHWx.Dsch2{sumlen+jj});
            iHWx.Dsch2t{sumlen+jj} = iHWx.Dsch2{sumlen+jj}';
            iHWx.Dsch1{sumlen+jj} = 1/tau2;
        end
        sumlen = sumlen + length(cone.size);
    elseif strcmp(cone.type, 'q')
        if trans.socp_formula == 1 || trans.socp_formula == 2
            % scale and add identity
            iHWx.Dsch1{sumlen+1} = - iHWx.Dsch1{sumlen+1} / sig;
            iHWx.Dsch2{sumlen+1} = - iHWx.Dsch2{sumlen+1} / sig;
            iHWx.shift{p} = - iHWx.shift{p} / sig + 1 / sig + tau2;
            % inverse
            iHWx.Dsch1{sumlen+1} = - iHWx.Dsch1{sumlen+1} ./ (iHWx.shift{p} .* (iHWx.shift{p} + iHWx.Dsch1{sumlen+1}));
            iHWx.Dsch2{sumlen+1} = - iHWx.Dsch2{sumlen+1} ./ (iHWx.shift{p} .* (iHWx.shift{p} + iHWx.Dsch2{sumlen+1}));
            iHWx.shift{p} = 1 ./ iHWx.shift{p};
            % multiply by tmptau and subtract by sigma I
            iHWx.Dsch1{sumlen+1} = iHWx.Dsch1{sumlen+1} * tmptau;
            iHWx.Dsch2{sumlen+1} = iHWx.Dsch2{sumlen+1} * tmptau;
            iHWx.shift{p} = iHWx.shift{p} * tmptau - sig;
        elseif trans.socp_formula == 3
            iHWx.Dsch3{p} = sig * tmptau * vec_c;
            iHWx.shift3{p} = (tmptau ./ Lambda1 - 1) * sig;
            iHWx.P3{p} = Ftz.par.P3{p} ./ Lambda1;
        end
        sumlen = sumlen + 1;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b') || strcmp(cone.type, 'b2l')
        iHWx.Dsch2{sumlen+1} = (iHWx.Dsch2{sumlen+1}*sig)./(1+sig*tau2-iHWx.Dsch2{sumlen+1});
        iHWx.Dsch2t{sumlen+1} = iHWx.Dsch2{sumlen+1}';
        iHWx.Dsch1{sumlen+1} = 1/tau2;
        sumlen = sumlen + 1;
    end
end


%     dx= DPhi_sdpnalN4T(blk,iHWx,rhs2,sig,tau);
dx= DPhimit(K,iHWx,rhs2);
for k = 1:trans.nblock
    dx{k} = dx{k} + sig * rhs2{k};
end

% for k = 1:trans.nblock
%     dx{k} = dx{k} + sig * rhs2{k};
% end
% % original code
% rhs2 = rhs2 - Fx ;
% rhs2 = rhs2 / (1 + sig * tau2);
% % iHWx: T
% iHWx = D;
% for p =1: length(K)
%     cone = K{p};
%     if strcmp(cone.type, 's')
%         iHWx.Dsch2{p} = (sig * D.Dsch2{p} )./ (1 + sig * tau2 - D.Dsch2{p} );
%         iHWx.Dsch1{p} = 1 / tau2;
%         iHWx.shift{p} = sig;
%     end
% end
% dx_= DPhi(K, iHWx, rhs2);

% norm(dx - dx_);



%% check residual of newton system wrt [dx; dy]
D_lmut = @(x) DPhimit(K, Ftz.par, x);
Atdy = ATmap(dy);
N1dy = sig * Amap(D_lmut(Atdy)) + tau1 * dy;
N2dx = Amap(D_lmut(dx));
N3dy = - D_lmut(Atdy) ;
N4dx = (1 / sig + tau2) * dx - 1 / sig * D_lmut(dx);
resy = N1dy + N2dx - rhsy;
resx = N3dy + N4dx - rhsx;
trans.newton_res = sqrt(norm(resy) ^ 2 + norm(resx) ^ 2) / (1 + sqrt(norm(rhsy) ^ 2 + norm(rhsx) ^ 2)) ;
% 
% 
% trans.newton_res = -1;


end

function [dy, relres, flag] = solve_AHAt(method, trans,matvecfname,K,At,rhs1,iHWy,L,pcgTol,CG_maxit,Lchol)
    if strcmp(method, 'iterative')
        [dy,relres,flag] = psqmr2(matvecfname,K,At,rhs1,iHWy,L,pcgTol,CG_maxit,Lchol);

    elseif strcmp(method, 'direct')
        % define default system_opt
        if ~isfield(trans, 'system_opt')
            trans.system_opt = 2;
        end

        if trans.system_opt == 0
            % 0: a basic implementation, support 's' cone, recommended for very sparse A
            [dy, relres, flag] = solve_direct(K, At, rhs1, iHWy, Lchol, trans);
        elseif trans.system_opt == 2
            % 2: for SOCP only, but more efficient, handle dense&sparse rows&columns
            [dy, relres, flag] = solve_direct_opt2(K, At, rhs1, iHWy, Lchol, trans);
        else
            error('Unknown system_opt');
        end
    end
end