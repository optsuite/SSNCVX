%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-18 19:42:12
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-18 19:59:34
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%

function [dy,dzb,dx,dq,trans]= gendirectionxz(Ftz,K,At,trans,cgopts)

Amap  = trans.Amap;
ATmap = trans.ATmap;
sig = trans.sigma;
tau1 = trans.NEWT.tau1;
tau2 = trans.NEWT.tau2;
tau3 = trans.NEWT.tau3;
tau4 = trans.NEWT.tau4;
CG_maxit = cgopts.CG_maxit;
pcgTol = cgopts.CG_tol;
Lchol = trans.Lchol;
Fx = Ftz.FX;
Fy = Ftz.FY;
FZ = Ftz.FZ;
Fq = Ftz.Fq;

D = Ftz.par;
iHW = D;
tmptau = - (1 + sig*tau3);
for p =1: length(K)
    cone = K{p};
    if strcmp(cone.type, 's')
        iHW.Dsch2{p} = sig * D.Dsch2{p} ./ (1 + sig * tau3 - D.Dsch2{p}); %iHW: \tilde{\Sigma}
        iHW.Dsch1{p} = 1 / tau3 ;
    elseif strcmp(cone.type, 'q')
        iHW.Dsch1{p} = zeros(size(D.Dsch1{p}));
        iHW.Dsch2{p} = zeros(size(D.Dsch2{p}));
        iHW.shift{p} = zeros(size(D.shift{p}));
        idx1 = (D.dd{p}(:, 2) > 0);   % both eigenvalues are positive
        idx2 = (D.dd{p}(:, 1) >= 0 & D.dd{p}(:, 2) <= 0); % one eigenvalue is positive
        idx3 = (D.dd{p}(:, 1) < 0); % both eigenvalues are negative
        iHW.Dsch1{p}(idx1) = 0;
        iHW.Dsch2{p}(idx1) = 0;
        iHW.shift{p}(idx1) = 1 / tau3;
        iHW.Dsch1{p}(idx2) = - sig * tmptau ./ (D.shift{p}(idx2) + tmptau) .* (D.Dsch1{p}(idx2) ./ (D.Dsch1{p}(idx2) + D.shift{p}(idx2) + tmptau) );
        iHW.Dsch2{p}(idx2) = - sig * tmptau ./ (D.shift{p}(idx2) + tmptau) .* (D.Dsch2{p}(idx2) ./ (D.Dsch2{p}(idx2) + D.shift{p}(idx2) + tmptau) );
        iHW.shift{p}(idx2) = - sig * D.shift{p}(idx2) ./ ( D.shift{p}(idx2) + tmptau) ;
        iHW.Dsch1{p}(idx3) = 0;
        iHW.Dsch2{p}(idx3) = 0;
        iHW.shift{p}(idx3) = 0;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
        iHW.Dsch2{p} = sig * D.Dsch2{p} ./ (1 + sig * tau3 - D.Dsch2{p});
    end
end

tmp1 = DPhi(K, iHW, Fx);
N24rhsx = Amap(tmp1);
rhs1 = -Fy + N24rhsx;

% invDtmptau = ProjJac_ops(D, K, 'affine_inv', 1, tmptau); % (D + tmptau * I)^{-1}
% iHWy = ProjJac_ops(invDtmptau, K, 'affine', -sig*tmptau^2, sig*tmptau);
%
% tq = Fq{1};
% tq(Ftz.Dh{1}==0) = 0;
% tmp2 = 1 / tau4 * tq;
% rsh2 = tmp1 - tmp2 - FZ;
%
% iHWy.sig = 1;

% % iHWy.ttau = sig2 + 1/tau4;
% iHWy.Dh4 = (sig + 1/tau4)*Ftz.Dh{1};
tq = Fq;
for p = 1:length(K)
tq{p}(Ftz.Dh{p}==0) = 0;
end
tmp2 = tq/tau4;
rsh2 = tmp1-tmp2-FZ;


iHWy = Ftz.par;
iHWy.sig = 1;
for k = 1:trans.nblock
    iHWy.Dsch2{k} = (iHWy.Dsch2{k}.*(1+sig*tau3)*sig)./(1+sig*tau3-iHWy.Dsch2{k});
    iHWy.Dsch1{k} = (1+sig*tau3)/tau3;
end
iHWy.epsilon1 = tau1;
iHWy.epsilon2 = tau2;
iHWy.Dh4 = (sig + 1/tau4)*Ftz.Dh{1};

L = struct();

[dy,dzb,relres,flag] = psqmrplus(@matvec_yz,K,At,rhs1,rsh2,iHWy,L, pcgTol, CG_maxit, Lchol);

trans.cgres = [trans.cgres relres(end)];
trans.flag = flag;
trans.cgiter = length(relres);

dx = ATmap(dy) + dzb;
dx = DPhi(K,D,dx);
dx =  -(dx - Fx)/tmptau;

iHWx = D;
for p = 1:length(K)
    iHWx.Dsch2{p} = sig * D.Dsch2{p} ./ (1 + sig * tau3 - D.Dsch2{p}); %iHW: \tilde{\Sigma}
    iHWx.Dsch1{p} = 1 / tau3 ;
end
%     dx= DPhi_sdpnalN4T(blk,iHWx,rhs2,sig,tau);

tmp2= DPhi(K,iHWx,dx);
for k = 1:length(K)
    dx{k} = tmp2{k} + sig * dx{k};
end

% iHWx = ProjJac_ops(D, K, 'affine_inv', - 1 / (sig * tmptau), - 1 / sig); % - (tmptau * D^{\tau_2})^{-1}
%
% %     dx= DPhi_sdpnalN4T(K,iHWx,rhs2,sig,tau);
% dx= DPhi(K,iHWx,dx);
% dx = tmp2 + sig * dx;

tzb = Ftz.Dh .* dzb;
tmp3 = - (tzb + Fq);
dq = (tmp3*sig/(1+tau4*sig) + (Ftz.Dh .* tmp3)*1/(tau4*(1+tau4*sig)));
%%

D_lmut = @(x) DPhi(K, Ftz.par, x);
Atdy = ATmap(dy);
dzbtmp{1} = dzb; 
tmpz = DPhi(K, Ftz.par, dzbtmp);
tmpx = DPhi(K, Ftz.par, dx);
N11dy = sig * Amap(D_lmut(Atdy)) + tau1 * dy;
N12dz = sig* Amap(tmpz);
N13dx = Amap(tmpx);

N21dy = sig*D_lmut(Atdy); 
N22dz = sig*tmpz + sig*Ftz.Dh.*dzb + tau1 * dzb;
N23dx = tmpx;
N24dq = -Ftz.Dh.* dzb;

N31dy = N21dy/-sig;
N33dx = (1 / sig + tau3) * dx - 1 / sig * D_lmut(dx);

% N31dy =  N21dx/-sig ;
N44dq = (1 / sig + tau4) * dq - 1 / sig * Ftz.Dh.*dq;
resy = N11dy + N12dz + N13dx + Fy;
resz = N21dy + N22dz + N23dx + N24dq + FZ;
resx = N31dy - N23dx + N33dx + Fx;
resq = -N24dq + N44dq + Fq;
trans.newton_res = sqrt(norm(resy) ^ 2 + norm(resz) ^ 2 + norm(resx) ^ 2 + norm(resq) ^ 2) / (1 + sqrt(norm(Fy) ^ 2 + norm(Fx) ^ 2 + norm(FZ) ^ 2 +norm(Fq) ^ 2 )) ;
end