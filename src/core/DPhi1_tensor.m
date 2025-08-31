%% Dphi1_tensor: use the iterative method to solve the linear equation
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
function iHWu= DPhi1_tensor(params)
    iHWu.par =params.par;
for k = 1:ceil((n3+1)/2)
    iHWu.par.tmpDsch11(:,:,k) = sig./(1+sig*tau4-trans.par.Dsch11(:,:,k));
    iHWu.par.tmpHsch11(:,:,k) = sig./(1+sig*tau4-trans.par.Hsch11(:,:,k));
    iHWu.par.tmpDsch12(:,:,k) = sig./(1+sig*tau4-trans.par.Dsch12(:,:,k));
    iHWu.par.tmpHsch12(:,:,k) = sig./(1+sig*tau4-trans.par.Hsch12(:,:,k));
    iHWu.par.tmpDsch13(:,:,k) = sig./(1+sig*tau4-trans.par.Dsch13(:,:,k));
end

for k = 1:ceil((n3+1)/2)
    iHWu.par.D2sch11(:,:,k) = 0.5*(trans.par.Dsch11(:,:,k) + trans.par.Hsch11(:,:,k));
    iHWu.par.H2sch11(:,:,k) = 0.5*(trans.par.Dsch11(:,:,k) - trans.par.Hsch11(:,:,k));
    iHWu.par.D2sch12(:,:,k) = 0.5*(trans.par.Dsch12(:,:,k) + trans.par.Hsch12(:,:,k));
    iHWu.par.H2sch12(:,:,k) = 0.5*(trans.par.Dsch12(:,:,k) - trans.par.Hsch12(:,:,k));
    iHWu.par.D2sch13(:,:,k) = trans.par.Dsch13(:,:,k);
end

for k = 1:ceil((n3+1)/2)
    iHWu.par.tmp2Dsch11(:,:,k) = 0.5*(iHWu.par.tmpDsch11(:,:,k) + iHWu.par.tmpHsch11(:,:,k));
    iHWu.par.tmp2Hsch11(:,:,k) = 0.5*(iHWu.par.tmpDsch11(:,:,k) - iHWu.par.tmpHsch11(:,:,k));
    iHWu.par.tmp2Dsch12(:,:,k) = 0.5*(iHWu.par.tmpDsch12(:,:,k) + iHWu.par.tmpHsch12(:,:,k));
    iHWu.par.tmp2Hsch12(:,:,k) = 0.5*(iHWu.par.tmpDsch12(:,:,k) - iHWu.par.tmpHsch12(:,:,k));
    iHWu.par.tmp2Dsch13(:,:,k) = iHWu.par.tmpDsch13(:,:,k);
end

iHWu.par.tDsch11(:,:,1) = (iHWu.par.D2sch11(:,:,1)).*iHWu.par.tmp2Dsch11(:,:,1) +  iHWu.par.H2sch11(:,:,1).*iHWu.par.tmp2Hsch11(:,:,1);
iHWu.par.tHsch11(:,:,1) = (iHWu.par.D2sch11(:,:,1)).*iHWu.par.tmp2Hsch11(:,:,1) +  iHWu.par.H2sch11(:,:,1).*iHWu.par.tmp2Dsch11(:,:,1);
iHWu.par.tDsch12(:,:,1) = (iHWu.par.D2sch12(:,:,1)).*iHWu.par.tmp2Dsch12(:,:,1) +  iHWu.par.H2sch12(:,:,1).*iHWu.par.tmp2Hsch12(:,:,1);
iHWu.par.tHsch12(:,:,1) = (iHWu.par.D2sch12(:,:,1)).*iHWu.par.tmp2Hsch12(:,:,1) +  iHWu.par.H2sch12(:,:,1).*iHWu.par.tmp2Dsch12(:,:,1);

end