function [out, params] = getdirectionall(F,params,cgopts)
%% getdirectionall: obtain the update direction of all priaml and the dual variables.
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
tau1 = params.NEWT.tau1;
tau2 = params.NEWT.tau2;
tau3 = params.NEWT.tau3;
tau4 = params.NEWT.tau4;
CG_maxit = cgopts.CG_maxit;
pcgTol = cgopts.CG_tol;
Lchol = params.Lchol;
if params.Aboxflag == 1 && params.Axeqb == 0
    taux1 = params.NEWT.taux1;
    D1 = params.D.D1;
end

if params.boxflag == 1
    Fx3 = F.Fx3;
    taux3 = params.NEWT.taux3;
    D3 = params.D.D3;
else
    Fx3 = 0;
end

if params.fnonsmooth
    taux2 = params.NEWT.taux2;
    D2 = params.D.D2;
end
taux4 = params.NEWT.taux4;

sigma = params.sigma;
Fy = F.Fy;
Fz = F.Fz;
Fr = F.Fr;
Fv = F.Fv;
if params.Axeqb == 0
    Fx1 = F.Fx1;
end
Fx2 = F.Fx2;
Fx4 = F.Fx4;
params.cgres = 1;


D4 = params.D.D4;
iHW = D4;

tmptau = - (1 + sigma*taux4);
for p =1: params.nblock
    cone = params.pblk{p};
    if strcmp(cone.type, 's')
        iHW{p}.D12 = sigma * D4{p}.D12 ./ (1 + sigma * taux4 - D4{p}.D12); %iHW: \tilde{\Sigma}
        iHW{p}.D11 = 1 / taux4 ;
    elseif strcmp(cone.type, 'q')
        iHW{p}.D1 = - iHW{p}.D1 / sigma;
        iHW{p}.D2 = - iHW{p}.D2 / sigma;
        iHW{p}.shift = - iHW{p}.shift / sigma + 1 / sigma + taux4;
        % inverse
        iHW{p}.D1 = - iHW{p}.D1 ./ (iHW{p}.shift .* (iHW{p}.shift + iHW{p}.D1));
        iHW{p}.D2 = - iHW{p}.D2 ./ (iHW{p}.shift .* (iHW{p}.shift + iHW{p}.D2));
        iHW{p}.shift = 1 ./ iHW{p}.shift;
        % multiply by D (i.e. Ftz.par)
        temp_Dsch1 = D4{p}.D1 .* iHW{p}.D1 + D4{p}.shift .* iHW{p}.D1 + D4{p}.D1 .* iHW{p}.shift;
        temp_Dsch2 = D4{p}.D2 .* iHW{p}.D2 + D4{p}.shift .* iHW{p}.D2 + D4{p}.D2 .* iHW{p}.shift;
        iHW{p}.D1 = temp_Dsch1;
        iHW{p}.D2 = temp_Dsch2;
        iHW{p}.shift = D4{p}.shift .* iHW{p}.shift;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
        iHW{p} = sigma * D4{p} ./ (1 + sigma * taux4 - D4{p});
    elseif strcmp(cone.type, 'l1') || strcmp(cone.type, 'linfty') || strcmp(cone.type, 'max') || strcmp(cone.type, 'box') || strcmp(cone.type, 'topk') || strcmp(cone.type, 'l1topk')|| strcmp(cone.type, 'l1con') || strcmp(cone.type, 'linftycon') || strcmp(cone.type, 'l1linfty') || strcmp(cone.type, 'huber')
        rr1 = params.D.D4{p};
        
        iHW{p} = sigma * D4{p} ./ (1 + sigma * taux4 - D4{p});
    elseif strcmp(cone.type, 'l2')
        if isfield(iHW{p},'type') && iHW{p}.type == 2
            tmpcoe = 1/(taux4 + iHW{p}.coefficient/iHW{p}.nrmx/iHW{p}.sigma);
            iHW{p}.coe1 = tmpcoe*(1 - iHW{p}.coefficient/iHW{p}.nrmx);
            iHW{p}.coe2 = tmpcoe*(1/taux4/iHW{p}.sigma + 1)*iHW{p}.coefficient/iHW{p}.nrmx;
        end
    elseif strcmp(cone.type, 'l1l2')
        for j = 1:cone.size(2)
            if isfield(iHW{p},'type') && iHW{p}.type == 2
                tmpcoe = 1/(taux4 + iHW{p,j}.coefficient/iHW{p,j}.nrmx/iHW{p}.sigma);
                iHW{p,j}.coe1 = tmpcoe*(1 - iHW{p,j}.coefficient/iHW{p}.nrmx);
                iHW{p,j}.coe2 = tmpcoe*(1/taux4/iHW{p}.sigma + 1)*iHW{p,j}.coefficient/iHW{p}.nrmx;
            else
                iHW{p,j}.coe1 = 1/(taux4 + 1);
                iHW{p,j}.coe2 = 0;
            end
        end

    elseif strcmp(cone.type, 'nuclear')  || strcmp(cone.type, 'l2l2')
        iHW{p}.tmpDsch11(:,:) = sigma./(1+sigma*taux4-D4{p}.Dsch11(:,:));
        iHW{p}.tmpHsch11(:,:) = sigma./(1+sigma*taux4-D4{p}.Hsch11(:,:));
        iHW{p}.tmpDsch12(:,:) = sigma./(1+sigma*taux4-D4{p}.Dsch12(:,:));
        iHW{p}.tmpHsch12(:,:) = sigma./(1+sigma*taux4-D4{p}.Hsch12(:,:));
        iHW{p}.tmpDsch13(:,:) = sigma./(1+sigma*taux4-D4{p}.Dsch13(:,:));

        iHW{p}.tmp2Dsch11(:,:) = 0.5*(iHW{p}.tmpDsch11(:,:) + iHW{p}.tmpHsch11(:,:));
        iHW{p}.tmp2Hsch11(:,:) = 0.5*(iHW{p}.tmpDsch11(:,:) - iHW{p}.tmpHsch11(:,:));
        iHW{p}.tmp2Dsch12(:,:) = 0.5*(iHW{p}.tmpDsch12(:,:) + iHW{p}.tmpHsch12(:,:));
        iHW{p}.tmp2Hsch12(:,:) = 0.5*(iHW{p}.tmpDsch12(:,:) - iHW{p}.tmpHsch12(:,:));
        iHW{p}.tmp2Dsch13(:,:) = iHW{p}.tmpDsch13(:,:);

        iHW{p}.D2sch11(:,:) = 0.5*(D4{p}.Dsch11(:,:) + D4{p}.Hsch11(:,:));
        iHW{p}.H2sch11(:,:) = 0.5*(D4{p}.Dsch11(:,:) - D4{p}.Hsch11(:,:));
        iHW{p}.D2sch12(:,:) = 0.5*(D4{p}.Dsch12(:,:) + D4{p}.Hsch12(:,:));
        iHW{p}.H2sch12(:,:) = 0.5*(D4{p}.Dsch12(:,:) - D4{p}.Hsch12(:,:));
        iHW{p}.D2sch13(:,:) = D4{p}.Dsch13(:,:);

        iHW{p}.tDsch11(:,:) = (iHW{p}.D2sch11(:,:)).*iHW{p}.tmp2Dsch11(:,:) +  iHW{p}.H2sch11(:,:).*iHW{p}.tmp2Hsch11(:,:);
        iHW{p}.tHsch11(:,:) = (iHW{p}.D2sch11(:,:)).*iHW{p}.tmp2Hsch11(:,:) +  iHW{p}.H2sch11(:,:).*iHW{p}.tmp2Dsch11(:,:);
        iHW{p}.tDsch12(:,:) = (iHW{p}.D2sch12(:,:)).*iHW{p}.tmp2Dsch12(:,:) +  iHW{p}.H2sch12(:,:).*iHW{p}.tmp2Hsch12(:,:);
        iHW{p}.tHsch12(:,:) = (iHW{p}.D2sch12(:,:)).*iHW{p}.tmp2Hsch12(:,:) +  iHW{p}.H2sch12(:,:).*iHW{p}.tmp2Dsch12(:,:);
        iHW{p}.tDsch13(:,:) = iHW{p}.D2sch13(:,:).*iHW{p}.tmp2Dsch13(:,:);
    elseif strcmp(cone.type, 'fused')
        rr1 = params.D.D4{p}.rr1;
        AP = params.B(:,rr1);
        Ph = params.D.D4{p}.h(rr1);
        V1 = AP(:,Ph==1);
        numblk1 = length(params.D.D4{p}.nzcolidx);
        if (numblk1 > 0)
            PU1 = params.D.D4{p}.PU(:,params.D.D4{p}.nzcolidx);
        else
            PU1 = [];
        end
        tmpfactor = 1/taux4;
        tmpfactor2 = 1/taux4;
        iHW{p}.tmpfactor = tmpfactor;
        iHW{p}.tmpfactor2 = tmpfactor2;
    end
end

tmp4 = params.DPhiP(iHW, Fx4);

if isfield(params.D,'D1')
    tx1 = Fx1;
    for p = 1:params.nblock
        if iscell(params.D.D1)
        tx1{p}(params.D.D1{p}==0) = 0;
        else
        tx1(params.D.D1==0) = 0;
        end
    end
    tmp1 = tx1/taux1;
    rhs.rhsy = params.Amap(tmp4) - tmp1 - Fy;
elseif ~(isempty(params.lb) && isempty(params.ub))
    rhs.rhsy = -Fy + params.Amap(tmp4);
else
    rhs.rhsy = 0;
end
for p = 1:params.nblock
    if isfield(params,'f')
        if isfield(params.D,'D2')
            tx2 = Fx2;
            if  strcmp(params.f{p}.type,'l1')  || strcmp(params.f{p}.type, 'linfty') || strcmp(params.f{p}.type, 'max') || strcmp(params.f{p}.type, 'box') || strcmp(params.f{p}.type, 'topk') || strcmp(params.f{p}.type, 'l1topk')|| strcmp(params.f{p}.type, 'l1con') || strcmp(params.f{p}.type, 'linftycon') || strcmp(params.f{p}.type, 'huber')
                tx2{p}(params.D.D2{p}==0) = 0;
                tmp2 = tx2/taux2;
            elseif  strcmp(params.f{p}.type,'l2')
                D2 = params.D.D2;
                if isfield( D2,'type') &&  D2.type == 2
                    tmpcoe = 1/(taux2 +  D2.coefficient/ D2.nrmx / D2.sigma);
                    D2.coe1 = tmpcoe*(1 -  D2.coefficient/ D2.nrmx);
                    D2.coe2 = tmpcoe*(1/taux2*D2.sigma + 1)* D2.coefficient/ D2.nrmx;
                    tmp2 = D2.coe1*tx2 + D2.coe2*D2.X'*tx2*D2.X;
                else
                    tmp2 = tx2;
                end
            elseif strcmp(params.f{p}.type,'l2con')
                D2 = params.D.D2{p};
                if isfield( D2,'type') &&  D2.type == 2
                    tmpcoe = taux2 + (D2.nrmx - D2.sigma*D2.coefficient)/D2.sigma*D2.nrmx;
                    tmpcoe = tmpcoe^(-1);
                    tmpcoe2 = D2.sigma*D2.coefficient/D2.nrmx;
                    D2.coe1 = tmpcoe2*tmpcoe;
                    D2.coe2 = -tmpcoe2*tmpcoe;
                    tmp2 = D2.coe1*tx2 + D2.coe2*D2.X'*tx2{p}*D2.X;
                else
                    tmp2 = tx2;
                end
            elseif strcmp(cone.type, 'l1l2')
                D2 = params.D.D2;
                for j = 1:cone.size(2)
                    if isfield( D2,'type') &&  D2.type == 2
                        tmpcoe = 1/(taux2 +  D2.coefficient/ D2.nrmx);
                        D2.coe1 = tmpcoe*(1 -  D2.coefficient* D2.sigma/ D2.nrmx);
                        D2.coe2 = tmpcoe*(1/taux2 +  D2.sigma)* D2.coefficient/ D2.nrmx;
                        tmp2 = D2.coe1*tx2 + D2.coe2*D2.X'*tx2*D2.X;
                    else
                        tmp2 = tx2;
                    end
                end
            else
                tmp2{p} = params.f{p}.DPhi2(params.D.D2{p},tx2{p},taux2);
            end
            % rhs.rhsz = params.Bmap(tmp4)-tmp2-Fz;
            rhs.rhsz{p,1} = params.Bt{p}'*tmp4{p} -tmp2{p}-Fz{p};
        else
            if ~strcmp(params.pblk{p}.type,'fused') && ~strcmp(params.pblk{p}.type,'l1') || ~strcmp(params.f{p}.type,'l1')
                if isfield(params,'Bt')
                rhs.rhsz{p,1} = params.Bt{p}'*(tmp4{p})-Fz{p};
                else
                rhs.rhsz = params.Bmap(tmp4{p})-Fz;
                end
            else
                AP = params.B(:,rr1);
                rhs.rhsz = AP * tmp4{p}(rr1) - Fz; % AP * (Fxz.x(rr1).*Ph)*tempcoe2^2 = AP(:,Ph==1) * Fxz.x(Phindex)
            end
        end
    else
        rhs.rhsz =  zeros_like(F.Fz);
    end
end

if isfield(params.D,'D3')
    tx3 = Fx3;
    for i = 1:length(tx3)
        tx3{i}(params.D.D3{i}==0) = 0;
    end
    tmp3 = tx3/taux3;
    % tmp3 = sigma * params.D.D3{p} ./ (1 + sigma * taux3 - params.D.D3{p}).*tx3;
    rhs.rhsr = tmp4-tmp3-Fr;
else

    rhs.rhsr = zeros_like(F.Fr);
end

if params.Qflag == 1

    rhs.rhsv = -params.Qcmap(tmp4)-Fv;
    % end
else
    rhs.rhsv =  zeros_like(F.Fv);
end


iHWy = params.D;
iHWy.sigma = params.sigma;

for k = 1:params.nblock
    cone = params.pblk{k};
    if strcmp(cone.type, 's')
        iHWy.D12{k,1} = (iHWy.D4{k}.D12.*(1+sigma*taux4)*sigma)./(1+sigma*taux4-iHWy.D4{k}.D12);
        iHWy.D11{k,1} = (1+sigma*taux4)/taux4;
    elseif strcmp(cone.type,'q')
        iHWy.D11{k,1} = iHWy.D4{k}.D1 .* iHW{k}.D1 + iHWy.D4{k}.shift .* iHW{k}.D1 + iHWy.D4{k}.D1 .* (iHW{k}.shift + sigma);
        iHWy.D12{k,1} = iHWy.D4{k}.D2 .* iHW{k}.D2 + iHWy.D4{k}.shift .* iHW{k}.D2 + iHWy.D4{k}.D2 .* (iHW{k}.shift + sigma);
        iHWy.qshift = iHWy.D4{k}.shift .* (iHW{k}.shift + sigma);
    elseif strcmp(cone.type, 'l')
        iHWy.D11{k,1} = (iHWy.D4{k}.*(1+sigma*taux4)*sigma)./(1+sigma*taux4-iHWy.D4{k});
    elseif strcmp(cone.type, 'u')
        iHWy.D11{k,1} = (iHWy.D4{k}.*(1+sigma*taux4)*sigma)./(1+sigma*taux4-iHWy.D4{k});
    elseif strcmp(cone.type, 'b2l')
        % iHWy.b2lD = (iHWy.D4{k}.*(1+sigma*taux4)*sigma)./(1+sigma*taux4-iHWy.D4{k});
        iHWy.D11{k,1} = (iHWy.D4{k}.*(1+sigma*taux4)*sigma)./(1+sigma*taux4-iHWy.D4{k});
    elseif strcmp(cone.type, 'l1') || strcmp(cone.type, 'linfty') || strcmp(cone.type, 'max') || strcmp(cone.type, 'box') || strcmp(cone.type, 'topk') || strcmp(cone.type, 'l1topk')|| strcmp(cone.type, 'l1con') || strcmp(cone.type, 'linftycon') || strcmp(cone.type, 'l1linfty') || strcmp(cone.type, 'huber')
        iHWy.tDD = (1+sigma*taux4)/taux4;
        iHWy.D11{k,1} = iHWy.tDD*iHWy.D4{k};
    elseif strcmp(cone.type, 'l2')
        tmpcoe = 1/(taux4 + iHW{p}.coefficient/iHW{p}.nrmx/iHWy.D4{k}.sigma);
        iHWy.D11{k,1} = (tmpcoe + iHWy.D4{k}.sigma)*(1 - iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx);
        iHWy.D12{k,1} = tmpcoe*(2*iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx + iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx/iHWy.D4{k}.sigma/taux4 - iHW{p}.coefficient^2/iHWy.D4{k}.nrmx^2  ) + iHWy.D4{k}.sigma*iHW{p}.coefficient/iHWy.D4{k}.nrmx;
        % tmpcoe2 = iHWy.D4{k}.sigma*iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx;
        % iHWy.D11{k,1} =  tmpcoe*(1 - 2*tmpcoe2 + (tmpcoe2)^2) + iHWy.D4{k}.sigma*(1 - tmpcoe2);
        % iHWy.D12{k,1} =  tmpcoe*(2*tmpcoe2 - (tmpcoe2)^2 + iHWy.D4{k}.coefficient/taux4/iHWy.D4{k}.nrmx)  + iHWy.D4{k}.sigma*tmpcoe2;
    elseif strcmp(cone.type, 'l2con')
        D2 = params.D.D2{k};
        tmpcoe = taux2 + (iHWy.D4{k}.nrmx - iHWy.D4{k}.sigma*iHWy.D4{k}.coefficient)/iHWy.D4{k}.sigma*iHWy.D4{k}.nrmx;
        tmpcoe = tmpcoe^(-1);
        tmpcoe2 = iHWy.D4{k}.sigma*iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx;
        iHWy.D11{k,1}= tmpcoe2*tmpcoe;
        iHWy.D12{k,1}= -tmpcoe2*tmpcoe;
    elseif strcmp(cone.type, 'l1l2')
        for j = 1: cone.size(2)
            tmpcoe = 1/(taux4 + iHW{p}.coefficient/iHW{p}.nrmx/iHWy.D4{k}.sigma);
            % tmpcoe = 1/(taux4 + iHW{p}.coefficient/iHW{p}.nrmx);
            % tmpcoe2 = iHWy.D4{k}.sigma*iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx;
            % iHWy.tD4coe1(j) =  tmpcoe*(1 - 2*tmpcoe2 + (tmpcoe2)^2) + iHWy.D4{k}.sigma*(1 - tmpcoe2);
            % iHWy.tD4coe2(j) =  tmpcoe*(2*tmpcoe2 - (tmpcoe2)^2 + iHWy.D4{k}.coefficient/taux4/iHWy.D4{k}.nrmx)  + iHWy.D4{k}.sigma*tmpcoe2;
            % iHWy.D11{k,1}(j) =  tmpcoe*(1 - 2*tmpcoe2 + (tmpcoe2)^2) + iHWy.D4{k}.sigma*(1 - tmpcoe2);
            % iHWy.D12{k,1}(j) =  tmpcoe*(2*tmpcoe2 + (tmpcoe2)^2 + iHWy.D4{k}.coefficient/taux4/iHWy.D4{k}.nrmx)  + iHWy.D4{k}.sigma*tmpcoe2;
            iHWy.D11{k,1}(j) = (tmpcoe + iHWy.D4{k}.sigma)*(1 - iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx);
            iHWy.D12{k,1}(j) = tmpcoe*(2*iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx + iHWy.D4{k}.coefficient/iHWy.D4{k}.nrmx/iHWy.D4{k}.sigma/taux4 - iHW{p}.coefficient^2/iHWy.D4{k}.nrmx^2  ) + iHWy.D4{k}.sigma*iHW{p}.coefficient/iHWy.D4{k}.nrmx;
        end
    elseif strcmp(cone.type, 'nuclear') || strcmp(cone.type, 'l2l2')
        iHWy.sigma = sigma;

        iHWy.tmpDsch12{k,1} = zeros(size(iHWy.D4{k}.Dsch12));
        iHWy.tmpHsch12{k,1} = zeros(size(iHWy.D4{k}.Hsch12));
        iHWy.tmpDsch11{k,1} = zeros(size(iHWy.D4{k}.Dsch11));
        iHWy.tmpHsch11{k,1} = zeros(size(iHWy.D4{k}.Hsch11));

        iHWy.tmpDsch11{k,1}(:,:) = sigma./(1+sigma*taux4-iHWy.D4{k}.Dsch11(:,:));
        iHWy.tmpHsch11{k,1}(:,:) = sigma./(1+sigma*taux4-iHWy.D4{k}.Hsch11(:,:));
        iHWy.tmpDsch12{k,1}(:,:) = sigma./(1+sigma*taux4-iHWy.D4{k}.Dsch12(:,:));
        iHWy.tmpHsch12{k,1}(:,:) = sigma./(1+sigma*taux4-iHWy.D4{k}.Hsch12(:,:));
        iHWy.tmpDsch13{k,1}(:,:) = sigma./(1+sigma*taux4-iHWy.D4{k}.Dsch13(:,:));


        iHWy.tmp2Dsch11{k,1}(:,:) = 0.5*(iHWy.tmpDsch11{k,1}(:,:) + iHWy.tmpHsch11{k,1}(:,:));
        iHWy.tmp2Hsch11{k,1}(:,:) = 0.5*(iHWy.tmpDsch11{k,1}(:,:) - iHWy.tmpHsch11{k,1}(:,:));
        iHWy.tmp2Dsch12{k,1}(:,:) = 0.5*(iHWy.tmpDsch12{k,1}(:,:) + iHWy.tmpHsch12{k,1}(:,:));
        iHWy.tmp2Hsch12{k,1}(:,:) = 0.5*(iHWy.tmpDsch12{k,1}(:,:) - iHWy.tmpHsch12{k,1}(:,:));
        iHWy.tmp2Dsch13{k,1}(:,:) = iHWy.tmpDsch13(:,:);

        iHWy.D2sch11{k,1}(:,:) = 0.5*(iHWy.D4{k}.Dsch11(:,:) + iHWy.D4{k}.Hsch11(:,:));
        iHWy.H2sch11{k,1}(:,:) = 0.5*(iHWy.D4{k}.Dsch11(:,:) - iHWy.D4{k}.Hsch11(:,:));
        iHWy.D2sch12{k,1}(:,:) = 0.5*(iHWy.D4{k}.Dsch12(:,:) + iHWy.D4{k}.Hsch12(:,:));
        iHWy.H2sch12{k,1}(:,:) = 0.5*(iHWy.D4{k}.Dsch12(:,:) - iHWy.D4{k}.Hsch12(:,:));
        iHWy.D2sch13{k,1}(:,:) = iHWy.D4{k}.Dsch13(:,:);

        iHWy.tmp3Dsch11{k,1}(:,:) = iHWy.tmp2Dsch11{k,1}(:,:).*iHWy.D2sch11{k,1}(:,:) + iHWy.tmp2Hsch11{k,1}(:,:).*iHWy.H2sch11{k,1}(:,:);
        iHWy.tmp3Hsch11{k,1}(:,:) = iHWy.tmp2Dsch11{k,1}(:,:).*iHWy.H2sch11{k,1}(:,:) + iHWy.tmp2Hsch11{k,1}(:,:).*iHWy.D2sch11{k,1}(:,:);
        iHWy.tmp3Dsch12{k,1}(:,:) = iHWy.tmp2Dsch12{k,1}(:,:).*iHWy.D2sch12{k,1}(:,:) + iHWy.tmp2Hsch12{k,1}(:,:).*iHWy.H2sch12{k,1}(:,:);
        iHWy.tmp3Hsch12{k,1}(:,:) = iHWy.tmp2Dsch12{k,1}(:,:).*iHWy.H2sch12{k,1}(:,:) + iHWy.tmp2Hsch12{k,1}(:,:).*iHWy.D2sch12{k,1}(:,:);
        iHWy.tmp3Dsch13{k,1}(:,:) = iHWy.tmp2Dsch13{k,1}(:,:).*iHWy.D2sch13{k,1}(:,:);

        iHWy.D3sch11{k,1}(:,:) = (iHWy.D2sch11{k,1}(:,:)).*iHWy.tmp3Dsch11{k,1}(:,:) + iHWy.H2sch11{k,1}(:,:).*iHWy.tmp3Hsch11{k,1}(:,:);
        iHWy.H3sch11{k,1}(:,:) = (iHWy.D2sch11{k,1}(:,:)).*iHWy.tmp3Hsch11{k,1}(:,:) + iHWy.H2sch11{k,1}(:,:).*iHWy.tmp3Dsch11{k,1}(:,:);
        iHWy.D3sch12{k,1}(:,:) = (iHWy.D2sch12{k,1}(:,:)).*iHWy.tmp3Dsch12{k,1}(:,:) + iHWy.H2sch12{k,1}(:,:).*iHWy.tmp3Hsch12{k,1}(:,:);
        iHWy.H3sch12{k,1}(:,:) = (iHWy.D2sch12{k,1}(:,:)).*iHWy.tmp3Hsch12{k,1}(:,:) + iHWy.H2sch12{k,1}(:,:).*iHWy.tmp3Dsch12{k,1}(:,:);
        iHWy.D3sch13{k,1}(:,:) = iHWy.D2sch13{k,1}(:,:).*iHWy.tmp3Dsch13{k,1}(:,:);

        iHWy.tDsch11{k,1}(:,:) = iHWy.D3sch11{k,1}(:,:)+ sigma*iHWy.D2sch11{k,1}(:,:);
        iHWy.tHsch11{k,1}(:,:) = iHWy.H3sch11{k,1}(:,:)+ sigma*iHWy.H2sch11{k,1}(:,:);
        iHWy.tDsch12{k,1}(:,:) = iHWy.D3sch12{k,1}(:,:)+ sigma*iHWy.D2sch12{k,1}(:,:);
        iHWy.tHsch12{k,1}(:,:) = iHWy.H3sch12{k,1}(:,:)+ sigma*iHWy.H2sch12{k,1}(:,:);
        iHWy.tDsch13{k,1}(:,:) = iHWy.D3sch13{k,1}(:,:)+ sigma*iHWy.D2sch13{k,1}(:,:);
    elseif strcmp(cone.type, 'fused')
        tempcoe1 = tmpfactor + sigma;
        tempcoe1 = sqrt(tempcoe1);
        tempcoe2 =  1/taux4 + sigma;
        tempcoe2 = sqrt(tempcoe2);
        iHWy.lenP = length(Ph);
        iHWy.PU2 =  PU1*tempcoe1; %%weighted
        iHWy.V1 = V1*tempcoe2; %%weighted
        iHWy.PU = params.D.D4{k}.PU;
        iHWy.AP =AP;
        iHWy.Ph = Ph;
        iHWy.PU1 = PU1;
        iHWy.numblk1 = numblk1;
        iHWy.info_u = params.D.D4{k}.info;
    end
end

if isfield(iHWy,'D1')
    iHWy.tD1 = (1+sigma*taux1)/taux1*iHWy.D1;
end

if isfield(iHWy,'D3')
    iHWy.tD3 = (1+sigma*taux3)/taux3*iHWy.D3;
end
for p =1:params.nblock
    if isfield(params,'f') && params.fnonsmooth == 1
        if strcmp(params.f{p}.type,'l1')  || strcmp(params.f{p}.type, 'linfty') || strcmp(params.f{p}.type, 'max') || strcmp(params.f{p}.type, 'box') || strcmp(params.f{p}.type, 'topk') || strcmp(params.f{p}.type, 'l1topk')|| strcmp(params.f{p}.type, 'l1con') || strcmp(params.f{p}.type, 'linftycon') || strcmp(params.f{p}.type, 'l1linfty') || strcmp(params.f{p}.type, 'huber')
            iHWy.D21{p} = (1+sigma*taux2)/taux2*iHWy.D2{p};
        elseif  strcmp(params.f{p}.type,'l2')
            D2 = params.D.D2{p};
            tmpcoe = 1/(taux2 + D2.coefficient/D2.nrmx);
            iHWy.D11{p} = (tmpcoe + iHWy.D2{k}.sigma)*(1 - iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx);
            iHWy.D12{p} = tmpcoe*(2*iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx + iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx/iHWy.D2{k}.sigma/taux4 - iHWy.D2{k}.coefficient^2/iHWy.D2{k}.nrmx^2  ) + iHWy.D2{k}.sigma*iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx;
            % tmpcoe2 = iHWy.D2{k}.sigma*iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx;
            % iHWy.D21{p} =  tmpcoe*(1 - 2*tmpcoe2 + (tmpcoe2)^2) + iHWy.D2{k}.sigma*(1 - tmpcoe2);
            % iHWy.D22{p} =  tmpcoe*(2*tmpcoe2 - (tmpcoe2)^2 + iHWy.D2{k}.coefficient/taux2/iHWy.D2{k}.nrmx)  + iHWy.D2{k}.sigma*tmpcoe2;
        elseif strcmp(params.f{p}.type,'l2con')
            D2 = params.D.D2{p};
            tmpcoe = taux2 + (iHWy.D2{k}.nrmx - iHWy.D2{k}.sigma*iHWy.D2{k}.coefficient)/iHWy.D2{k}.sigma*iHWy.D2{k}.nrmx;
            tmpcoe = tmpcoe^(-1);
            tmpcoe2 = iHWy.D2{k}.sigma*iHWy.D2{k}.coefficient/iHWy.D2{k}.nrmx;
            iHWy.D21{p}= tmpcoe2^2*tmpcoe + iHWy.D2{k}.sigma*tmpcoe2;
            iHWy.D22{p} = -tmpcoe2^2*tmpcoe - iHWy.D2{k}.sigma*tmpcoe2;
        end
    elseif isfield(params,'f') && params.fnonsmooth == 0
        if strcmp(params.f{p}.type,'square')
            iHWy.tD2 = params.f{p}.coefficient *2;
        elseif strcmp(params.f{p}.type,'exp')
            iHWy.tD2 = iHWy.g2;
        elseif strcmp(params.f{p}.type,'logdet')
            iHWy.tD2 = iHWy.g2;
        end
    end
end


iHWy.epsilon1 = tau1;
iHWy.epsilon2 = tau2;
iHWy.epsilon3 = tau3;
iHWy.epsilon4 = tau4;
L = struct();

[outcg,relres,flag] = linearsolveall_ssncvx(@matvec_ssncvx,rhs,iHWy,L, pcgTol, CG_maxit, params);

%  Aty = params.ATmap(outcg.dy)+ params.BTmap(outcg.dz{1})  - params.Qmap(outcg.dv) ;
% tmp = iHWy.tD4{p}.*Aty;
% % tmp = iHWy.tD4coe1 *Aty + iHWy.tD4coe2 *iHWy.D4{p}.X'*Aty*iHWy.D4{p}.X;
% Ax1 = sigma*(tmp'*params.At{p})';
% Ax2 = sigma*(tmp'*params.Bt{p})';
% Ax3 = sigma*tmp;
% Ax4 = -sigma*(tmp'*params.Qt{p})';
%
% Ax1 = Ax1 + iHWy.epsilon1 * outcg.dy;
% Ax2 = Ax2 + outcg.dz .* (iHWy.tD2 + iHWy.epsilon2);
% Ax3 = Ax3 + outcg.dr .* (iHWy.tD3 + iHWy.epsilon3);
% Ax4 = Ax4 + outcg.dv *  iHWy.epsilon4 + params.Qmap(outcg.dv);
%
% norm(Ax1 - rhs.rhsy)
% norm(Ax2 - rhs.rhsz)
% norm(Ax3 - rhs.rhsr)
% norm(Ax4 - rhs.rhsv)

params.cgres = [params.cgres relres(end)];
params.flag = flag;
params.cgiter = length(relres);



for p = 1:length(params.pblk)
    cone = params.pblk{p};
    if ~strcmp(cone.type, 'fused') && ~strcmp(cone.type, 's') && ~strcmp(cone.type, 'l1') || ~(length(params.pblk) == 1) || ~params.fflag  ||  ~strcmp(params.f{p}.type,'square')
        dx4 = params.ATmap(outcg.dy) + params.BTmap(outcg.dz) + params.idmap(outcg.dr) - params.Qcmap(outcg.dv);
    elseif strcmp(cone.type, 's')
        dx4 = params.ATmap(outcg.dy) + params.BTmap(outcg.dz) + params.idmap(outcg.dr) - params.Qcmap(outcg.dv);
    elseif strcmp(cone.type, 'fused') || strcmp(cone.type, 'l1')
        if isfield(params,'Bt')
            dx4{p} = zeros(size(params.Bt{p},1),1);
        else
            dx4 = zeros(size(params.B,2),1);
        end
        if strcmp(cone.type, 'l1') && isfield(params,'Bt')
        AP = params.Bt{p}(rr1,:)';
            else
        AP = params.B(:,rr1);
        end
 
        
        dx4{p}(rr1) = AP'*outcg.dz{p};
    end
    if strcmp(cone.type, 'nuclear') || strcmp(cone.type, 'l2l2')
        params.D.D4{p}.tDsch11(:,:) = 0.5*(params.D.D4{p}.Dsch11(:,:) + params.D.D4{p}.Hsch11(:,:));
        params.D.D4{p}.tHsch11(:,:) = 0.5*(params.D.D4{p}.Dsch11(:,:) - params.D.D4{p}.Hsch11(:,:));
        params.D.D4{p}.tDsch12(:,:) = 0.5*(params.D.D4{p}.Dsch12(:,:) + params.D.D4{p}.Hsch12(:,:));
        params.D.D4{p}.tHsch12(:,:) = 0.5*(params.D.D4{p}.Dsch12(:,:) - params.D.D4{p}.Hsch12(:,:));
        params.D.D4{p}.tDsch13(:,:) = params.D.D4{p}.Dsch13(:,:);
    end
end
dx4 = params.DPhiP(params.D.D4,dx4);



for p = 1:length(params.pblk)
    cone = params.pblk{p};
    if strcmp(cone.type, 's')
        iHWx =params.D.D4;
            iHWx{p}.D12 = sigma * iHWx{p}.D12 ./ (1 + sigma * taux4 - iHWx{p}.D12); %iHW: \tilde{\Sigma}
            iHWx{p}.D11 = 1 / taux4 ;
        rhs2{p} =  -(dx4{p} - Fx4{p})/tmptau;

        tmp = DPhi_all({iHWx{p}},rhs2{p},{params.pblk{p}});
        dx4{p} = tmp{1} + sigma * rhs2{p};
    elseif strcmp(cone.type, 'q')
        iHWx =params.D.D4;
        iHWx{p}.D1 = - iHWx{p}.D1 / sigma;
        iHWx{p}.D2 = - iHWx{p}.D2 / sigma;
        iHWx{p}.shift = - iHWx{p}.shift / sigma + 1 / sigma + taux4;
        % inverse
        iHWx{p}.D1 = - iHWx{p}.D1 ./ (iHWx{p}.shift .* (iHWx{p}.shift + iHWx{p}.D1));
        iHWx{p}.D2 = - iHWx{p}.D2 ./ (iHWx{p}.shift .* (iHWx{p}.shift + iHWx{p}.D2));
        iHWx{p}.shift = 1 ./ iHWx{p}.shift;
        % multiply by tmptau and subtract by sigma I
        iHWx{p}.D1 = - iHWx{p}.D1 * tmptau;
        iHWx{p}.D2 = - iHWx{p}.D2 * tmptau;
        iHWx{p}.shift = - iHWx{p}.shift * tmptau - sigma;
        % rhs2 =  -(dx4 - Fx4)/tmptau;
        rhs2{p} =  -(dx4{p} - Fx4{p})/tmptau;

        % dx4_tmp= params.DPhiP(iHWx,rhs2);
        tmp = DPhi_all({iHWx{p}},rhs2{p},{params.pblk{p}});
        dx4{p} = tmp{1} + sigma * rhs2{p};
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
        iHWx =params.D.D4;
        iHWx{p} = sigma * iHWx{p} ./ (1 + sigma * taux4 - iHWx{p}); %iHW: \tilde{\Sigma}
        rhs2 =  -(dx4 - Fx4)/tmptau;
        dx4_tmp= params.DPhiP(iHWx,rhs2);
        dx4{p} = dx4_tmp{p} + sigma * rhs2{p};
    elseif strcmp(cone.type, 'l1') || strcmp(cone.type, 'linfty') || strcmp(cone.type, 'max') || strcmp(cone.type, 'box') || strcmp(cone.type, 'topk') || strcmp(cone.type, 'l1topk') || strcmp(cone.type, 'l1con') || strcmp(cone.type, 'linftycon') || strcmp(cone.type, 'l1linfty') || strcmp(cone.type, 'huber')
        dx4{p} =  -(dx4{p} - Fx4{p})/tmptau;
        iHWxtmp = params.D.D4;
        iHWx{p} = 1/taux4*iHWxtmp{p};
        tmp2 =  iHWx{p}.*dx4{p};
        dx4{p} = tmp2 + sigma * dx4{p};
    elseif  strcmp(cone.type, 'l2')
        dx4{p} = (dx4{p} - Fx4{p});
            iHWx{p} = params.D.D4{p};
            tmpcoe = 1/(taux4 + iHW{p}.coefficient/iHW{p}.nrmx/iHWx{p}.sigma);
            iHWx{p}.coe1 = tmpcoe;
            iHWx{p}.coe2 = iHWx{p}.coefficient/iHWx{p}.sigma/taux4/iHWx{p}.nrmx*tmpcoe;
        % dx4= params.DPhiP(iHWx,dx4);
        tmp2= DPhi_all(iHWx,dx4{p},{params.pblk{p}});
        dx4{p}=  tmp2{1};
    elseif  strcmp(cone.type, 'l1l2')
        dx4{p} = (dx4{p} - Fx4{p});
            for j = 1:cone.size(2)
                iHWx{p,j} = params.D.D4{p,j};
                tmpcoe = 1/(taux4 + iHW{p,j}.coefficient/iHW{p,j}.nrmx);
                iHWx{p,j}.coe1 = tmpcoe;
                iHWx{p,j}.coe2 = iHWx{p,j}.coefficient/taux4/iHWx{p,j}.nrmx*tmpcoe;
            end
        tmp2= DPhi_all(iHWx,dx4{p},{params.pblk{p}});
        dx4{p}=  tmp2{1};
    elseif strcmp(cone.type, 'nuclear') || strcmp(cone.type, 'l2l2')
        dx4{p} = (dx4{p} - Fx4{p});
        iHWx =params.D.D4;
        iHWx{p}.Dsch11(:,:) = (sigma)./(1+sigma*taux4-D4{p}.Dsch11(:,:));
        iHWx{p}.Hsch11(:,:) = (sigma)./(1+sigma*taux4-D4{p}.Hsch11(:,:));
        iHWx{p}.Dsch12(:,:) = (sigma)./(1+sigma*taux4-D4{p}.Dsch12(:,:));
        iHWx{p}.Hsch12(:,:) = (sigma)./(1+sigma*taux4-D4{p}.Hsch12(:,:));
        iHWx{p}.Dsch13(:,:) = (sigma)./(1+sigma*taux4-D4{p}.Dsch13(:,:));

        iHWx{p}.tDsch11(:,:) = 0.5*(iHWx{p}.Dsch11(:,:) + iHWx{p}.Hsch11(:,:))-sigma./(1+sigma*taux4);
        iHWx{p}.tHsch11(:,:) = 0.5*(iHWx{p}.Dsch11(:,:) - iHWx{p}.Hsch11(:,:));
        iHWx{p}.tDsch12(:,:) = 0.5*(iHWx{p}.Dsch12(:,:) + iHWx{p}.Hsch12(:,:))-sigma./(1+sigma*taux4);
        iHWx{p}.tHsch12(:,:) = 0.5*(iHWx{p}.Dsch12(:,:) - iHWx{p}.Hsch12(:,:));
        iHWx{p}.tDsch13(:,:) = iHWx{p}.Dsch13(:,:)  -sigma./(1+sigma*taux4);

        tmp2= DPhi_all({iHWx{p}},dx4{p},{params.pblk{p}});
        dx4{p} = tmp2{1} + 1/(tau4 + 1/sigma)*dx4{p};
    elseif strcmp(cone.type, 'fused')
        iHWx =params.D.D4;
        tmpdx4 = dx4;
        dx =   tmpdx4{p} - Fx4{p} ;
        tmp2 = dx;%  hat{D}-1 has indetityiHWx{p}.
        if (numblk1 > 0)
            dx4{p}(rr1) = iHWx{p}.PU*(iHWx{p}.PU'*dx(rr1))./(1+taux4 * sigma)./(taux4)+ iHWx{p}.Ph.*dx(rr1)./(taux4*(1+taux4*sigma)) ;
        else
            dx4{p}(rr1) = iHWx{p}.Ph.*dx(rr1)./(taux4*(1+taux4*sigma)) ;
        end
        dx4 = dx4 + sigma/(1+taux4*sigma)*tmp2;
    end
end


if isfield(params.D,'D1')
    tzb = params.D.D1 .* outcg.dy;
    tmp4 = - (tzb + Fx1);
    dx1 = (tmp4*sigma/(1+taux1*sigma) + (params.D.D1 .* tmp1)*1/(taux1*(1+taux1*sigma)));
else
    dx1 = 0;
end

if isfield(params.D,'D2')
    for p = 1:params.nblock
        if strcmp(params.f{p}.type,'l1')  || strcmp(params.f{p}.type, 'linfty') || strcmp(params.f{p}.type, 'max') || strcmp(cone.type, 'box') || strcmp(params.f{p}.type, 'topk') || strcmp(params.f{p}.type, 'l1topk')|| strcmp(params.f{p}.type, 'l1con') || strcmp(params.f{p}.type, 'linftycon')
            tzb = params.D.D2.* outcg.dz;
            tmp2 = - (tzb + Fx2);
            dx2 = (tmp2*sigma/(1+taux2*sigma) + (params.D.D2 .* tmp2)*1/(taux2*(1+taux2*sigma)));
        elseif strcmp(params.f{p}.type,'l2')
            iHWx{p} = params.D.D2{p};
            if iHWx{p}.type == 1
                dx2 = outcg.dz;
            else
                dx2 = iHWx{p}.coe1*outcg.dz{p} + iHWx{p}.coe2*iHWx{p}.X'*outcg.dz{p}*iHWx{p}.X;
            end
            dx2 = (dx2 - Fx2);
            tmpcoe = 1/(taux2 + iHWx{p}.coefficient/iHWx{p}.nrmx);
            iHWx{p}.coe1 = tmpcoe;
            iHWx{p}.coe2 = iHWx{p}.coefficient/taux2/iHWx{p}.nrmx*tmpcoe;
            if iHWx{p}.type == 1
                dx2 = outcg.dz;
            else
                dx2{p} = iHWx{p}.coe1*dx2{p} + iHWx{p}.coe2*iHWx{p}.X'*dx2{p}*iHWx{p}.X;
            end
        elseif strcmp(params.f{p}.type,'l2con')
            iHWx{p} = params.D.D2{p};
            if iHWx{p}.type == 1
                dx2 = outcg.dz;
            else
                dx2 = iHWx{p}.coe1*outcg.dz{p} + iHWx{p}.coe2*iHWx{p}.X'*outcg.dz{p}*iHWx{p}.X;
            end
            dx2 = (dx2 - Fx2);
            tmpcoe = taux2 + (iHWy.D2{p}.nrmx - iHWy.D2{p}.sigma*iHWy.D2{p}.coefficient)/iHWy.D2{p}.sigma*iHWy.D2{p}.nrmx;
            tmpcoe = tmpcoe^(-1);
            tmpcoe2 = iHWy.D2{p}.coefficient/iHWy.D2{p}.nrmx;
            iHWx{p}.coe1 = tmpcoe;
            iHWx{p}.coe2 = tmpcoe*tmpcoe2/(1 - tmpcoe*tmpcoe2);
            if iHWx{p}.type == 1
                dx2 = outcg.dz;
            else
                dx2{p} = iHWx{p}.coe1*dx2{p} + iHWx{p}.coe2*iHWx{p}.X'*dx2{p}*iHWx{p}.X;
            end
        else
            dx2{p} = params.f{p}.DPhi3(params.D.D2{p},iHWy.D2{p},outcg.dz{p},Fx2{p},taux2);
        end
    end
else
    dx2 = 0;
end

if isfield(params.D,'D3')
    tzb = params.D.D3 .* outcg.dr;
    tmp3 = - (tzb + Fx3);
    dx3 = (tmp3*sigma/(1+taux3*sigma) + (params.D.D3 .* tmp3)*1/(taux3*(1+taux3*sigma)));
else
    dx3 = 0;
end
%%
out.dy = outcg.dy;
out.dz = outcg.dz;
out.dr = outcg.dr;
out.dv = outcg.dv;

out.dx1 = dx1;
out.dx2 = dx2;
out.dx3 = dx3;
out.dx4 = dx4;
%%
%
% D_lmut = @(x) params.DPhiP( params.D.D4, x);
% Atdy = params.ATmap(outcg.dy);
% Btdz = params.BTmap(outcg.dz{1});
% QTdv{1} = params.Qmap(outcg.dv);
% % dzbtmp{1} = dr;
% summ = Atdy + Btdz +  - QTdv ; %outcg.dr
% tmp_test1 = params.DPhiP( params.D.D4, summ);
% tmpx = params.DPhiP( params.D.D4, dx4);
% N1 = sigma * params.Amap(tmp_test1);
% N2 = sigma * params.Bmap(tmp_test1{1});
% N3 = sigma * tmp_test1;
% N4 = -sigma * params.Qmap(tmp_test1);
% N18 =  params.Amap(tmpx);
% N28 =  params.Bmap(tmpx{1});
% N26 = -params.D.D2.*dx2;
% N38 =  tmpx;
% N48 =  -params.Qmap(tmpx);
% % params.D.D3 = 0;
% % N37 = - params.D.D3.* dx3;
%
% N66 = (1 / sigma + taux2) * dx2 - 1 / sigma * params.D.D2.*dx2;
% if ~exist('taux3','var')
%     taux3 = 0;
% end
% % N77 = (1 / sigma + taux3) * dx3 - 1 / sigma * params.D.D3.*dx3;
% N88 = (1 / sigma + taux4) * dx4 - 1 / sigma * params.DPhiP(params.D.D4,dx4);
%
% N62 =  params.D.D2.* outcg.dz{1};
% % N73 =  params.D.D3.* outcg.dr;
% N81{1} = -params.ATmap(outcg.dy);
% N81 = params.DPhiP( params.D.D4, N81);
% N82{1} = -params.BTmap(outcg.dz{1});
% N82 = params.DPhiP( params.D.D4, N82);
% % N83 = - params.DPhiP( params.D.D4, outcg.dr);
% N83 = 0;
% N84 = params.DPhiP( params.D.D4, QTdv);
% % N84 = - N48;
%
% % N11dy = sigma * Amap(D_lmut(Atdy)) + tau1 * dy;
% % N12dz = sigma* Amap(tmpz);
% % N13dx = Amap(tmpx);
% %
% % N21dy = sigma*D_lmut(Atdy);
% % N22dz = sigma*tmpz + sigma*Ftz.Dh.*dzb + tau1 * dzb;
% % N23dx = tmpx;
% % N24dq = -Ftz.Dh.* dzb;
% %
% % N51dy = N21dy/-sigma;
% % N53dx = (1 / sigma + tau3) * dx - 1 / sigma * D_lmut(dx);
% %
% % N31dy =  N21dx/-sigma ;
% % N44dq = (1 / sigma + tau4) * dq - 1 / sigma * Ftz.Dh.*dq;
% resy = N1 + N18 + tau1 * outcg.dy + Fy;
% norm(resy)
% resz = N2 + N26 + N28 + (tau2 + params.D.D2 ) * outcg.dz + Fz;
% norm(resz)
% % resr = N3 + sigma*params.D.D3.* outcg.dr + N38 +  N37 + (tau3) * outcg.dr + Fr;
% % norm(resr);
% % resv = N4 + N48 +  (tau4) * outcg.dv + params.Qmap(outcg.dv) + Fv;
% % norm(resv)
%
% resx2 = N62 + N66 + Fx2;
% norm(resx2)
% % resx3 = N73 + N77 + Fx3;
% % norm(resx3)
% resx4 = N81 + N82 + N83 + N84 + N88  + Fx4;
% norm(resx4)
% 1;
% % params.D.D3 = []
% newton_res = sqrt(norm(resy) ^ 2 + norm(resz) ^ 2 )  ;
% params.newton_res = newton_res;
end