function out = DPhi_all(iHW, d, pblk)
%% DPhi_all: compuate the D(d) operation, D is the generalized Jacboian operator 
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
plength = length(pblk);
for i = 1:plength
    cone = pblk{i};
    if strcmp(cone.type, 'l1') || strcmp(cone.type, 'linfty') || strcmp(cone.type, 'max') || strcmp(cone.type, 'topk') || strcmp(cone.type, 'l1topk') ||  strcmp(cone.type, 'linftycon') ||  strcmp(cone.type, 'box') ||  strcmp(cone.type, 'l1con') || strcmp(cone.type, 'huber')
        if iscell(d)
            out{i,1} = iHW{i}.*d{i};
        else
            out{i,1} = iHW{i}.*d;
        end
    elseif strcmp(cone.type, 'l1linfty')
        if iscell(d)
            out{i,1} = iHW{i}.*d{i};
        else
            out{i,1} = iHW{i}.*d;
        end
    elseif strcmp(cone.type, 's')
        n = sum(cone.size);
        if iscell(d)
            outtmp =  d{i};
        else
            outtmp = d;
        end
        Y{i} = sparse(n,n);
        rr = size(iHW{i}.P1,2);
        if (rr > 0 && rr < n)
            if (rr <= n/2)
                tmp0 = iHW{i}.P1'*outtmp;                  % tmp0: U = Q_{\alpha}^\top H
                tmp1 = (tmp0*iHW{i}.P1)*iHW{i}.P1';      % iHW{i}.P1: Q_{\alpha}
                tmp2 = iHW{i}.D12.*(tmp0*iHW{i}.P2);  % iHW{i}.D12: \nu_{\alpha \bar{\alpha}}; iHW{i}.P2: Q_{\bar{\alpha}}
                tmp2 = tmp2*iHW{i}.P2';                  % temp2: \nu_{\alpha \bar{\alpha}} \circ (U Q_{\bar{\alpha}}) Q_{\bar{\alpha}}^\top
                tmp3 = iHW{i}.P1*(0.5*iHW{i}.D11*tmp1 + tmp2);
                Y{i} = tmp3+tmp3';
            else
                tmp0 = iHW{i}.P2'*outtmp;
                tmp1 = (tmp0*iHW{i}.P2)*iHW{i}.P2';
                tmp2 = (iHW{i}.D11-iHW{i}.D12').*(tmp0*iHW{i}.P1);
                tmp2 = tmp2*iHW{i}.P1';
                tmp3 = iHW{i}.P2*(0.5*iHW{i}.D11*tmp1 + tmp2);
                Y{i} = iHW{i}.D11*outtmp-tmp3-tmp3';
            end
        elseif (rr == n)
            Y{i} = iHW{i}.D11* outtmp;
        end
        out{i,1} = Y{i};
    elseif strcmp(cone.type, 'q')
        n = sum(cone.size);
        if iscell(d)
            outtmp = d{i};
        else
            outtmp = d;
        end
        tmp1 = iHW{i}.Q1' * outtmp .* iHW{i}.D1;
        tmp1 = repelem(tmp1, cone.size, 1) .* iHW{i}.P1;
        tmp2 = iHW{i}.Q2' * outtmp .* iHW{i}.D2;
        tmp2 = repelem(tmp2, cone.size, 1) .* iHW{i}.P2;
        tmp3 = repelem(iHW{i}.shift, cone.size, 1) .* outtmp;
        out{i,1} = tmp1 + tmp2 + tmp3;
    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
        if iscell(d)
            outtmp = d{i};
        else
            outtmp = d;
        end
        out{i,1} = iHW{i} .* outtmp;
    elseif strcmp(cone.type, 'fused')
        if iscell(d)
            outtmp =  d{i};
        else
            outtmp = d;
        end
        if norm(double(iHW{i}.rr1)) ~= 0
        rr1 = iHW{i}.rr1;
        tmp = outtmp(rr1);
        outtmp(rr1) =   iHW{i}.Ph.*tmp*iHW{i}.tmpfactor + (iHW{i}.PU*iHW{i}.PU')*tmp*iHW{i}.tmpfactor2;
        out{i,1} = outtmp;
        else
        out{i,1} = zeros(size(iHW{i}.rr1));
        end
    elseif strcmp(cone.type, 'l2') || strcmp(cone.type, 'l2con')
        if iHW{i}.type == 1
            if iscell(d)
                out{i,1} = d{i};
            else
                out{i,1} = d;
            end
        else
            if iscell(d)
                out{i,1} = iHW{i}.coe1*d{i} + iHW{i}.coe2*iHW{i}.X'*d{i}*iHW{i}.X;
            else
                out{i,1} = iHW{i}.coe1*d + iHW{i}.coe2*iHW{i}.X'*d*iHW{i}.X;
            end
        end
    elseif strcmp(cone.type, 'l1l2')
        if iscell(d)
            [m,n] = size(d{i});
            tmpd = d{i};
        else
            [m,n] = size(d);
            tmpd = d;
        end
        outtmp = zeros(m,n);
        for j = 1:n
            if iHW{i,j}.type == 1
                outtmp(:,j) = tmpd(:,j);
            else
                outtmp(:,j) = iHW{i,j}.coe1*tmpd(:,j) + iHW{i,j}.coe2*iHW{i,j}.X'*tmpd(:,j)*iHW{i,j}.X;
            end
        end
        out{i,1} = outtmp;
    elseif strcmp(cone.type, 'nuclear') || strcmp(cone.type, 'l2l2')
        if ~iscell(d)
            tmp{i} = d;
        else
            tmp = d;
        end
        [n1,n2,~] = size(tmp{i});
        rr = iHW{i}.trank;
        out{i,1} = zeros(size(tmp{i}));
        %     Aty = Atyfun_sdpnal(pblk,At{i},y1) + y2;
        if (rr > 0 )
            % iy = fft(y,[],3);
            % tmp1 = compute_svd(params.U1T,iy,params.V1,1);
            tmp1 = iHW{i}.U1T*tmp{i}*iHW{i}.V1;
            if rr < n1/2
                tmp2 = iHW{i}.U1T*tmp{i}*iHW{i}.V2;
                tmp3 = iHW{i}.U2T*(tmp{i}*iHW{i}.V1);
                % tmp2 = compute_svd(iHW{i}ams.U1T,iy,iHW{i}ams.V2,1);
                % tmp3 = compute_svd(iHW{i}ams.U2T,iy,iHW{i}ams.V1,2);
            else
                tmp2 = iHW{i}.U1T*(tmp{i}*iHW{i}.V2);
                tmp3 = iHW{i}.U2T*tmp{i}*iHW{i}.V1;
                % tmp2 = compute_svd(iHW{i}ams.U1T,iy,iHW{i}ams.V2,2);
                % tmp3 = compute_svd(iHW{i}ams.U2T,iy,iHW{i}ams.V1,1);
            end
            tmp11 = tmp1.* iHW{i}.tDsch11 + tmp1'.* iHW{i}.tHsch11;
            tmp12 = tmp2.* iHW{i}.tDsch12 + tmp3'.* iHW{i}.tHsch12;
            tmp21 = tmp3.* iHW{i}.tDsch12'+ tmp2'.* iHW{i}.tHsch12';
            H11 = iHW{i}.U1*tmp11*iHW{i}.V1T;
            if rr <n1/2
                H12 = iHW{i}.U1*(tmp12*iHW{i}.V2T);
                H21 = iHW{i}.U2*tmp21*iHW{i}.V1T;
            else
                H12 = iHW{i}.U1*tmp12*iHW{i}.V2T;
                H21 = iHW{i}.U2*(tmp21*iHW{i}.V1T);
            end
            if n1~=n2
                tmp4 = iHW{i}.tDsch13.*iHW{i}.U1T*tmp{i}*iHW{i}.V3;
                H13 =  iHW{i}.U1*(tmp4*iHW{i}.V3T);
            else
                H13 = zeros(size(H21));
            end
            out{i,1} = H11 + H12 + H21 + H13;
        end
    else % user deined
        if iscell(d)
        out{i,1} = pblk.DPhi(iHW{i},d{i});
        else
        out{i,1} = pblk.DPhi(iHW{i},d);
        end
    end
end
end