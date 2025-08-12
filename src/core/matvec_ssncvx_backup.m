function [out] = matvec_ssncvx( At, par, q, params)

if (nargin < 5); AL = []; end
if isempty(AL); existAL = 0; else; existAL = 1; end
N = length(q.q1);
AL = params.Lchol;
%%
yorg1 = q.q1;
yorg2 = q.q2;
yorg3 = q.q3;
yorg4 = q.q4;

y1 = q.q1;
y2 = q.q2;
y3 = q.q3;
y4 = q.q4;
if ~AL.isidentity  && (existAL)
    if strcmp(AL.matfct_options,'chol')
        y1(AL.perm) = AL.R \ y1;
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        y1(AL.perm) = mexbwsolve(AL.Rt,y1);
    end
end

%% perform A * (par.sigma * par.Dsch + par.epsilon * identity) * A' * y
Ax1 = zeros(N,1);
Ax2 = MatCell(params.z0); %modified
Ax3 = MatCell(params.r0);
Ax4 = MatCell(params.v0);

for p = 1: length(params.pblk)
    cone = params.pblk{p};
    n = sum(cone.size);
    if strcmp(cone.type,'s')
        tmpq4 = params.Qcmap(y4);
        rr = size(par.D4{p}.P1, 2);
        tmp = zeros(cone.size);
        % isfield
        Aty = Atyfun(params.K{p}, params.At{p}, y1) + params.BTmap(y2{p})  + params.idmap(y3) - tmpq4;
        % Aty = Aty{1};
        if (rr > 0 && rr < n)
            if (rr <= n/2)
                tmp0 = par.D4{p}.P1'*Aty;
                tmp1 = (tmp0*par.D4{p}.P1)*par.D4{p}.P1';
                tmp2 = par.D12{p}.*(tmp0*par.D4{p}.P2);
                tmp2 = tmp2*par.D4{p}.P2t;
                tmp3 = par.D4{p}.P1*(0.5*par.D11{p}*tmp1 + tmp2);
                tmp =  tmp3+tmp3';
            else
                tmp0 = par.D4{p}.P2t*Aty;
                tmp1 = (tmp0*par.D4{p}.P2)*par.D4{p}.P2t;
                tmp2 = (par.D11{p}-par.D12{p}').*(tmp0*par.D4{p}.P1);
                tmp2 = tmp2*par.D4{p}.P1t;
                tmp3 = par.D4{p}.P2*(0.5*par.D11{p}*tmp1 + tmp2);
                tmp =  par.D11{p}*Aty-tmp3-tmp3';
            end
            % Ax1 = Ax1 + par.sigma*AXfun(params.K{p}, params.At{p}, tmp);
        elseif (rr == n)
            tmp = Aty;
            % Ax1 = Ax1 + par.sigma*AXfun(params.K{p}, params.At{p}, tmp);
        end
        % Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
        % Ax3 = Ax3 + par.sigma*params.idmap(tmp);
        % Ax4 = Ax4 - params.Qmap(par.sigma*tmp );
    elseif strcmp(cone.type, 'q')
        if (~isempty(par.D12))
            Aty = At{p} * y1;
            P1blkdiag = blk_spdiag(par.P1{p}, cone.size);
            P2blkdiag = blk_spdiag(par.P2{p}, cone.size);
            tmp = repelem(par.shift{p}, cone.size, 1) .* Aty + ...
                P1blkdiag * spdiag(par.D11) * (P1blkdiag' * Aty)+ ...
                P2blkdiag * spdiag(par.D12) * (P2blkdiag' * Aty);
            Ax1 = Ax1 + par.sigma*(tmp'*At{p})';
        end
    elseif strcmp(cone.type, 'l')
        tmpq4 = params.Qcmap(y4);
        Aty = params.At{p}*y1 + params.BTmap(y2)  + params.idmap(y3) - tmpq4;
        tmp = par.lD.*Aty;
    elseif strcmp(cone.type,'u')
        tmpq4 = params.Qcmap(y4);
        Aty = params.At{p}*y1 + params.BTmap(y2)  + params.idmap(y3) - tmpq4;
        tmp = par.uD.*Aty;
    elseif strcmp(cone.type,'u2b')
        tmpq4 = params.Qcmap(y4);
        Aty = params.At{p}*y1 + params.BTmap(y2)  + params.idmap(y3) - tmpq4;
        tmp = par.b2lD.*Aty;
    elseif strcmp(cone.type,'l1') || strcmp(cone.type, 'linfty') || strcmp(cone.type, 'max') || strcmp(cone.type, 'box') || strcmp(cone.type, 'topk') || strcmp(cone.type, 'l1topk') || strcmp(cone.type, 'l1con') || strcmp(cone.type, 'linftycon')  || strcmp(cone.type, 'l1linfty') || strcmp(cone.type, 'huber')
        % Aty = Atyfun(cone, params.At{p}, y1)+ params.BTmap(y2) + params.idmap(y3) - params.Qmap(y4) ;
        tmpq4 = params.Qcmap(y4);
        Aty = params.ATmap(y1)+ params.BTmap(y2) + params.idmap(y3) - tmpq4;
        % if iscell(Aty)
        % Aty = Aty{1};
        % end
        tmp = par.tD4{p}.*Aty;
        % Ax1 = Ax1 + par.sigma*params.Amap(tmp);
        % Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
        % Ax3 = Ax3 + par.sigma*params.idmap(tmp);
        % Ax4 = Ax4 - params.Qcmap(par.sigma*tmp );
    elseif strcmp(cone.type,'l2') 
        tmpq4 = params.Qcmap(y4);
        Aty = params.ATmap(y1)+ params.BTmap(y2) + params.idmap(y3) - tmpq4;
        if par.D4{p}.type == 1
            tmp = zeros(size(Aty));
        else
            tmp = par.tD4coe1 *Aty + par.tD4coe2 *par.D4{p}.X'*Aty{p}*par.D4{p}.X;
        end
        % Ax1 = Ax1 + par.sigma*params.Amap(tmp);
        % Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
        % Ax3 = Ax3 + par.sigma*params.idmap(tmp);
        % Ax4 = Ax4 - params.Qcmap(par.sigma*tmp );
    elseif strcmp(cone.type,'l1l2')
        tmpq4 = params.Qcmap(y4);
        tmp2 = y2{1};
        Aty = params.ATmap(y1)+ params.BTmap(tmp2) + params.idmap(y3) - tmpq4;
        tmp = zeros(size(Aty));
        for j = 1:cone.size(2)
            if par.D4{1,j}.type ~= 1
                tmp(:,j) = par.tD4coe1(j) *Aty(:,j) + par.tD4coe2(j) *par.D4{p,j}.X'*Aty(:,j)*par.D4{p,j}.X;
            end
        end
        % Ax1 = Ax1 + par.sigma*params.Amap(tmp);
        % Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
        % Ax3 = Ax3 + par.sigma*params.idmap(tmp);
        % Ax4 = Ax4 - params.Qmap(par.sigma*tmp );
    elseif strcmp(cone.type,'nuclear') || strcmp(cone.type, 'l2l2')
        tmpq4 = params.Qcmap(y4);
        tmp2 = params.BTmap(y2);
        Aty = params.ATmap(y1)+ tmp2{1} + params.idmap(y3) - tmpq4;
        [n1,n2,~] = size(Aty);
        rr = par.D4{p}.trank;
        tmp = zeros(size(Aty));
        if (rr > 0 )
            tmp1 = par.D4{p}.U1T*Aty*par.D4{p}.V1;
            if rr < n1/2
                tmp2 = par.D4{p}.U1T*Aty*par.D4{p}.V2;
                tmp3 = par.D4{p}.U2T*(Aty*par.D4{p}.V1);
            else
                tmp2 = par.D4{p}.U1T*(Aty*par.D4{p}.V2);
                tmp3 = par.D4{p}.U2T*Aty*par.D4{p}.V1;
            end
            tmp11 = tmp1.* par.tDsch11 + tmp1'.* par.tHsch11;
            tmp12 = tmp2.* par.tDsch12 + tmp3'.* par.tHsch12;
            tmp21 = tmp3.* par.tDsch12'+ tmp2'.* par.tHsch12';
            H11 = par.D4{p}.U1*tmp11*par.D4{p}.V1T;
            if rr <n1/2
                H12 = par.D4{p}.U1*(tmp12*par.D4{p}.V2T);
                H21 = par.D4{p}.U2*tmp21*par.D4{p}.V1T;
            else
                H12 = par.D4{p}.U1*tmp12*par.D4{p}.V2T;
                H21 = par.D4{p}.U2*(tmp21*par.D4{p}.V1T);
            end
            if n1~=n2
                tmp4 = par.D4{p}.tDsch13.*par.D4{p}.U1T*tmp{p}*par.D4{p}.V3;
                H13 =  par.D4{p}.U1*(tmp4*par.D4{p}.V3T);
            else
                H13 = zeros(size(H21));
            end
            tmp = H11 + H12 + H21 + H13;
        end
        % Ax1 = Ax1 + par.sigma*params.Amap(tmp);
        % Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
        % Ax3 = Ax3 + par.sigma*params.idmap(tmp);
        % Ax4 = Ax4 - params.Qmap(par.sigma*tmp );
    end
    if strcmp(cone.type,'s')
        Ax1 = Ax1 + par.sigma*AXfun(params.K{p}, params.At{p}, tmp);
    elseif strcmp(cone.type,'l')
        Ax1 = Ax1+ par.sigma*(tmp'*At{p})';
    else
        Ax1 = Ax1 + par.sigma*params.Amap(tmp);
    end
    Ax2 = Ax2 + par.sigma*params.Bmap(tmp);
    Ax3 = Ax3 + par.sigma*params.idmap(tmp);
    if iscell(tmp)
        Ax4 = Ax4 - params.Qcmap(par.sigma*tmp );
    else
        Ax4 = Ax4 - params.Qmap(par.sigma*tmp );
    end
end
if (existAL) && ~AL.isidentity
    if strcmp(AL.matfct_options,'chol')
        Ax1 = AL.Rt \ Ax1(AL.perm);
    elseif strcmp(AL.matfct_options,'spcholmatlab')
        Ax1 = mexfwsolve(AL.R,Ax1(AL.perm));
    end
end



if isfield(par,'tD1')
    out.Ax1 = Ax1 + yorg1 .* (par.tD1 + par.epsilon1);
else
    out.Ax1 = Ax1 + par.epsilon1 * yorg1;
end

if isfield(par,'tD2')
    for p = 1: length(params.pblk)
        if strcmp(params.f{p}.type,'l1')  || strcmp(params.f{p}.type, 'linfty') || strcmp(params.f{p}.type, 'max') || strcmp(params.f{p}.type, 'topk') || strcmp(params.f{p}.type, 'l1topk')|| strcmp(params.f{p}.type, 'l1con') || strcmp(params.f{p}.type, 'linftycon') || strcmp(params.f{p}.type, 'box')
            out.Ax2 = Ax2 + yorg2 .* (par.tD2 + par.epsilon2);
        elseif strcmp(params.f{p}.type,'l2') || strcmp(params.f{p}.type,'l2con')
            if par.D2{p}.type == 1
                out.Ax2 = Ax2 + par.epsilon2 * yorg2   + par.tD2.coe2 *par.D2{p}.X'*Aty*par.D2{p}.X;
            else
                out.Ax2 = Ax2 + par.epsilon2 * yorg2;
            end
        elseif strcmp(params.f{p}.type,'logdet')
            out.Ax2 = Ax2 + par.epsilon2 * yorg2 +  params.f{p}.coefficient*par.tD2*yorg2*par.tD2;
        else
            out.Ax2 = Ax2 + yorg2 .*  (par.tD2 + par.epsilon2) ;
        end
    end
else
    out.Ax2 = Ax2 + yorg2 .*  (par.epsilon2) ;
end


if isfield(par,'tD3')
    out.Ax3 = Ax3 + yorg3*  par.epsilon3 + par.tD3.*yorg3 ;
else
    out.Ax3 = Ax3 + yorg3 .*   par.epsilon3;
end

if isfield(params,'Qt')
    out.Ax4 = Ax4 + yorg4 *  par.epsilon4 + tmpq4;
end

%%************************************************************************
