
function [F,S,x4, params] = compute_gradient(y,z,r,v,x1,x2,x3,x4,params)
%% compute_gradient: compuate the gradient of the saddle point problem
%%
%% Copyright (c) 2025 by
%% Zhanwang Deng, Tao Wei, Jirui Ma, Zaiwen Wen
%%
F.res = 0;
if params.Aboxflag == 1
    if params.Axeqb
        tmp1 = params.b;
    else
        tmp1 = x1.var - params.sigma*y.var;
        [tmp1, D1] = params.P2box(tmp1);
        params.D.D1 = D1;
    end
else
    tmp1 = [];
end

if params.fnonsmooth
    tmp2 = x2.var - params.sigma*z.var;
    for i =1: params.nblock
        if isfield(params.f{i},'shift')
            tmp2{i} = tmp2{i}  - params.sigma*params.f{i}.shift;
        end
    end
else
    tmp2 = 0;
end

if params.boxflag == 1
    tmp3k = x3.var - params.sigma*r.var;
    [tmp3, D3]= params.P1box(tmp3k);
    params.D.D3 = D3;
end
if params.Qflag == 1
    Qv = params.Qcmap(v.var);
else
    Qv = 0;
end


tmp40 = x4.var + params.sigma*(params.ATmap(y.var) + params.BTmap(z.var) - Qv + params.idmap( r.var) - params.C );
for i = 1:params.nblock
    if isfield(params.pblk{i},'shift')
        tmp40{i} = tmp40{i} - params.pblk{i}.shift;
    end
end
% end
%%

if isfield(params,'f') && params.fnonsmooth
    [tmp2, D2] = params.prox_f(tmp2,params);
    params.D.D2 = D2;
elseif isfield(params,'f') && ~params.fnonsmooth
    tmp2 = params.dual_gradientf(z.var);
    for i = 1:params.nblock
        if strcmp(params.f{i}.type,'square')
            params.D.g2 = 0.5/params.f{i}.coefficient;
        elseif  strcmp(params.f{i}.type,'exp')
            params.D.g2 = 1./max(-z.var{i},1e-7);
        elseif  strcmp(params.f{i}.type,'log')
            params.D.g2 = params.f{i}.coefficient/(z.var{i}.^2);
        elseif  strcmp(params.f{i}.type,'logdet')
            params.D.g2 = inv(z.var{i});
        else %
            params.D.g2 = params.D2(z.var{i});
        end
    end
end



[tmp4, D4] = params.prox_p(tmp40,params);
S.var = (tmp4 - tmp40)/params.sigma;
for i = 1:params.nblock
if isfield(params.pblk{i},'shift')
    tmp4{i} = tmp4{i} + params.pblk{i}.shift;
    end
end

params.D.D4 = D4;


if params.Aboxflag == 1
    F.Fy = params.Amap(tmp4) - tmp1;
    F.Fyres = norm(F.Fy);
    F.res = F.res + F.Fyres;
else
    F.Fy = 0;
    F.Fyres = 0;
end

if isfield(params,'f') && ~isempty(params.f)
    F.Fz = params.Bmap(tmp4) - tmp2;
    for p = 1:params.nblock
        if isfield(params.f{p},'shift') && params.fnonsmooth
            F.Fz{p,1} = F.Fz{p,1} - params.f{p}.shift ;
        end

    end
    F.Fzres = norm(F.Fz);
    F.res = F.res + F.Fzres;
else
    F.Fz = zeros_like(z.var);
    F.Fzres = 0;
end

if params.boxflag == 1
    F.Fr = tmp4 -tmp3;
    F.Frres = norm(F.Fr);
    F.res = F.res + F.Frres;
else
    % F.Fr = {zeros(1,1)};
    F.Fr = zeros_like(r.var);
    F.Frres = 0;
end

if params.Qflag == 1
    F.Fv = params.Qcmap(v.var - tmp4);
    if ~iscell(F.Fv)
        F.Fv = {F.Fv};
    end
    F.Fvres = norm(F.Fv);
    F.res = F.res + F.Fvres;
else
    F.Fv = zeros_like(v.var);
    F.Fvres = 0;
end

if params.Aboxflag == 1 && params.Axeqb == 0
    F.Fx1 = (x1.var - tmp1)/params.sigma;
    F.Fx1res = norm(F.Fx1);
    F.res = F.res + F.Fx1res;
else
    F.Fx1 = zeros_like(x1.var);
    F.Fx1res = [];
end

if params.fnonsmooth
    F.Fx2 = x2.var - tmp2;
    for i = 1:params.nblock
        if isfield(params.f{i},'shift')
            F.Fx2{i} = (F.Fx2{i}- params.f{i}.shift)/params.sigma;
        end
    end
    F.Fx2res = norm(F.Fx2);
    F.res = F.res + F.Fx2res;
else
    F.Fx2 = zeros_like(x2.var);
    F.Fx2res = 0;
end

if params.boxflag == 1
    F.Fx3 = (x3.var - tmp3)/params.sigma;
    F.Fx3res = norm(F.Fx3);
    F.res = F.res + F.Fx3res;
else
    F.Fx3 = [];
    F.Fx3res = 0;
end


F.Fx4 = (x4.var - tmp4 )/params.sigma;

F.Fx4res = norm(F.Fx4);
F.res = F.res + F.Fx4res;





end