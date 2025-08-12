%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-11 21:45:32
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% Now we support only affine and affine_inv 
% Given a ProjJac object pj. pj can be seen as a mappping.
% affine means 
%   out = scale * pj  + shift * identity
% affine_inv means
%   out = (scale * pj  + shift * identity)^{-1}
% out is also a ProjJac object

function out = ProjJac_ops(pj, K, operand, scale, shift)
    assert(strcmp(operand, 'affine') || strcmp(operand, 'affine_inv'))
    
    % check pj and K are compatible
    % to be implemented

    out = ProjJac(length(K));

    % copy spectral decomposition 
    out.dd = pj.dd;
    out.P1 = pj.P1;
    out.P2 = pj.P2;
    out.posidx = pj.posidx;

    % transformation on Dsch1, Dsch2 and shift
    if strcmp(operand, 'affine')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 's')
                out.Dsch1{p} = scale * pj.Dsch1{p};
                out.Dsch2{p} = scale * pj.Dsch2{p};
                out.Dsch2t{p} = out.Dsch2{p}';
                out.shift{p} = scale * pj.shift{p} + shift;
            elseif strcmp(cone.type, 'q')    
                out.Dsch1{p} = scale * pj.Dsch1{p};
                out.Dsch2{p} = scale * pj.Dsch2{p};
                out.shift{p} = scale * pj.shift{p} + shift;
            elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
                out.Dsch2{p} = scale * pj.Dsch2{p} + shift;
            end
        end
    elseif strcmp(operand, 'affine_inv')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 's')
                out.shift{p} = 1 ./ (scale * pj.shift{p}  + shift);
                out.Dsch1{p} = 1 ./ (scale * ( pj.Dsch1{p} + pj.shift{p} ) + shift) - out.shift{p};
                out.Dsch2{p} = 1 ./ (scale * ( pj.Dsch2{p} + pj.shift{p} ) + shift) - out.shift{p};
                out.Dsch2t{p} = out.Dsch2{p}';
            elseif strcmp(cone.type, 'q')    
                out.Dsch1{p} = scale * pj.Dsch1{p};
                out.Dsch2{p} = scale * pj.Dsch2{p};
                out.shift{p} = scale * pj.shift{p} + shift;  % note that pj.shift{p} is a vector
                idx1 = (pj.dd{p}(:, 2) >= 0);   % both eigenvalues are positive
                idx2 = (pj.dd{p}(:, 1) > 0 & pj.dd{p}(:, 2) < 0); % one eigenvalue is positive
                idx3 = (pj.dd{p}(:, 1) <= 0); % both eigenvalues are negative
                out.Dsch1{p}(idx1) = 0;
                out.Dsch2{p}(idx1) = 0;
                out.shift{p}(idx1) = 1 ./ out.shift{p}(idx1);
                out.Dsch1{p}(idx2) = - 1 ./ out.shift{p}(idx2) .* out.Dsch1{p}(idx2) ./ ( out.Dsch1{p}(idx2) + out.shift{p}(idx2) );
                out.Dsch2{p}(idx2) = - 1 ./ out.shift{p}(idx2) .* out.Dsch2{p}(idx2) ./ ( out.Dsch2{p}(idx2) + out.shift{p}(idx2) );
                out.shift{p}(idx2) = 1 ./ out.shift{p}(idx2);
                out.Dsch1{p}(idx3) = 0;
                out.Dsch2{p}(idx3) = 0;
                out.shift{p}(idx3) = 1 ./ out.shift{p}(idx3);
            elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
                out.Dsch2{p} = 1 ./ (scale * pj.Dsch2{p} + shift);
                out.Dsch2t{p} = out.Dsch2{p}';
            end
        end
    end
    % copy precond
    out.precond = pj.precond;

    %% copy all other fileds from pj
    pj_fields = fieldnames(pj);
    for i = 1: length(pj_fields)
        field = pj_fields{i};
        if ~ isfield(out, field)
            out.(field) = pj.(field);
        end
    end
    

end