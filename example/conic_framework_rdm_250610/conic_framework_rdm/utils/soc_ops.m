function out = soc_ops(X, operand, cone_size)
% % some useful operations for varibles in a second order cone
% or concatenateed by several second order cones
% example:
%   size(x) ==  [14, 1]
%   cone_size = [2 3 4 5]
%   here x is concatenated by four cones
% if cone_size is missing, then X is treated as a single second order cone
    assert(nargin >= 2);
    assert(size(X, 2) == 1);

    assert( strcmp(operand, 'norm') || strcmp(operand, 'det') || strcmp(operand, 'sqrtdet') || ...
        strcmp(operand, 'inv') || strcmp(operand, 'Q') || strcmp(operand, 'Qinv') || strcmp(operand, 'norm_xbar'))

    % single block operations
    if nargin == 2 
        if strcmp(operand, 'det')
            out = X(1)^2 - norm(X(2:end))^ 2;
        elseif strcmp(operand, 'norm')
            out = norm(X);
        elseif strcmp(operand, 'norm_xbar')
            out = norm(X(2:end));
        elseif strcmp(operand, 'sqrtdet')
            out = sqrt(X(1)^2 - sum(X(2:end) .^ 2));
        elseif strcmp(operand, 'inv')
            out = 1 / soc_ops(X, 'det') * [X(1); -X(2:end)];
        elseif strcmp(operand, 'Q')
            n = size(X, 1);
            xbar = X(2:end); 

            detx = soc_ops(X, 'det');
            out = [norm(X)^2, 2 * X(1) * xbar';
                2 *X(1) * xbar, detx * eye(n-1) + 2 * xbar * xbar'];
        elseif strcmp(operand, 'Qinv')
            n = size(X, 1);
            xbar = X(2:end); 
            detx = soc_ops(X, 'det');
            out = 1 / ( detx ^ 2) * ...
                [norm(X)^2, - 2 *X(1) * xbar';
                - 2 * X(1) * xbar, detx * eye(n-1) + 2 * xbar * xbar'];    
        end
        return
    end

    % multi block operations
    if iscell(cone_size)
        assert(size(cone_size, 1) == 1)
        cone_size = cone_size{2};
    end

    if strcmp(operand, 'norm') || strcmp(operand, 'det') || strcmp(operand, 'sqrtdet') || ...
        strcmp(operand, 'inv') || strcmp(operand, 'norm_xbar')
        out = mexcone_q_ops(X, operand, cone_size);
    elseif strcmp(operand, 'Q') 
        n = size(X, 1);
        detX = mexcone_q_ops(X, 'det', cone_size);
        ind_head = cumsum(cone_size) - cone_size + 1;
        temp = ones(n, 1);
        temp(ind_head) = -1;
        X1 = blk_spdiag(X, cone_size);
        out = spdiags(repelem(detX, cone_size, 1) .* temp, 0, n, n) + 2 * (X1 * X1');
    elseif strcmp(operand, 'Qinv')
        n = size(X, 1);
        detX = mexcone_q_ops(X, 'det', cone_size);
        Xinv = mexcone_q_ops(X, 'inv', cone_size);
        ind_head = cumsum(cone_size) - cone_size + 1;
        temp = ones(n, 1);
        temp(ind_head) = -1;
        Xinv1 = blk_spdiag(Xinv, cone_size);
        out = spdiags(repelem(1 ./ detX, cone_size, 1) .* temp, 0, n, n) + 2 * (Xinv1 * Xinv1');
    end
        
end