%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-09 16:57:12
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function print_std_model_info(fid, model)

    if nargin < 2
        fid = 1;
    end

    fprintf(fid, "[Matrices A]\n");
    At = model.At;
    m = 0; n = 0;
    normA = 0;
    for p = 1: length(At)
        normAp = norm(At{p}, 'fro');
        normA = normA + normAp ^2;
        fprintf(fid, "\tblock %i: %i x %i, norm  = %.4e, sparcity = %.4e\n", ...
            p, size(At{p}, 2), size(At{p}, 1), normAp, nnz(At{p}) / numel(At{p}));
        if p == 1
            m = size(At{p}, 2);
        end
        n = n + size(At{p}, 1);
    end
    normA = sqrt(normA);
    fprintf(fid, "\ttotal: %i x %i, norm = %.4e, sparcity = %.4e\n", ...
        m, n, normA, sparcity(model.At));

    K = model.K;
    fprintf(fid, "[Cone]\n");
    print_tol = 10;
    for p = 1: length(K)
        cone = K{p};
        fprintf(fid, "\tblock %i: '%s' ", p, cone.type);
        if length(cone.size) <= print_tol  % print all if not too many
            fprintf(fid, "[ ");
            for j = 1: length(cone.size)
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "]\n");
        else    % print first and last few if too many
            fprintf(fid, "[ ");
            n_head = ceil(print_tol / 2); n_tail = print_tol - n_head;
            for j = 1: n_head
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "... ");
            for j = length(cone.size) - n_tail + 1: length(cone.size)
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "]\n");
            fprintf(fid, "\t\t(Cartesian product of %i cones with sizes ranging from %i to %i)\n", ...
                numel(cone.size), min(cone.size), max(cone.size));
        end  
    end 

    
end