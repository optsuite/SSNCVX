function out = blk_average(A, block_size, dim, broadcast)
    %% out = blk_average(A, block_size, dim, broadcast)
    % block average for a matrix along the given dimension
    %
    % input arguments 
    %     A: a sparse or dense matrix of size [m, n]
    %     block_size: block informations, a vector with k integer elements, i.e., block_size = [q_1, ..., q_k]
    %          also support SDPT3 type format block_size = {'q', [q_1, ..., q_k]}
    %     dim: 1 or 2
    %         1: compute the average vertly, i.e., 
    %             it should hold that q_1 + ... + q_k == m
    %             it will return a sparse matrix U of size [k, n] such that
    %             the j-th row of U equals the average of the q_1 + ... + q_j + 1 -th to q1 + ...+q_{j+1} -th 
    %             rows of A.
    %         2 : compute the average horizonally, i.e, 
    %             it should hold that q_1 + ... + q_k == n
    %             it will return a sparse matrix U of size [m, k] such that
    %             the j-th column of U equals the average of the q_1 + ... + q_j + 1 -th to q1 + ...+q_{j+1} -th 
    %             columns of A.
    %         defaulted dim = 1
    %     broadcast: bool, decide whether broadcast such that the final output has the same size as input
    %         false: if dim == 1 then size(output) == [k, n]
    %                if dim == 2 then size(output) == [m, k]
    %         true:  size(output) == size(A)
    %         defaulted broadcast = false
    % output
    %     a matrix with its sparcity dependent on A and its size dependent on the input parameters
    %
    % example
    %     >> A = eye(10); block_size = [1 2 3 4];
    %     >> blk_average(A, block_size, 1) 
    %
    %       ans =
    %           1.0000    0         0         0         0         0         0         0         0         0
    %           0    0.5000    0.5000         0         0         0         0         0         0         0
    %           0         0         0    0.3333    0.3333    0.3333         0         0         0         0
    %           0         0         0         0         0         0    0.2500    0.2500    0.2500    0.2500
    %
    %     >> blk_average(A, block_size, 2, true)
    %
    %       ans =
    %
    %           1.0000    0         0         0         0         0         0         0         0         0
    %           0    0.5000    0.5000         0         0         0         0         0         0         0
    %           0    0.5000    0.5000         0         0         0         0         0         0         0
    %           0         0         0    0.3333    0.3333    0.3333         0         0         0         0
    %           0         0         0    0.3333    0.3333    0.3333         0         0         0         0
    %           0         0         0    0.3333    0.3333    0.3333         0         0         0         0
    %           0         0         0         0         0         0    0.2500    0.2500    0.2500    0.2500
    %           0         0         0         0         0         0    0.2500    0.2500    0.2500    0.2500
    %           0         0         0         0         0         0    0.2500    0.2500    0.2500    0.2500
    %           0         0         0         0         0         0    0.2500    0.2500    0.2500    0.2500


    if nargin < 3
        dim = 1;
    end
    if nargin < 4
        broadcast = false;
    end

    if iscell(block_size)
        block_size = cell2mat(block_size(:, 2)');
    end

    assert(dim == 1 || dim == 2)
    [m, n] = size(A);
    if dim == 1
        assert(sum(block_size, 'all') == m);
        coeff = repelem(reshape(1 ./ block_size, 1, length(block_size)), 1, block_size);
        out = blk_spdiag(coeff, block_size) * A;
        if broadcast
            out = repelem(out, block_size, 1);
        end
    else
        assert(sum(block_size, 'all') == n);
        coeff = repelem(reshape(1 ./ block_size, length(block_size), 1), block_size, 1);
        out = A * blk_spdiag(coeff, block_size);
        if broadcast
            out = repelem(out, 1, block_size);
        end
    end

end