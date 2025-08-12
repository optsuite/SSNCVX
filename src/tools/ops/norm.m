

function out = norm(A)
    norms = cellfun(@(x) computeNorm(x), A);
    out = sqrt(sum(norms));
end

function y = computeNorm(x)

    if issparse(x)
        y = norm(x, 'fro')^2;
    else
        y = sum(x(:)' * x(:));
    end
end