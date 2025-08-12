function C = dot_ssn(A, B)
% if isempty(A) || isempty(B) || ( ~nnz(B) ) || ( ~nnz(A) )
if isempty(A) || isempty(B) || isequal(B,0) || isequal(A,0)
    C = 0;
elseif iscell(A) && iscell(B)
    C = full(sum(cellfun(@(x, y) computeinprod(x,y), A, B)));
elseif ~iscell(A) && ~iscell(B)
    C = A(:)'*B(:);
elseif ~iscell(A) && iscell(B)
    C = 0;
    for i = 1:length(B)
        C = C + A(:)'*B{i}(:);
    end
elseif iscell(A) && ~iscell(B)
    C = 0;
    for i = 1:length(A)
        C = C + A{i}(:)'*B(:);
    end
end
end
function n = computeinprod(x,y)
if issparse(x)
    if isequal(x,0)
        n = 0;
    else
    n = sum(x .* y,'all');
    end
else
    n = sum(x(:)' * y(:));
end
end