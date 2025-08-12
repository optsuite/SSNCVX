%%-------------------------------------------------------------
%% projection on the box [L,U]
%%-------------------------------------------------------------

function [PX,PK] = projBox(X,L,U)
tiny = 1e-14;
PX = X;
plength = length(X);
if iscell(X)
    for i = 1:plength
        if ~isempty(L)
            if iscell(L)
            PX{i,1} = max(X{i},L{i}+tiny);
            else
            PX{i,1} = max(X{i},L+tiny);
            end
        end
        if ~isempty(U)
            if iscell(U)
            PX{i,1} = min(PX{i},U{i}-tiny);
            else
            PX{i,1} = min(PX{i},U-tiny);
            end
        end
        if nargout == 2
            PK{i,1} = abs(PX{i} - X{i}) < 1e-16;
        end
    end
else
    if ~isempty(L)
        PX = max(X,L+tiny);
    end
    if ~isempty(U)
        PX = min(PX,U-tiny);
    end
    if nargout == 2
        PK = abs(PX - X) < 1e-16;
    end
end
%%-------------------------------------------------------------

