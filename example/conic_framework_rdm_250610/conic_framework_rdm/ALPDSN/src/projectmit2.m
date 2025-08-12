%%
% usage
% [Y, par] = project(K, X);
% [Y, par] = project(K, X, L, U);  % Not recommend. In this case, this Jacobian is computed at the projection onto K hence not accurate .
% [Y, par] = project(K, X);
%%
%ouput
%  Y: MatCell, the projection of X onto the area
%  par: ProjJac, the Jacobian of the projection operator. par is the shortname of partial.
%      detailed decsription of each field of par is in ProjJac.m
%%************************************************************************

function [Y, par] = projectmit2(K, X, L, U)

if nargin < 3
    L = [];
end

if nargin < 4
    U = [];
end

tol = 1e-15; % % must be small as it will affect gap
addtol = 1e-6;
Y = MatCell(length(X));

n_sub_cones = sum(cellfun(@(cone) strcmp(cone.type, 's') * sum(cone.size) + (1 - strcmp(cone.type, 's')) * 1, K));
if nargout > 1
    par = ProjJac(n_sub_cones);
end
sumlen = 0;
%%
for p = 1:length(K)
    cone = K{p};

    if strcmp(cone.type, 's')
        n = sum(cone.size);
        Prow_index1 = [];
        Prow_index2 = [];
        Pcol_index1 = [];
        Pcol_index2 = [];
        %         P2row_index1 = [];
        %         P2row_index2 = [];
        %         P2col_index1 = [];
        %         P2col_index2 = [];
        Pvalue1 = [];
        Pvalue2 = [];
        Srow_index1 = [];
        Srow_index2 = [];
        Scol_index1 = [];
        Scol_index2 = [];
        %         S2row_index1 = [];
        %         S2row_index2 = [];
        %         S2col_index1 = [];
        %         S2col_index2 = [];
        Svalue1 = [];
        Svalue2 = [];
        cumsum_size = [0 cumsum(cone.size)] ;
        cumsumsquare = (cone.size.^2 );
        cumsumsquare = [0 cumsum(cumsumsquare )];
        n2 = sum(cone.size.^2);
        row_index = zeros(n2,1);
        col_index = zeros(n2,1);
        value_index = zeros(n2,1);
        if nargout == 1
            %                 if (length(cone.size) ==  1) || (sum(cone.size) < 100)
            %                     [Yp] = project_sdp(X{p}, tol);
            nsize = sum(cone.size);
%             Yp = sparse(nsize,nsize);
            for jj = 1: length(cone.size)
                cj1 = cumsum_size(jj)+1;
                cj2 = cumsum_size(jj+1);
                cjs1 = cumsumsquare(jj) + 1;
                cjs2 = cumsumsquare(jj+1);
                indx = cumsum_size(jj)+1:cumsum_size(jj+1);
                [aa] = project_sdp(X{p}(indx,indx), tol);
                [row, col] = meshgrid(cj1:cj2, cj1:cj2) ;
                row_index(cjs1:cjs2 ) = row(:);
                col_index(cjs1:cjs2) = col(:);
                value_index(cjs1:cjs2) =  aa(:);
            end

            Yp = sparse(row_index,col_index,value_index,n,n);
        else
            %                 if (length(cone.size) == 1) || (sum(cone.size) < 100)
            %                     [Yp, P, Dsch2, dd, posidx] = project_sdp(X{p}, tol);
            %             nsize = sum(cone.size);
            %             Yp = sparse(nsize,nsize);
            %             par.Q{p} = sparse(n,n);
            %             par.Sigma{p} = sparse(n,n);
            tolrkpso = 1;
            tolrkneg = 1;
            cone2 = cone;
            for jj = 1: length(cone.size)
                indx = cumsum_size(jj)+1:cumsum_size(jj+1);
                cj1 = cumsum_size(jj)+1;
                cj2 = cumsum_size(jj+1);
                cjs1 = cumsumsquare(jj) + 1;
                cjs2 = cumsumsquare(jj+1);
                indx = cumsum_size(jj)+1:cumsum_size(jj+1);
                [aa, P, Dsch2, dd, posidx] = project_sdp(X{p}(indx,indx), tol);
                [row, col] = meshgrid(cj1:cj2, cj1:cj2) ;
                row_index(cjs1:cjs2 ) = row(:);
                col_index(cjs1:cjs2) = col(:);
                value_index(cjs1:cjs2) =  aa(:);
                Dsch2 = max(addtol, Dsch2);
                Dsch1 = 1;
                shift = 0;
                if ~isempty(posidx)
                    P1 = P(:, posidx);
                    P2 = P(:, posidx(end)+1:cone.size(jj));
                else
                    P1 = [];
                    P2 = P;
                end
                cj1 = cumsum_size(jj)+1;
                cj2 = cumsum_size(jj+1);
                cjs1 = cumsum_size(jj) + 1;
                cjs2 = cumsum_size(jj+1);
                rkpos = length(posidx);

                %                 cone2.size(jj) = rkpos;
                %                 rkneg = cone.size(jj) - length(posidx);
                %                 [P2row1, P2col1] = ndgrid(cj1:cj2, tolrkpso:tolrkpso+rkpos-1) ;
                %                 [P2row2, P2col2] = ndgrid(cj1:cj2, tolrkneg:tolrkneg+rkneg-1) ;
                %                 [S2row1, S2col1] = ndgrid(tolrkpso:tolrkpso+rkpos-1, tolrkpso:tolrkpso+rkpos-1) ;
                %                 [S2row2, S2col2] = ndgrid(tolrkneg:tolrkneg+rkneg-1, tolrkneg:tolrkneg+rkneg-1) ;
                %                 tolrkpso = tolrkpso + rkpos;
                %                 tolrkneg = tolrkneg + rkneg;
                %
                %                 P2row_index1 = [P2row_index1; P2row1(:)];
                %                 P2row_index2 = [P2row_index2; P2row2(:)];
                %                 P2col_index1 = [P2col_index1; P2col1(:)];
                %                 P2col_index2 = [P2col_index2; P2col2(:)];
                %                 Pvalue1 = [Pvalue1; P1(:)];
                %                 Pvalue2 = [Pvalue2; P2(:)];
                %                 Srow_index1 = [Srow_index1; S2row1(:)];
                %                 Srow_index2 = [Srow_index2; S2row2(:)];
                %                 Scol_index1 = [Scol_index1; S2col1(:)];
                %                 Scol_index2 = [Scol_index2; S2col2(:)];

                [Prow1, Pcol1] = ndgrid(cj1:cj2, cj1:cj1+rkpos-1) ;
                [Prow2, Pcol2] = ndgrid(cj1:cj2, cj1+rkpos:cj2) ;
                [Srow1, Scol1] = ndgrid(cj1:cj1+rkpos-1, cj1:cj1+rkpos-1);
                [Srow2, Scol2] = ndgrid(cj1:cj1+rkpos-1, cj1+rkpos:cj2) ;
                Prow_index1 = [Prow_index1; Prow1(:)];
                Prow_index2 = [Prow_index2; Prow2(:)];
                Pcol_index1 = [Pcol_index1; Pcol1(:)];
                Pcol_index2 = [Pcol_index2; Pcol2(:)];
                Pvalue1 = [Pvalue1; P1(:)];
                Pvalue2 = [Pvalue2; P2(:)];
                Srow_index1 = [Srow_index1; Srow1(:)];
                Srow_index2 = [Srow_index2; Srow2(:)];
                Scol_index1 = [Scol_index1; Scol1(:)];
                Scol_index2 = [Scol_index2; Scol2(:)];
                if rkpos > 0
                    Svalue1 = [Svalue1; Dsch1*ones(length(Srow1(:)),1)];
                end

                Svalue2 = [Svalue2; Dsch2(:)];
                %                 par.Q{p}(indx,indx) = [P1,P2];
                %                 par.Sigma{p}(indx,indx) = [ ones(length(posidx),length(posidx)), Dsch2;Dsch2', sparse(length(setdiff(1:cone.size(jj), posidx)),length(setdiff(1:cone.size(jj), posidx)) ) ];
                par.P1{jj+sumlen} = P1;
                par.P1t{jj+sumlen} = P1';
                par.P2{jj+sumlen} = P2;
                par.P2t{jj+sumlen} = P2';
                par.dd{jj+sumlen} = dd;
                par.posidx{jj+sumlen} = posidx;
                par.Dsch2{jj+sumlen} = Dsch2;
                %             par.Dsch2t{p} = Dsch2';
                par.Dsch1{jj+sumlen} = Dsch1;
            end
            Yp = sparse(row_index,col_index,value_index,n,n);
            par.shift{p} = shift;
            par.Q1{p} = sparse(Prow_index1,Pcol_index1,Pvalue1,n,n);
            par.Q2{p} = sparse(Prow_index2,Pcol_index2,Pvalue2,n,n);
            par.Sigma1{p} = sparse(Srow_index1,Scol_index1,Svalue1,n,n);
            par.Sigma2{p} = sparse(Srow_index2,Scol_index2,Svalue2,n,n);
            %                 par.cone = cone2;

            %                 else
            %                     [Yp, P, Dsch2, dd, posidx] = project_sdp(X{p}, tol);
            %                 end



        end

    elseif strcmp(cone.type, 'q')

        if nargout == 1
            [Yp] = mexprojection_cone_q(full(X{p}), cone.size);
        else
            [Yp, dd, Dsch1, Dsch2, P1, P2, shift] = mexprojection_cone_q(full(X{p}), cone.size);

            par.P1{sumlen+1} = P1;
            par.P1t{sumlen+1} = P1';
            par.P2{sumlen+1} = P2;
            par.P2t{sumlen+1} = P2';
            par.dd{sumlen+1} = dd;
            par.posidx{sumlen+1} = find(dd > tol);
            par.Dsch1{sumlen+1} = Dsch1;
            par.Dsch2{sumlen+1} = Dsch2;
            par.shift{p} = shift;

            ncols = length(cone.size);
            nrows = sum(cone.size);
            Prow_index = 1:nrows;
            Pcol_index = repelem(1:ncols, cone.size);
            Pvalue1 = P1(:);
            Pvalue2 = P2(:);
            par.Q1{p} = sparse(Prow_index, Pcol_index, Pvalue1, nrows, ncols);
            par.Q2{p} = sparse(Prow_index, Pcol_index, Pvalue2, nrows, ncols);
            par.Sigma1{p} = sparse([]);
            par.Sigma2{p} = sparse([]);
        end

    elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'b2l')
        n = sum(cone.size);
        Yp = zeros(n, 1);
        posidx = find(X{p} > tol);

        if ~isempty(posidx)
            Yp(posidx) = abs(X{p}(posidx));
        end

        if nargout > 1
            P1 = []; P2 = [];
            %Dsch2 = addtol*ones(n,1);  % purtubation
            Dsch2 = zeros(n, 1);
            Dsch1 = [];
            posidx = find(X{p} > tol);

            if ~isempty(posidx)
                Dsch2(posidx) = ones(length(posidx), 1);
            end

            dd = X{p};
            shift = 0;
            par.P1{sumlen+1} = P1;
            par.P1t{sumlen+1} = P1';
            par.P2{sumlen+1} = P2;
            par.P2t{sumlen+1} = P2';
            par.dd{sumlen+1} = dd;
            par.posidx{sumlen+1} = posidx;
            par.Dsch2{sumlen+1} = Dsch2;
            %             par.Dsch2t{p} = Dsch2';
            par.Dsch1{sumlen+1} = Dsch1;
            par.shift{p} = shift;
            par.Sigma1{p} = sparse([]);
            par.Sigma2{p} = sparse([]);
        end

    elseif strcmp(cone.type, 'u')
        Yp = X{p};

        if nargout > 1
            P1 = []; P2 = [];
            %%***** perturbation *****
            Dsch2 = ones(cone.size, 1);
            Dsch1 = [];
            posidx = [];
            dd = [];
            shift = 0;
            par.P1{sumlen+1} = P1;
            par.P1t{sumlen+1} = P1';
            par.P2{sumlen+1} = P2;
            par.P2t{sumlen+1} = P2';
            par.dd{sumlen+1} = dd;
            par.posidx{sumlen+1} = posidx;
            par.Dsch2{sumlen+1} = Dsch2;
            %             par.Dsch2t{p} = Dsch2';
            par.Dsch1{sumlen+1} = Dsch1;
            par.shift{p} = shift;
            par.Sigma1{p} = sparse([]);
            par.Sigma2{p} = sparse([]);
        end

    elseif strcmp(cone.type, 'b')
        n = sum(cone.size);
        Yp = X{p};
        loweridx = find(X{p} < cone.params(:,1) + tol);
        upperidx = find(X{p} > cone.params(:,2) - tol);
        Yp(loweridx) = cone.params(loweridx, 1);
        Yp(upperidx) = cone.params(upperidx, 2);

        if nargout > 1
            P1 = []; P2 = [];
            %Dsch2 = addtol*ones(n,1);  % purtubation
            Dsch2 = zeros(n, 1);
            Dsch1 = [];
            posidx = setdiff(1:n, union(loweridx, upperidx));

            if ~isempty(posidx)
                Dsch2(posidx) = ones(length(posidx), 1);
            end

            dd = X{p};
            shift = 0;
            par.P1{sumlen+1} = P1;
            par.P1t{sumlen+1} = P1';
            par.P2{sumlen+1} = P2;
            par.P2t{sumlen+1} = P2';
            par.dd{sumlen+1} = dd;
            par.posidx{sumlen+1} = posidx;
            par.Dsch2{sumlen+1} = Dsch2;
            %             par.Dsch2t{p} = Dsch2';
            par.Dsch1{sumlen+1} = Dsch1;
            par.shift{p} = shift;
            par.Sigma1{p} = sparse([]);
            par.Sigma2{p} = sparse([]);
        end
    end
    sumlen = sumlen + strcmp(cone.type, 's') * length(cone.size) + (1 - strcmp(cone.type, 's')) * 1;
    if ~isempty(L)
        Yp = max(Yp, L{p});
    end

    if ~isempty(U)
        Yp = min(Yp, U{p});
    end

    Y{p} = Yp;

    %         if nargout > 1
    %             par.P1{sumlen} = P1;
    %             par.P1t{sumlen} = P1';
    %             par.P2{sumlen} = P2;
    %             par.P2t{sumlen} = P2';
    %             par.dd{sumlen} = dd;
    %             par.posidx{sumlen} = posidx;
    %             par.Dsch2{sumlen} = Dsch2;
    % %             par.Dsch2t{p} = Dsch2';
    %             par.Dsch1{sumlen} = Dsch1;
    %             par.shift{sumlen} = shift;
    %         end

end

end

%%***************************************************************************

function [Y, V, Dsch2, d, posidx] = project_sdp(X, tol)
n = length(X);
%     exist_mexeig = exist('mexeig', 'file');
exist_mexeig = 0;
X(abs(X) < 1e-14) = 0;
X = 0.5 * (X + X');

if (exist_mexeig == 3)
    %[V,D] = mexeig(full(X));
    [V, D] = eig(full(X));
else
    [V, D] = eig(full(X));
end

d = diag(D);
[d, idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:, idx);
posidx = find(d > tol);

if isempty(posidx)
    Y = sparse(n, n);
    Dsch2 = [];
elseif (length(posidx) == n)
    Y = X;
    Dsch2 = [];
else
    r = length(posidx); s = n - r;
    negidx = [r + 1:n];
    dp = abs(d(posidx));
    dn = abs(d(negidx));
    Vtmp = V(:, posidx) * diag(sqrt(dp));
    Y = Vtmp * Vtmp';
    Y = 0.5 * (Y + Y');
    Dsch2 = (dp * ones(1, s)) ./ (dp * ones(1, s) + ones(r, 1) * dn');
end

end

function [Y, V, Dsch2, d, posidx] = project_sdp_multiblock(X, tol)

n = length(X);
exist_mexeig = exist('mexeig', 'file');
X(abs(X) < 1e-14) = 0;

X = 0.5 * (X + X');
pm = symamd(X);
Xperm = X(pm, pm);
[t, ~] = etree(Xperm);
idx0 = find(t == 0);
len0 = length(idx0);
Dsub = zeros(n, 1);
d = zeros(n, 1);
Vsub = cell(n, 1);
offset = 1;

for it = 1:len0
    idx = offset:idx0(it);

    if (exist_mexeig == 3)
        %[Vsub0,Dsub0] = mexeig(full(Xperm(idx,idx)));
        [Vsub0, Dsub0] = eig(full(Xperm(idx, idx)));
    else
        [Vsub0, Dsub0] = eig(full(Xperm(idx, idx)));
    end

    Dsub0 = diag(Dsub0);
    Vsub{it} = Vsub0;
    Dsub(idx) = Dsub0;
    offset = idx0(it) + 1;
end

V = blkdiag(Vsub{:});
V(pm, pm) = V;
d(pm) = Dsub;

[d, idx] = sort(real(d));
idx = idx(n:-1:1); d = d(n:-1:1);
V = V(:, idx);
posidx = find(d > tol);

if isempty(posidx)
    Y = sparse(n, n);
    Dsch2 = [];
elseif (length(posidx) == n)
    Y = X;
    Dsch2 = [];
else
    r = length(posidx); s = n - r;
    negidx = [r + 1:n];
    dp = abs(d(posidx));
    dn = abs(d(negidx));
    Vtmp = V(:, posidx) * diag(sqrt(dp));
    Y = Vtmp * Vtmp';
    Y = 0.5 * (Y + Y');
    Dsch2 = (dp * ones(1, s)) ./ (dp * ones(1, s) + ones(r, 1) * dn');
end

end