%%***********************************************************************
%% AXfun: compute AX(k) = <Ak,X>, k = 1:m
%%
%% AX = AXfun(K,At,X);
%%***********************************************************************

function AX = AXfun(K,At,X)

if isa(At, 'MatCell')
   m = size(At{1},2);
   AX = zeros(m,1);
   for p = 1:length(K)
      cone = K{p};
      if strcmp(cone.type,'s')
         if (isempty(At{p}))
            AX = [];
            AXtmp = [];
         elseif length(cone.size) == 1
            AXtmp = (mexsvec(cone.blk, X{p})'*At{p})';
         else
            AXtmp = (mexsvec(cone.blk, sparse(X{p}))'*At{p})';
         end
      else
         if (isempty(At{p}))
            AX = [];
            AXtmp = [];
         else
            AXtmp = (X{p}'*At{p})';
         end
      end
      AX = AX + AXtmp;
   end
else
   if isa(K, 'Cone')
      cone = K{1};
   else
      assert (isa(K, 'BasicCone'));
      cone = K;
   end
   if strcmp(cone.type,'s')
      if (isempty(At))
         AX =[];
      elseif (length(cone.size)==1)
         AX = (mexsvec(cone.blk, X)'*At)';
      else
         AX = (mexsvec(cone.blk, sparse(X))'*At)';
      end
   else
      if (isempty(At))
         AX =[];
      else
         AX = (X'*At)';
      end
   end
end
%%*********************************************************
