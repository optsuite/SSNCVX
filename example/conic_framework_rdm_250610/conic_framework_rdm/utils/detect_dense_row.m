%% dense column detection
% This function is copied from the  SDPT3 checkdense.m

function  [idxden, idxsp] = detect_dense_row(A)
   [m,n] = size(A);
   idxden = []; 
   nzratio = 1;
   if (m > 1000); nzratio = 0.20; end
   if (m > 2000); nzratio = 0.10; end
   if (m > 5000); nzratio = 0.05; end
   if (nzratio < 1)
      nzrowA = reshape(full(sum(spones(A), 2)), 1, []); 
      idxden = find(nzrowA > nzratio*n); 
      % if (length(idxden) > max(200,0.1*n))
      %    idxden = [];
      % end 
   end

   idxsp = setdiff([1:m], idxden);
end
 