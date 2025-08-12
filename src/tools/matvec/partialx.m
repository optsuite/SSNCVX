function aimage=partialx(image)
    [m, n]=size(image);
    aimage=zeros(m,n);
%     for i =1:m
%         for j=1:n-1
%             aimage(i,j)=image(i,j+1)-image(i,j);
     aimage(:,1:end-1) = image(:,2:end) - image(:,1:end-1);
%         end
%     end
end