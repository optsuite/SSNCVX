function aimage=partialxt(image)
    [m,n]=size(image);
    aimage=zeros(m,n);
%     for i = 1:m
%         for j =2:n-1
%             aimage(i,j)=image(i,j-1)-image(i,j);
%         end
%         aimage(i,1)=-image(i,1);
%         aimage(i,n)=image(i,n-1);
%     end
    aimage(:,2:end-1) = image(:,1:end-2)-image(:,2:end-1);
        aimage(:,1)=-image(:,1);
        aimage(:,n)=image(:,n-1);    
end