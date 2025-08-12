function aimage=partialyt(image)
    [m,n]=size(image);
    aimage=zeros(m,n);
%     for j = 1:n
%         for i =2:m-1
%             aimage(i,j)=image(i-1,j)-image(i,j);
%         end
%         aimage(1,j)=-image(1,j);
%         aimage(n,j)=image(n-1,j);
%     end
    aimage(2:end-1,:) = image(1:end-2,:)-image(2:end-1,:);
        aimage(1,:)=-image(1,:);
        aimage(m,:)=image(m-1,:);    
end