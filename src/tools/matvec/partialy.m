function aimage=partialy(image)
    [m, n]=size(image);
    aimage=zeros(m,n);
%     for i =1:m-1
%         for j=1:n
%             aimage(i,j)=image(i+1,j)-image(i,j);
%         end
%     end
    aimage(1:end-1,:)=image(2:end,:)-image(1:end-1,:);
end