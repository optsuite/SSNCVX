function aimage = DToperater(image)
[m,~] = size(image);
m = m/2;
 a1 = partialxt(image(1:m,:));
 a2 = partialyt(image(m+1:end,:));
 aimage = a1 +a2;
end