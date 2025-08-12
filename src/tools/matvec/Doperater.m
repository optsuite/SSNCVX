function aimage = Doperater(image)
    a1 = partialx(image);
    a2 = partialy(image);
    aimage = [a1; a2];
end
