function [ poly ] = coorRangeToPolyShape(coorRange)
%COORRANGETOPOLYSHAPE Convert the input coordinate range [minX, minY; maxX,
%maxY] to a rectangle polyshape.
% Note that the output vertices are arragned in a clockwise order. Yaguang
% Zhang, Purdue, 05/11/2021
minX = coorRange(1,1);
minY = coorRange(1,2);
maxX = coorRange(2,1);
maxY = coorRange(2,2);

vertices = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; minX, minY];
poly = polyshape(vertices(:,1), vertices(:,2));

end
% EOF