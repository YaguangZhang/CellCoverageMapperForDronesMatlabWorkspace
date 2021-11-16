function [ rect ] = ptPairToRect(ptPair)
%PTPAIRTORECT Convert a pair of points [minX, minY, maxX, maxY] to a
%rectangle polygon.
%
% Yaguang Zhang, Purdue, 11/15/2021

minX = ptPair(1);
minY = ptPair(2);
maxX = ptPair(3);
maxY = ptPair(4);

% Note that the output vertices are arragned in a clockwise order.
rect = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; minX, minY];

end
% EOF