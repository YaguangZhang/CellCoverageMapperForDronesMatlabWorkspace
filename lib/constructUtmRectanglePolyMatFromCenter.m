function [utmPolyMat, utmZone] ...
    = constructUtmRectanglePolyMatFromCenter( ...
    centerLatLon, sideLengthsInM)
%CONSTRUCTUTMRECTANGLEPOLYMATFROMCENTER Convert the center of the rectangle
%and the side lengths sideLengthsInM = [lengthAlongX, lengthAlongY] to a
%polygon matrix representing a rectangle defined by the input range.
%
% Yaguang Zhang, Purdue, 11/05/2019

halfSideLenX = sideLengthsInM(1)/2;
halfSideLenY = sideLengthsInM(2)/2;
[cenX, cenY, utmZone] = deg2utm(centerLatLon(1), centerLatLon(2));
minX = cenX-halfSideLenX;
maxX = cenX+halfSideLenX;
minY = cenY-halfSideLenY;
maxY = cenY+halfSideLenY;

% Note that the output vertices are arragned in a clockwise order.
utmPolyMat = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; minX, minY];

end
% EOF