function [ utmPolyMat, utmZone ] ...
    = constructUtmRectanglePolyMat( rangeLatLonPts )
%CONSTRUCTUTMRECTANGLEPOLYMAT Convert rangeLatLonPts = [minLat, minLon;
%maxLat, maxLon] to a polygon matrix representing a rectangle defined by
%the input range.
%
% Yaguang Zhang, Purdue, 09/19/2019

[minX, minY, utmZone] ...
    = deg2utm(rangeLatLonPts(1,1), rangeLatLonPts(1,2));
[maxX, maxY, secUtmZone] ...
    = deg2utm(rangeLatLonPts(2,1), rangeLatLonPts(2,2));

assert(strcmp(utmZone, secUtmZone), ...
    'GPS range points are not in the same UTM zone!');

% Note that the output vertices are arragned in a clockwise order.
utmPolyMat = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; minX, minY];

end
% EOF