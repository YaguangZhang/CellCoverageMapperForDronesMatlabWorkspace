function [ utmPolyMat, utmZone ] ...
    = constructUtmRectanglePolyMat( rangeLatLonPts, deg2utm_speZone )
%CONSTRUCTUTMRECTANGLEPOLYMAT Convert rangeLatLonPts = [minLat, minLon;
%maxLat, maxLon] to a polygon matrix representing a rectangle defined by
%the input range.
%
% Update 20220901: Add optional input.
%   - deg2utm_speZone
%     Fct to convert GPS (lat, lon) to UTM (x, y) with zone information
%     imbeded. If present, the conversion will be done using this function
%     instead and the "same UTM zone" check will be skipped.
%
% Yaguang Zhang, Purdue, 09/19/2019

if exist('deg2utm_speZone', 'var')
    [minX, minY] = deg2utm_speZone( ...
        rangeLatLonPts(1,1), rangeLatLonPts(1,2));
    [maxX, maxY] = deg2utm_speZone( ...
        rangeLatLonPts(2,1), rangeLatLonPts(2,2));
else
    [minX, minY, utmZone] ...
        = deg2utm(rangeLatLonPts(1,1), rangeLatLonPts(1,2));
    [maxX, maxY, secUtmZone] ...
        = deg2utm(rangeLatLonPts(2,1), rangeLatLonPts(2,2));

    assert(strcmp(utmZone, secUtmZone), ...
        'GPS range points are not in the same UTM zone!');
end

% Note that the output vertices are arragned in a clockwise order.
utmPolyMat = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; minX, minY];

end
% EOF