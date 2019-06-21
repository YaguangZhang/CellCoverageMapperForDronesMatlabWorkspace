function [ eles ] = genUtmEleGetter(elevDataLats, elevDataLons, ...
    elevDataEles, xs, ys, utmZone)
%GENUTMELEGETTER A helper function to generate the function
%   [eles] = getEleFromXYFct(xs, ys).
%
% Inputs:
%   - elevDataLats, elevDataLons, elevDataEles
%     Three matrices for the longitude, latitude and elevation of the area
%     of interest. The elevation for the location
%           (elevDataLats(m,n), elevDataLons(m,n))
%     is elevDataLons(m,n). Also, elevDataLats and elevDataLons should be
%     generated via meshgrid from grid vectors that are strictly
%     increasing.
%   - xs, ys, utmZone
%     The UTM coordinates (xs,ys), in the form of column vectors, and their
%     zone string for the target locations where we will estimate the
%     output elevations according to the know elevation data.
%
% Output:
%   - eles
%     A column vector for the elevations for locations (xs, ys).
%
% Note that all the locations should be in the same UTM zone.
%
% Yaguang Zhang, Purdue, 06/13/2019

[lats, lons] = utm2deg(xs, ys, repmat(utmZone, length(xs), 1));

eles = interp2(elevDataLons, ...
    elevDataLats, elevDataEles, lons, lats);
end
% EOF