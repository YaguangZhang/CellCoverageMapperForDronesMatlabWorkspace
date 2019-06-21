function [ lidarZs ] = genRasterLidarZGetter( ...
    getLiDarZFromStatePlaneXYFct, fctLonLatToLidarStatePlaneXY, ...
    xs, ys, utmZone)
%GENRASTERLIDARZGETTER A helper function to generate the function
%   [zs] = getLiDarZFromXYFct(xs, ys).
%
% Inputs:
%   - getLiDarZFromStatePlaneXYFct
%     A function to get raster LiDAR z data from state plane coordinates
%     (spXs, spYs):
%           lidarZs = getLiDarZFromStatePlaneXYFct(spXs, spYs).
%   - fctLonLatToLidarStatePlaneXY
%     A function to convert (lons, lats) to state plane coordinates (spXs,
%     spYs):
%           [spXs, spYs] = fctLonLatToLidarStatePlaneXY(lons, lats).
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
[spXs, spYs] = fctLonLatToLidarStatePlaneXY(lons, lats);
lidarZs = getLiDarZFromStatePlaneXYFct(spXs, spYs);

end
% EOF