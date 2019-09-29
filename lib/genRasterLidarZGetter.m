function [ lidarZs ] = genRasterLidarZGetter( ...
    getLiDarZFromStatePlaneXYFct, fctLonLatToLidarStatePlaneXY, ...
    xs, ys, utmZoneOrUtmToDegFct)
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
%   - xs, ys
%     The UTM coordinates (xs,ys), in the form of column vectors for the
%     target locations where we will estimate the output elevations
%     according to the know elevation data.
%   - utmZoneOrUtmToDegFct
%     Either the zone string for (xs, ys), or a function to convert UTM
%     (xs, ys) to (lat, lon).
%
% Output:
%   - eles
%     A column vector for the elevations for locations (xs, ys).
%
% Note that all the locations should be in the same UTM zone.
%
% Yaguang Zhang, Purdue, 06/13/2019

switch class(utmZoneOrUtmToDegFct)
    case 'function_handle'
        [lats, lons] = utmZoneOrUtmToDegFct(xs, ys);
    case 'char'
        [lats, lons] = utm2deg(xs, ys, ...
            repmat(utmZoneOrUtmToDegFct, length(xs), 1));
    otherwise
        error('Unsupported type of input utmZoneOrUtmToDegFct!');
end

[spXs, spYs] = fctLonLatToLidarStatePlaneXY(lons, lats);
lidarZs = getLiDarZFromStatePlaneXYFct(spXs, spYs);

end
% EOF