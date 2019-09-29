function [ eles ] = genUtmEleGetter(elevDataLats, elevDataLons, ...
    elevDataEles, xs, ys, utmZoneOrUtmToDegFct)
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

eles = interp2(elevDataLons, ...
    elevDataLats, elevDataEles, lons, lats);

end
% EOF