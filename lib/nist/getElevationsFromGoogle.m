function [ alts ] ...
    = getElevationsFromGoogle( lats, lons, GOOGLE_MAPS_API)
%GETELEVATIONSFROMGOOGLE A wrapper to use getElevations with possibly nan
%lats and lons.
%
% Inputs:
%   - lats, lons
%     Column vectors representing the GPS locations.
%   - GOOGLE_MAPS_API
%     Optinal. The Google Maps Elevation key to use.
%
% Yaguang Zhang, Purdue, 04/27/2018

alts = nan(length(lats), 1);
boolsNotNan = ~(isnan(lats) | isnan(lons));

if exist('GOOGLE_MAPS_API', 'var')
    altsValid = getElevations(lats(boolsNotNan), lons(boolsNotNan), ...
        'key', GOOGLE_MAPS_API);
else
    altsValid = getElevations(lats(boolsNotNan), lons(boolsNotNan));
end

alts(boolsNotNan) = altsValid;

end
%EOF