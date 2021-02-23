function [deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(utmZone)
%GENUTMCONVERTERSFORFIXEDZONE A helper to generate functions to
%convert/compute GPS coordinates to/from UTM (x, y) for a fixed zone.
%
% Yaguang Zhang, Purdue, 02/02/2021

% Convert GPS degrees to UTM coordinates for the specified zone.
utmstruct_speZone = defaultm('utm');
% Remove white space in the zone label.
utmstruct_speZone.zone ...
    = utmZone(~isspace(utmZone));
utmstruct_speZone.geoid = wgs84Ellipsoid;
utmstruct_speZone = defaultm(utmstruct_speZone);

deg2utm_speZone = @(lat, lon) mfwdtran(utmstruct_speZone, lat,lon);
utm2deg_speZone = @(x, y) minvtran(utmstruct_speZone, x, y);

end
% EOF