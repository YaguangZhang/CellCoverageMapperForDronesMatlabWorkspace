function [newLonLatAxis, weightForWidth] ...
    = extendLonLatAxisByFactor(oldLonLatAxis, factorExtend, simConfigs)
%EXTENDLONLATAXISBYFACTOR Extend a (lon, lat) axis by the input factor in
%UTM, and output the resultant axis in (lon, lat).
%
% Inputs:
%   - oldLonLatAxis
%     The axis (lonMin, lonMat, latMin, latMax) for the current (lon, lat)
%     plot that is to be extended.
%   - factorExtend
%     A float to indicate the factor for extension.
%   - simConfigs
%     A structure for the configuration parameters for the simulation. We
%     will need fields:
%       - deg2utm_speZone, utm2deg_speZone
%         The functions to convert GPS (lat, lon) to (and back from) UTM
%         (x, y) for the LiDAR data area.
%
% Outputs:
%   - newLonLatAxis
%     The output extended axis.
%   - weightForWidth
%     A parameter to adjust the (lon, lat) plot side weight for precise
%     figure size adjustment. This is to be used for
%     adjustFigSizeByContent.m.
%
% Yaguang Zhang, Purdue, 10/18/2019

[oldAxisMinX, oldAxisMinY] = simConfigs.deg2utm_speZone( ...
    oldLonLatAxis(3), oldLonLatAxis(1));
[oldAxisMaxX, oldAxisMaxY] = simConfigs.deg2utm_speZone( ...
    oldLonLatAxis(4), oldLonLatAxis(2));

axisXYToSet = extendAxisByFactor( ...
    [oldAxisMinX, oldAxisMaxX, oldAxisMinY, oldAxisMaxY], factorExtend);

[axisToSetMinLat, axisToSetMinLon] = simConfigs.utm2deg_speZone( ...
    axisXYToSet(1), axisXYToSet(3));
[axisToSetMaxLat, axisToSetMaxLon] = simConfigs.utm2deg_speZone( ...
    axisXYToSet(2), axisXYToSet(4));

newLonLatAxis = [axisToSetMinLon, axisToSetMaxLon, ...
    axisToSetMinLat, axisToSetMaxLat];

weightForWidth = (axisXYToSet(2)-axisXYToSet(1)) ...
    ./(newLonLatAxis(2)-newLonLatAxis(1)) ...
    ./((axisXYToSet(4)-axisXYToSet(3))...
    ./(newLonLatAxis(4)-newLonLatAxis(3)));
end
% EOF