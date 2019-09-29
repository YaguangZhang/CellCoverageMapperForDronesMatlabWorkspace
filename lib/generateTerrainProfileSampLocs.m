function [ sampLocs ] ...
    = generateTerrainProfileSampLocs(startLocXY, endLocXY, ...
    maxAllowedResInM, minNumSamps)
%GENERATETERRAINPROFILESAMPLOCS Get a list of locations along the TX and RX
%for the terrain profile.
%
% Inputs:
%   - startLocXY, endLocXY
%     The UTM (x,y) coordinates for the start point and end point of the
%     terrain profile, repectively.
%   - maxAllowedResInM
%     The maximum distance between adjacent samples in meter.
%   - minNumSamps
%     The minimum number of resulting sample points.
%
% Outpus:
%   - sampLocs
%     A matrix with each row being a sample location for the terrain
%     profile between the input start point and end point; these locations
%     are ordered from the start point (the first row) to the end point
%     (the last row).
%
% Yaguang Zhang, Purdue, 09/17/2019

% Reference - the elevation profile needed by the extended Hata model
%   - elev
%     An array containing elevation profile between Tx & Rx, where:
%       - elev(1) = numPoints - 1
%         (for both Matlab eHata lib and C++ eHata lib; note, numPoints is
%         the number of points between Tx & Rx)
%       - elev(2) = distance between points (in meters).
%         (thus, elev(1)*elev(2)=distance between Tx & Rx)
%       - elev(3) = Tx elevation
%         (in meters)
%       - elev(numPoints+2) = Rx elevation
%         (in meters)

distStartToEnd = norm(startLocXY-endLocXY);
numSampLocs = max( ...
    ceil(distStartToEnd./maxAllowedResInM), ...
    minNumSamps);

eleProfXs = linspace(startLocXY(1), endLocXY(1), numSampLocs)';
eleProfYs = linspace(startLocXY(2), endLocXY(2), numSampLocs)';

sampLocs = [eleProfXs, eleProfYs];

end
% EOF