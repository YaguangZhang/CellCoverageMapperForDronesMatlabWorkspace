function [blockagePL, coveragePL, blockageDistInM] ...
    = computeBlockageAndCoveragePLs( ...
    lidarProfile, terrainProfile, ...
    startXYH, endXYH, ...
    simConfigs)
%COMPUTEBLOCKAGEANDCOVERAGEPLS Compute the path loss (in dB) for the
%blockage and coverage maps.
%
% For the blockage map, we will
%   (1) Find line-of-sight (LoS) locations
%       We check the clearance of the 1st Fresnel zone and locations with
%       enough clearance will be considered as LoS.
%   (2) Evaluate path loss for LoS locations
%       The free space path loss (FSPL) model is used to evaluate the path
%       loss for non-blocked locations.
%
% For the coverage map, we will
%   (1) For drone locations with distances over (and including) 1km from
%   the cell tower
%       We use NTIA eHata model to evaluate the path loss.
%   (2) For drone locations with distances from the cell tower strictly
%   under 1km
%       If the path loss from the blockage map is valid for this drone
%       location, we will weight that value with the eHata value (please
%       see the code below for more details); otherwise, NaN will be
%       assigned for the path loss of that location.
%
% Yaguang Zhang, Purdue, 09/18/2019

% Treat the higher location as the TX location. 
%   Note: altitude = elevation + height.
startAlt = terrainProfile(1) + startXYH(3);
endPtAlt = terrainProfile(end) + endXYH(3);

% Inputs for the eHata model.
if startAlt>=endPtAlt
    txXYH = startXYH;
    rxXYH = endXYH;
else
    % Assuming reciprocal channels, we will flip TX and RX for cases where
    % the mobile antenna is higher than the cellular tower antenna, because
    % eHata is designed for higher TXs.
    txXYH = endXYH;
    rxXYH = startXYH;
    lidarProfile = lidarProfile(end:-1:1);
    terrainProfile = terrainProfile(end:-1:1);
end

%% Blockage
% Make sure there is no obstacles blocking the LoS path too much according
% to the LiDAR data. Note that here we always conider the profiles with the
% cellular tower as the start point and the mobile device as the end point.

txXYAlt = [txXYH(1:2), terrainProfile(1)+txXYH(3)];
rxXYAlt = [rxXYH(1:2), terrainProfile(end)+rxXYH(3)];

[blockagePL, blockageDistInM] ...
    = computeBlockagePLAndDist(txXYAlt, rxXYAlt, lidarProfile, simConfigs);

%% Coverage

% Generate the elevation profile needed by the extended Hata model
% according to our terrain information.
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
numTerrianPoints = length(terrainProfile);

elev = nan(numTerrianPoints+2, 1);
elev(1) = numTerrianPoints-1;
distTxToRx = norm(txXYH(1:2)-rxXYH(1:2));
elev(2) = distTxToRx/(numTerrianPoints-1);
elev(3:end) = terrainProfile;

coveragePL ...
    = computeCoveragePL(txXYH, rxXYH, elev, simConfigs, lidarProfile);

end
% EOF