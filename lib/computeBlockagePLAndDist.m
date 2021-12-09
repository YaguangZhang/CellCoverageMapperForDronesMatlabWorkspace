function [ blockagePL, blockageDistInM ] ...
    = computeBlockagePLAndDist(txXYAlt, rxXYAlt, lidarProfile, simConfigs)
%COMPUTEBLOCKAGEPL Compute the path loss (in dB) for a pixel of the
%blockage map.
%
% For the blockage map, we will
%   (1) Find line-of-sight (LoS) locations
%       We check the clearance of the 1st Fresnel zone and locations with
%       enough clearance will be considered as LoS.
%   (2) Evaluate path loss for LoS locations
%       The free space path loss (FSPL) model is used to evaluate the path
%       loss for non-blocked locations.
%
% Inputs:
%   - txXYAlt, rxXYAlt
%     The UTM (x,y) with altitude (elevation + height) for the TX and the
%     RX, respectively.
%   - lidarProfile
%     The LiDAR pofile (LiDAR z in meters) between (and including) the TX
%     and the RX.
%   - simConfigs
%     A structure holding basic parameters of the simulation; we need here:
%       - simConfigs.LOS_FIRST_FRES_CLEAR_RATIO
%         The clearance ratio of the first Fresnel zone around the direct
%         path to ensure stable communication links: at least this ratio of
%         the first Fresnel zone needs to be clear for a path to be
%         considered as "line of sight" (LoS); we expect the value to be at
%         least 60% (80% is recommended; for solid system designs, it is
%         not rare to have the full 1st Fresnel zone clear).
%       - simConfigs.CARRIER_WAVELENGTH_IN_M
%         Carrier wavelength in meters.
% Output:
%   - blockagePL
%     The path loss for the blockage map. If the link is blocked, a NaN
%     will be returned.
%   - blockageDistInM
%     The distance (in meters) blocked according to the LiDAR profile.
%
% Yaguang Zhang, Purdue, 09/18/2019

distTxToRx = norm(txXYAlt(1:2)-rxXYAlt(1:2));
numPoints = length(lidarProfile);
blockagePL = nan;

% We only need to compare the points between the cellular tower (TX) and
% the mobile device (RX). These points are considered as the top points of
% the obstacles.
obsLidarProfDists = linspace(0, distTxToRx, numPoints)';
% Ignore the LiDAR data at the TX and the RX.
obsLidarProfDists = obsLidarProfDists(2:(end-1));
obsLidarProfileZs = lidarProfile(2:(end-1));

% Polynomial parameters for the LoS path.
parsLoSPath = polyfit([0; distTxToRx], ...
    [txXYAlt(3); rxXYAlt(3)], 1);
% LoS path height at the LiDAR profile sample locations.
losPathHs = polyval(parsLoSPath, obsLidarProfDists);

% Locations with LiDAR z points NOT below the direct path can be already
% treated as blocked.
boolsBelowDirectPath = obsLidarProfileZs<losPathHs;
boolsBlocked = ~boolsBelowDirectPath;

% For locations with LiDAR z below the direct path, we need to consider the
% first Fresnel zone to accurately determine the blockage distance.

% Distances between the celluar tower and the lidar z locations.
d1s = vecnorm( ...
    [obsLidarProfDists(boolsBelowDirectPath)'; ...
    obsLidarProfileZs(boolsBelowDirectPath)'-txXYAlt(3)]);
% Distances between the lidar z locations and the mobile device.
d2s = vecnorm( ...
    [distTxToRx - obsLidarProfDists(boolsBelowDirectPath)'; ...
    obsLidarProfileZs(boolsBelowDirectPath)'-rxXYAlt(3)]);
% First Fresnel zone radii for the lidar z locations.
firstFresRadii ...
    = (sqrt( ...
    (simConfigs.CARRIER_WAVELENGTH_IN_M .* d1s .* d2s)./(d1s + d2s)))';

% The distances between the lidar z locations and the TX-RX direct path.
distsToDirectPath ...
    = point_to_line_distance( ...
    [obsLidarProfDists(boolsBelowDirectPath), ...
    obsLidarProfileZs(boolsBelowDirectPath)], ...
    [0, txXYAlt(3)], ...
    [distTxToRx, rxXYAlt(3)]);

% The extra clearance ratio respective to the first Fresnel zone radius for
% LoS paths.
boolsBlocked(boolsBelowDirectPath) = distsToDirectPath ...
    <= simConfigs.LOS_FIRST_FRES_CLEAR_RATIO.*firstFresRadii;

% For simplicity, we will estimate blockageDistInM based on boolsBlocked
% (without considering the start and end points).
blockageDistInM = distTxToRx*(sum(boolsBlocked)/length(boolsBlocked));

if all(~boolsBlocked)
    % We will consider the 3D distance in FSPL computation.
    distTxToRx3D = norm(txXYAlt - rxXYAlt);
    blockagePL = fspl(distTxToRx3D, simConfigs.CARRIER_WAVELENGTH_IN_M);
end
end
% EOF