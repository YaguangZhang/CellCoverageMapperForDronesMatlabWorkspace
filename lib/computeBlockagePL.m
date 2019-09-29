function [ blockagePL ] = computeBlockagePL(txXYAlt, rxXYAlt, ...
    lidarProfile, simConfigs)
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
%         The clearance ratio of the first Fresnel zone for a LoS path: at
%         least this ratio of the first Fresnel zone needs to be clear for
%         a path to be considered as "line of sight" (LoS); we expect the
%         value to be larger or equal to 50%.
%       - simConfigs.CARRIER_WAVELENGTH_IN_M
%         Carrier wavelength in meters.
% Output:
%   - blockagePL
%     The path loss for the blockage map. If the link is blocked, a NaN
%     will be returned.
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
lbsLidarProfileZs = lidarProfile(2:(end-1));

% Polynomial parameters for the LoS path.
parsLoSPath = polyfit([0; distTxToRx], ...
    [txXYAlt(3); rxXYAlt(3)], 1);
% LoS path height at the LiDAR profile sample locations.
losPathHs = polyval(parsLoSPath, obsLidarProfDists);

if all(losPathHs>=lbsLidarProfileZs)
    % A 50-percent clearance is now ensured. We need to consider the first
    % Fresnel zone here to more accurately determine the LoS paths.
    
    % Distance between the celluar tower and the lidar z locations.
    d1s = vecnorm( ...
        [obsLidarProfDists(:)'; ...
        lbsLidarProfileZs(:)'-txXYAlt(3)]);
    % Distance between the lidar z locations and the mobile device.
    d2s = vecnorm( ...
        [distTxToRx - obsLidarProfDists(:)'; ...
        lbsLidarProfileZs(:)'-rxXYAlt(3)]);
    % First Fresnel zone radii for the lidar z locations.
    firstFresRadii ...
        = (sqrt( ...
        (simConfigs.CARRIER_WAVELENGTH_IN_M .* d1s .* d2s)./(d1s + d2s)))';
    
    % The distances between the lidar z locations and the TX-RX direct
    % path.
    distsToDirectPath ...
        = point_to_line_distance( ...
        [obsLidarProfDists(:), lbsLidarProfileZs(:)], ...
        [0, txXYAlt(3)], ...
        [distTxToRx, rxXYAlt(3)]);
    
    % The extra clearance ratio respective to the first Fresnel zone radius
    % for LoS paths.
    extraClearRatioVsR ...
        = (simConfigs.LOS_FIRST_FRES_CLEAR_RATIO-0.5)*2;
    
    if all(distsToDirectPath ...
            -extraClearRatioVsR*firstFresRadii>0)
        blockagePL = fspl(distTxToRx, simConfigs.CARRIER_WAVELENGTH_IN_M);
    end
end

end
% EOF