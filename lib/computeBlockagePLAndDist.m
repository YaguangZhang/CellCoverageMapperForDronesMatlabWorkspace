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

% Set this to be true to generate debugging figures.
FLAG_DEBUG = false;

distTxToRx = norm(txXYAlt(1:2)-rxXYAlt(1:2));
numPoints = length(lidarProfile);
blockagePL = nan;

calcFirstFresRadii = @(d1s, d2s) sqrt( ...
    (simConfigs.CARRIER_WAVELENGTH_IN_M .* d1s .* d2s)./(d1s + d2s));

% We only need to compare the points between the cellular tower (TX) and
% the mobile device (RX). These points are considered as the top points of
% the obstacles.
lidarProfDists = linspace(0, distTxToRx, numPoints)';
% Ignore the LiDAR data at the TX and the RX.
obsLidarProfDists = lidarProfDists(2:(end-1));
obsLidarProfileZs = lidarProfile(2:(end-1));

% % Polynomial parameters for the LoS path.
%  parsLoSPath = polyfit([0; distTxToRx], ...
%     [txXYAlt(3); rxXYAlt(3)], 1);
% % LoS path height at the LiDAR profile sample locations.
%  losPathHs = polyval(parsLoSPath, obsLidarProfDists);
%
% For better speed (we assume distTxToRx is not zero):
losPathHs = obsLidarProfDists./distTxToRx ...
    .*(rxXYAlt(3)-txXYAlt(3))+txXYAlt(3);

boolsBlocked = nan(numPoints-2, 1);
% (i) Locations with LiDAR z points NOT below the direct path can be
% already treated as blocked.
boolsBelowDirectPath = obsLidarProfileZs<losPathHs;
boolsBlocked(~boolsBelowDirectPath) = true;

% For locations with LiDAR z below the direct path, we need to consider the
% first Fresnel zone to accurately determine the blockage distance.

% (ii) LiDAR z points below the direct path that are far enough from the
% direct path are "clear".
boolsTBD = isnan(boolsBlocked);

% Max 1st Fresnel zone radius.
distTxToRx3D = norm(txXYAlt - rxXYAlt);
firstFresRadiiMax = calcFirstFresRadii(distTxToRx3D/2, distTxToRx3D/2);

% The distances between the lidar z locations and the TX-RX direct path.
curDistsToDirectPath = point_to_line_distance( ...
    [obsLidarProfDists(boolsTBD), ...
    obsLidarProfileZs(boolsTBD)], ...
    [0, txXYAlt(3)], ...
    [distTxToRx, rxXYAlt(3)]);

% Note: simConfigs.LOS_FIRST_FRES_CLEAR_RATIO.*firstFresRadius is treated
% as the minimum distance required for reliable wireless communication
% links.
boolsBlockedTBD = nan(size(curDistsToDirectPath));
boolsBlockedTBD(curDistsToDirectPath >= ...
    simConfigs.LOS_FIRST_FRES_CLEAR_RATIO.*firstFresRadiiMax) = false;
boolsBlocked(boolsTBD) = boolsBlockedTBD;
if FLAG_DEBUG
    boolsBlockedType2 = boolsBlocked; %#ok<UNRCH>
end

% (iii) Remaining points will essentially be inspected one by one.
boolsTBD = isnan(boolsBlocked);
curDistsToDirectPath = curDistsToDirectPath(isnan(boolsBlockedTBD));
boolsBlockedTBD = nan(size(curDistsToDirectPath));

% Distances between the Tx and P, the intersection point of (a) the direct
% path with (b) the line which goes through that profile point and is
% perpendicular to the direct path.
d1s = sqrt( sum( ...
    ([obsLidarProfDists(boolsTBD), obsLidarProfileZs(boolsTBD)] ...
    - [0, txXYAlt(3)]).^2, 2) - curDistsToDirectPath.^2);
% % Distances between P and the Rx.
%  d2s = sqrt( sum(([obsLidarProfDists(boolsTBD), ...
%     obsLidarProfileZs(boolsTBD)] ...
%      - [distTxToRx, rxXYAlt(3)]).^2, 2) - curDistsToDirectPath.^2);
d2s = distTxToRx - d1s;

% (iii-1) P NOT on line segment Rx-Tx: clear.
boolsBlockedTBD(d2s<0) = false;
boolsBlocked(boolsTBD) = boolsBlockedTBD;
if FLAG_DEBUG
    boolsBlockedType3_1 = boolsBlocked; %#ok<UNRCH>
end

% (iii-2) All remaining points' blockage labels are determined by the
% point-by-point Fresnel Zone radius.
boolsTBD = isnan(boolsBlocked);
boolsNanBlockedTBD = isnan(boolsBlockedTBD);
curDistsToDirectPath = curDistsToDirectPath(boolsNanBlockedTBD);
d1s = d1s(boolsNanBlockedTBD);
d2s = d2s(boolsNanBlockedTBD);

% First Fresnel zone radii for the lidar z locations.
firstFresRadii = calcFirstFresRadii(d1s, d2s);

% The extra clearance ratio respective to the first Fresnel zone radius for
% LoS paths.
boolsBlocked(boolsTBD) = curDistsToDirectPath ...
    < simConfigs.LOS_FIRST_FRES_CLEAR_RATIO.*firstFresRadii;

% For simplicity, we will estimate blockageDistInM based on boolsBlocked
% (without considering the start and end points).
blockageDistInM = distTxToRx*(sum(boolsBlocked)/length(boolsBlocked));

if all(~boolsBlocked)
    % We will consider the 3D distance in FSPL computation.
    blockagePL = fspl(distTxToRx3D, simConfigs.CARRIER_WAVELENGTH_IN_M);
end

if FLAG_DEBUG
    dirToSaveDebugFig = fullfile(pwd, '..', 'PostProcessingResults', ...
        'blockageDistPathDebugFigs'); %#ok<UNRCH>
    if ~exist(dirToSaveDebugFig, 'dir')
        mkdir(dirToSaveDebugFig)
    end

    absPathToSaveDebugFig = fullfile(dirToSaveDebugFig, ...
        ['path_timestampInMs_now_', ...
        num2str(floor(now*24*60*60*1000), '%d'), '_rand_', ...
        num2str(floor(rand*10^6), '%d')]);

    boolsBlocked = logical(boolsBlocked);
    boolsBlockedType3_1 = boolsBlockedType3_1 == true;
    boolsBlockedType2 = boolsBlockedType2 == true;

    % A side view of the path with profile points.
    hFigPath = figure('Position', [0, 0, 1000, 500]); hold on;
    hTx = plot(lidarProfDists(1), txXYAlt(3), 'vg');
    hRx = plot(lidarProfDists(end), rxXYAlt(3), 'ob');
    hLoS = plot([lidarProfDists(1), lidarProfDists(end)], ...
        [txXYAlt(3), rxXYAlt(3)], '--k');
    hProf = plot(lidarProfDists, lidarProfile, '.k', 'MarkerSize', 9);
    hBlocked3_2 = plot(obsLidarProfDists(boolsBlocked), ...
        obsLidarProfileZs(boolsBlocked), '.r');
    hBlocked3_1 = plot(obsLidarProfDists(boolsBlockedType3_1), ...
        obsLidarProfileZs(boolsBlockedType3_1), '.y');
    hBlocked2 = plot(obsLidarProfDists(boolsBlockedType2), ...
        obsLidarProfileZs(boolsBlockedType2), '.c');
    hBlocked1 = plot(obsLidarProfDists(~boolsBelowDirectPath), ...
        obsLidarProfileZs(~boolsBelowDirectPath), '.m');
    % axis equal;
    legend([hTx, hRx, hLoS, hProf, ...
        hBlocked1, hBlocked2, hBlocked3_1, hBlocked3_2], ...
        'Tx', 'Rx', 'LoS Path', 'Profile', ...
        'Blockage (i)', 'Blockage (ii)', ...
        'Blockage (iii-1)', 'Blockage (iii-2)');
    transparentizeCurLegends;

    saveas(hFigPath, [absPathToSaveDebugFig, '.jpg']);
    close(hFigPath);
end
end
% EOF