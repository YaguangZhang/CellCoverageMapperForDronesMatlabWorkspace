function [hFigRegRmsdWrtMean, gridXYPts, regRmsdsWrtMean] ...
    = plotRegionalRmsdWrtMean(lidarXYZ, regionalCircleRadiusInM, ...
    numberOfPixelsOnLongerSide, figVisible)
%PLOTREGIONALRMSDWRTMEAN Estimate the ground roughness via regional root
%mean square different (RMSD) on a grid constructed for the area covered by
%the LiDAR data.
%
% Inputs:
%   - lidarXYZ
%     A matrix with each row being (UTM x, UTM y, LiDAR z) for a sample of
%     a LiDAR grid dataset.
%   - regionalCircleRadiusInM
%     The radius in meter for the regional circular area where the LiDAR
%     data will be considered to compute the regional RMSD. For a location
%     of interest, we will fetch the LiDAR samples in this circular region,
%     compute the mean for LiDAR z, and use that as the reference value for
%     further get the RMSD.
%   - numberOfPixelsOnLongerSide
%     The number of grid points used in along the long side of the area of
%     interest.
%   - figVisible
%     A illustration plot for the resulting ground roughness will also be
%     generated. This flag constrols whether this plot will be visible.
%
% Outpus:
%   - hFigRegRmsdWrtMean
%     The handler to the ground roughness illustration figure.
%   - gridXYPts
%     A matrix with each row being the grid point coordinates (UTM x, UTM
%     y) where the ground roughness is estimated.
%   - regRmsdsWrtMean
%     A column vector with each element being the regional RMSD for the
%     corresponding grid piont.
%
% Yaguang Zhang, Purdue, 11/21/2019

%% Build grid for the resulting plot.

indicesLidarAreaConvHull = convhull(lidarXYZ(:,1), lidarXYZ(:,2));
UTM_X_Y_BOUNDARY_OF_INTEREST = [lidarXYZ(indicesLidarAreaConvHull,1), ...
    lidarXYZ(indicesLidarAreaConvHull,2)];

mapMinX = min(UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
mapMaxX = max(UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
mapMinY = min(UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
mapMaxY = max(UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

mapWidthInM = mapMaxX-mapMinX;
mapHeightInM = mapMaxY-mapMinY;

gridResolution = max([mapWidthInM, mapHeightInM]) ...
    ./numberOfPixelsOnLongerSide;

mapXLabels = constructAxisGrid( ...
    mean([mapMaxX, mapMinX]), ...
    floor((mapMaxX-mapMinX)./gridResolution), gridResolution);
mapYLabels = constructAxisGrid( ...
    mean([mapMaxY, mapMinY]), ...
    floor((mapMaxY-mapMinY)./gridResolution), gridResolution);
[mapXs,mapYs] = meshgrid(mapXLabels,mapYLabels);

% Discard map grid points out of the area of interest.
boolsMapGridPtsToKeep = inpolygon(mapXs(:), mapYs(:), ...
    UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

% The grid point locations to evaluate the regional RMSD.
gridXYPts = [mapXs(boolsMapGridPtsToKeep), ...
    mapYs(boolsMapGridPtsToKeep)];
numGridPts = sum(boolsMapGridPtsToKeep);

%% Compute regional RMSD.

regRmsdsWrtMean = nan(numGridPts, 1);
for idxGridPt = 1:numGridPts
    curGridPtXY = gridXYPts(idxGridPt, :);
    
    % First filter out Lidar points out of the maximum x and y ranges.
    minXToConsider = curGridPtXY(1)-regionalCircleRadiusInM;
    maxXToConsider = curGridPtXY(1)+regionalCircleRadiusInM;
    minYToConsider = curGridPtXY(2)-regionalCircleRadiusInM;
    maxYToConsider = curGridPtXY(2)+regionalCircleRadiusInM;
    
    boolsLidarPtsToConsider = (lidarXYZ(:,1)>=minXToConsider) ...
        & (lidarXYZ(:,1)<=maxXToConsider) ...
        & (lidarXYZ(:,2)>=minYToConsider) ...
        & (lidarXYZ(:,2)<=maxYToConsider);
    
    lidarXYZToConsider = lidarXYZ(boolsLidarPtsToConsider, :);
    
    % Further filter out Lidar points too far away.
    distsFromLidarPtsToCurGridPt = vecnorm( ...
        [lidarXYZToConsider(:,1)-curGridPtXY(1), ...
        lidarXYZToConsider(:,2)-curGridPtXY(2)]')';
    boolsLidarPtsToConsider = ...
        distsFromLidarPtsToCurGridPt<=regionalCircleRadiusInM;
    lidarXYZToConsider = lidarXYZToConsider(boolsLidarPtsToConsider, :);
    
    % Regional mean.
    regMean = mean(lidarXYZToConsider(:,3));
    regRmsdsWrtMean(idxGridPt) ...
        = mean((lidarXYZToConsider(:,3) - regMean).^2);
end

%% Plot.

hFigRegRmsdWrtMean = surfWithContour([gridXYPts, regRmsdsWrtMean], ...
    {'', '', '', '', 'Regional RMSD (m)', ''}, figVisible);

end
% EOF