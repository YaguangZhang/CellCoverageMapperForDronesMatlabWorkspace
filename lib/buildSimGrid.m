function [ gridXYPts, gridResolution ] = buildSimGrid( boundaryXYs, ...
    parameter, flagByResolution )
%BUILDSIMGRID Create a grid
%
% Input:
%   - boundaryXYs
%     A N-by-2 matrix for N (x, y) vertices for the boundary of interest.
%     Note that the last row should be the same as the first row
%     (indicating a closed polygon).
%   - parameter, flagByResolution
%     If flagByResolution is present and set to true, parameter will be the
%     grid resolution. Otherwise, parameter will be the number of points
%     along the longer side of the area of interest.
%
% Outputs:
%   - gridXYPts
%     A matrix for the output grid, with each row being one grid point in
%     the form of (x, y).
%   - gridResolution
%     The spatial resolution of the output grid.
%
% Yaguang Zhang, Purdue, 03/16/2022

if ~exist('flagByResolution', 'var')
    flagByResolution = false;
end

mapMinX = min(boundaryXYs(:,1));
mapMaxX = max(boundaryXYs(:,1));
mapMinY = min(boundaryXYs(:,2));
mapMaxY = max(boundaryXYs(:,2));

if flagByResolution
    gridResolution = parameter;
else
    mapWidthInM = mapMaxX-mapMinX;
    mapHeightInM = mapMaxY-mapMinY;

    gridResolution = max([mapWidthInM, mapHeightInM])./parameter;
end

mapXLabels = constructAxisGrid( ...
    mean([mapMaxX, mapMinX]), ...
    floor((mapMaxX-mapMinX)./gridResolution), gridResolution);
mapYLabels = constructAxisGrid( ...
    mean([mapMaxY, mapMinY]), ...
    floor((mapMaxY-mapMinY)./gridResolution), gridResolution);
[mapXs,mapYs] = meshgrid(mapXLabels,mapYLabels);

% Discard map grid points out of the area of interest.
boolsMapGridPtsToKeep = InPolygon(mapXs(:), mapYs(:), ...
    boundaryXYs(:,1), ...
    boundaryXYs(:,2));

gridXYPts = [mapXs(boolsMapGridPtsToKeep), ...
    mapYs(boolsMapGridPtsToKeep)];
end
% EOF