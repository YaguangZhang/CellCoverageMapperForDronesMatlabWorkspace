function [hCurBlockageMap, hCurHandleTxs] ...
    = plotBlockageMap(matRxLonLatWithPathLoss, cellAntLonLats, ...
    simConfigs, flagVisible, flagZoomIn, customFigSize)
%PLOTBLOCKAGEMAP Plot the blockage map.
%
% Generate a figure on Google map to show the blockage areas. We will use
% triangulation-based nearest neighbor interpolation for the (lon, lat)
% grid to get the blockage values.
%
% Inputs:
%   - matRxLonLatWithPathLoss
%     A matrix with each row being the (lon, lat, path loss) for a pixel of
%     the path loss map to be shown. The blocked areas should have a path
%     loss value of NaN.
%   - cellAntLoc
%     The (lon, lat) coordinates for the cellular antenna.
%   - simConfigs
%     The configuration struct for the simulation. We need fields:
%       - UTM_X_Y_BOUNDARY_OF_INTEREST
%         The UTM (x, y) polygon boundary vertices representing the area of
%         interest for generating the coverage maps; for presets
%         ExtendedTipp and IN, this is default to the range covered by the
%         availabe LiDAR data set for the corresponding area of interest
%         when it is empty.
%   - flagVisible
%     An optional boolean to control whether to show the plot or not.
%   - flagZoomIn
%     An optional boolean to control whether to zoom in to fit the area of
%     the input path loss map or not.
%   - customFigSize
%     An optional boolean to specify the desired figure size. The resultant
%     figure will start with this size and be adjusted according to the
%     figure content.
%
% Outputs:
%   - hCurBlockageMap
%     The handle to the resultant figure.
%   - hCurHandleTxs
%     The handle to the cellular tower.
%
% Yaguang Zhang, Purdue, 10/02/2019

legendBackgroundColor = ones(1,3).*0.8;

% For plotting.
colorTowers = 'w';
markerTowers = 'x';
markerSizeTowers = 6;
lineWidthTowers = 1;

% We will use the first color for clearance and the last color for
% blockage.
COLORMAP_TO_USE = 'jet';
% The location of the legend.
LEGEND_LOC = 'NorthEast';

if ~exist('customFigSize', 'var')
    customFigSize = [500, 500].*0.7;
end

% By default, show the plot.
if ~exist('flagVisible', 'var')
    flagVisible = true;
end

% By default, do not zoom in to the path loss map, so that all TXs will be
% shown.
if ~exist('flagZoomIn', 'var')
    flagZoomIn = false;
end

colorRange = [0,1];

hCurBlockageMap = figure('visible', flagVisible, ...
    'Position', [0,0,customFigSize]);
hCurBlockageMap.InvertHardcopy = 'off';
hCurBlockageMap.Color = 'none';
hold on;

% Area of interest.
[areaOfInterestLats, areaOfInterestLons] = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
plot(polyshape(areaOfInterestLons, areaOfInterestLats), ...
    'FaceColor', 'white');

if flagZoomIn
    extensionFactor = 0.05;
else
    % Extend the content by a constant factor in the UTM system.
    extensionFactor = 0.2;
    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'ExtendedTipp')
        extensionFactor = 0.25;
    end
end

[axisLonLatToSet, weightForWidth] ...
    = extendLonLatAxisByFactor( ...
    [min(areaOfInterestLons), max(areaOfInterestLons), ...
    min(areaOfInterestLats), max(areaOfInterestLats)], ...
    extensionFactor, simConfigs);

% TX.
hCurHandleTxs = plot3(cellAntLonLats(:,1), cellAntLonLats(:,2), ...
    ones(length(cellAntLonLats(:,1)), 1), ...
    markerTowers, ...
    'MarkerSize', markerSizeTowers, ...
    'Color', colorTowers, ...
    'LineWidth', lineWidthTowers);

% Plot simulation results.
xLabelToSet = 'Longitude';
yLabelToSet = 'Latitude';

set(gca, 'fontWeight', 'bold');

% Create meshgrid for surf. We will increase the grid density by a factor
% of 10 to better show the blockage areas.
upSampFactor = 10;
sufNumPtsPerSide = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE.*upSampFactor;
lons = matRxLonLatWithPathLoss(:, 1);
lats = matRxLonLatWithPathLoss(:, 2);
zs = matRxLonLatWithPathLoss(:, 3);

% Set blockage area to 1 and other areas to 0.
boolsBlocked = isnan(zs);
zs(boolsBlocked) = 1;
zs(~boolsBlocked) = 0;

% Find the ranges for the boundary of interet (BoI) to build a new grid for
% showing the results.
[latsBoI, lonsBoI] = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
lonMinBoI = min(lonsBoI);
lonMaxBoI = max(lonsBoI);
latMinBoI = min(latsBoI);
latMaxBoI = max(latsBoI);

[lonsNew, latsNew] = meshgrid( ...
    linspace(lonMinBoI, lonMaxBoI, sufNumPtsPerSide), ...
    linspace(latMinBoI, latMaxBoI, sufNumPtsPerSide));
zsNew = griddata(lons,lats,zs,lonsNew,latsNew,'Nearest');

% Ignore points out of the area of interest by seting the z values for them
% to NaN.
[in,on] = inpolygon(lonsNew(:), latsNew(:), lonsBoI, latsBoI);
boolsPtsToIgnore = ~(in|on);
if any(boolsPtsToIgnore)
    zsNew(boolsPtsToIgnore) = nan;
end

% Plot the blockage areas.
hRxs = surf(lonsNew,latsNew,zsNew, ...
    'FaceAlpha',0.5, 'EdgeColor', 'none');
curColormap = colormap(COLORMAP_TO_USE);
hClear = plot(polyshape(nan(3,2)), 'FaceColor', curColormap(1, :));
hBlocked = plot(polyshape(nan(3,2)), 'FaceColor', curColormap(end, :));
caxis(colorRange);
xticks([]); yticks([]);
xlabel(xLabelToSet); ylabel(yLabelToSet);

hLeg = legend([hCurHandleTxs, hClear, hBlocked], ...
    'Cell towers', 'Clear', 'Blocked', ...
    'Location', LEGEND_LOC); view(2);
set(hLeg, 'color', legendBackgroundColor);
set(hCurBlockageMap, 'Color', 'w');

adjustFigSizeByContent(hCurBlockageMap, ...
    axisLonLatToSet, 'height', weightForWidth*0.9);

plot_google_map('MapType', 'satellite');

end
% EOF