function [hCurPLMap, hCurHandleTxs] ...
    = plotPathLossMap(matRxLonLatWithPathLoss, cellAntLonLats, ...
    simConfigs, flagVisible, flagZoomIn, flagCmdToPlotPLs, customFigSize)
%PLOTPATHLOSSMAP Plot the path loss map.
%
% Generate a figure on Google map to show the path loss map.
%
% Inputs:
%   - matRxLonLatWithPathLoss
%     A matrix with each row being the (lon, lat, path loss) for a pixel of
%     the path loss map to be shown.
%   - cellAntLoc
%     The (lon, lat) coordinates for the cellular antenna.
%   - simConfigs
%     The configuration struct for the simulation. We need fields:
%       - ALLOWED_PATH_LOSS_RANGE_IN_DB
%         Color mapping limits for plot3k. A two element vector specifying
%         the values that map to the first and last colors. This is useful
%         for generating a series of plots with identical coloring. The
%         colormap but not the colorbar) is flipped upside down if
%         'ColorRange' is given as [max min] instead of [min max].
%       - UTM_X_Y_BOUNDARY_OF_INTEREST
%         The UTM (x, y) polygon boundary vertices representing the area of
%         interest for generating the coverage maps.
%       - utm2deg_speZone
%         Function to convert GPS degrees (lat, lon) from UTM (x, y).
%       - CURRENT_SIMULATION_TAG
%         A string label to identify this simulation.
%       - NUM_OF_PIXELS_FOR_LONGER_SIDE
%         The number of pixels for the longer side (width/height) of the
%         map; the number of pixels for the other side will be proportional
%         to its length.
%   - flagVisible
%     An optional boolean to control whether to show the plot or not.
%   - flagZoomIn
%     An optional boolean to control whether to zoom in to fit the area of
%     the input path loss map or not.
%   - flagCmdToPlotPLs
%     An optional string to control what command to use in plotting the
%     path loss map. Currently support: 'plot3k' and 'surf'.
%
%     If 'surf' is chosen, we will need a new meshgrid for the plot. In
%     this case, to obtain the z values for points in the new grid, we will
%     use triangulation-based linear interpolation (which may cause NaN
%     values for the edge points) for all grid points first, and use
%     triangulation-based nearest neighbor interpolation for the grid edge
%     points with NaN values.
%   - customFigSize
%     An optional boolean to specify the desired figure size. The resultant
%     figure will start with this size and be adjusted according to the
%     figure content.
%
% Outputs:
%   - hCurPLMap
%     The handle to the resultant figure.
%   - hCurHandleTxs
%     The handle to the cellular tower.
%
% Yaguang Zhang, Purdue, 10/02/2019

colorMap = 'jet';
legendBackgroundColor = ones(1,3).*0.8;
minPathLossInDbExpected = 80;
flagRiseTxToTop = true;

% For plotting.
colorTowers = 'w';
markerTowers = 'x';
markerSizeTowers = 6;
lineWidthTowers = 1;

% The location of the legend.
LEGEND_LOC = 'NorthEast';

% By default, show the plot.
if ~exist('flagVisible', 'var')
    flagVisible = true;
end

% By default, do not zoom in to the path loss map, so that all TXs will be
% shown.
if ~exist('flagZoomIn', 'var')
    flagZoomIn = false;
end

% We support: 'plot3k' and 'surf'(default).
if ~exist('flagCmdToPlotPLs', 'var')
    flagCmdToPlotPLs = 'surf';
end

if ~exist('customFigSize', 'var')
    customFigSize = [500, 500].*0.7;
end

hCurPLMap = figure('visible', flagVisible, ...
    'Position', [0,0,customFigSize]);
hCurPLMap.InvertHardcopy = 'off';
hCurPLMap.Color = 'none';

hold on;

if strcmpi(colorMap, 'hot')
    curCM = colormap('hot');
    cMtoSet = curCM(end:-1:1, :);
    colormap(cMtoSet);
else
    colormap(colorMap);
end

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

[axisToSet, weightForWidth] ...
    = extendLonLatAxisByFactor( ...
    [min(areaOfInterestLons), max(areaOfInterestLons), ...
    min(areaOfInterestLats), max(areaOfInterestLats)], ...
    extensionFactor, simConfigs);

% Check color range with simulation results.
colorRange = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB;
boolsPathLossesToShow = (matRxLonLatWithPathLoss(:,3)>= colorRange(1)) ...
    & (matRxLonLatWithPathLoss(:,3)<= colorRange(2));
if min(matRxLonLatWithPathLoss(boolsPathLossesToShow,3)) ...
        >minPathLossInDbExpected
    colorRange(1) = minPathLossInDbExpected;
end

% TX.
if flagRiseTxToTop
    zsToPlotTx = ones(length(cellAntLonLats(:,1)), 1).*(colorRange(2)+1);
else
    zsToPlotTx = zeros(length(cellAntLonLats(:,1)), 1);
end

hCurHandleTxs = plot3(cellAntLonLats(:,1), cellAntLonLats(:,2), ...
    zsToPlotTx, ...
    markerTowers, ...
    'MarkerSize', markerSizeTowers, ...
    'Color', colorTowers, ...
    'LineWidth', lineWidthTowers);

hRxs = nan;
if any(boolsPathLossesToShow)
    xLabelToSet = 'Longitude';
    yLabelToSet = 'Latitude';
    cLabelToSet = 'Path Loss (dB)';
    switch lower(flagCmdToPlotPLs)
        case 'plot3k'
            hRxs = plot3k( ...
                matRxLonLatWithPathLoss(boolsPathLossesToShow, :), ...
                'Labels', {'', ...
                xLabelToSet, yLabelToSet, ...
                '', cLabelToSet}, ...
                'ColorRange', colorRange);
        case 'surf'
            set(gca, 'fontWeight', 'bold');
            
            % Create meshgrid for surf.
            upSampFactor = 10;
            sufNumPtsPerSide = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE ...
                .*upSampFactor;
            lons = matRxLonLatWithPathLoss(:, 1);
            lats = matRxLonLatWithPathLoss(:, 2);
            zs = matRxLonLatWithPathLoss(:, 3);
            % Set points not to show to NaN.
            zs(~boolsPathLossesToShow) = nan;
            
            % Find the ranges for the boundary of interet (BoI) to build a
            % new grid for showing the results.
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
            zsNew = griddata(lons,lats,zs,lonsNew,latsNew);
            
            % Replace the resulting NaN grid data points with nearest
            % neighbor's value in the original dataset.
            ziNearest = griddata(lons,lats,zs,lonsNew,latsNew,'Nearest');
            boolsZsNewIsNan = isnan(zsNew);
            zsNew(boolsZsNewIsNan) = ziNearest(boolsZsNewIsNan);
            
            % For positions with NaN values in ziNearest, zsNew should also
            % have NaN.
            zsNew(isnan(ziNearest)) = nan;
            
            % Ignore points out of the area of interest by seting the z
            % values for them to NaN.
            [in,on] = inpolygon(lonsNew(:), latsNew(:), lonsBoI, latsBoI);
            boolsPtsToIgnore = ~(in|on);
            if any(boolsPtsToIgnore)
                zsNew(boolsPtsToIgnore) = nan;
            end
            
            hRxs = surf( lonsNew,latsNew,zsNew, ...
                'FaceAlpha',0.5, 'EdgeColor', 'none');
            caxis(colorRange); xlabel(xLabelToSet); ylabel(yLabelToSet);
            hCb = colorbar; % ylabel(hCb, cLabelToSet);
            title(hCb, cLabelToSet);
            [c,h] = contour(lonsNew,latsNew,zsNew);
            clabel(c,h);
    end
    
    % Put an empty title to avoid tightfig cutting out the clabel.
    title(' ');
end

hLeg = legend(hCurHandleTxs, 'Cell towers', 'Location', LEGEND_LOC);
set(hLeg, 'color', legendBackgroundColor);
set(hCurPLMap, 'Color', 'w');

view(2); xticks([]); yticks([]);
xticks([]); yticks([]);

adjustFigSizeByContent(hCurPLMap, axisToSet, ...
    'height', weightForWidth*1.05);

if strcmpi(flagCmdToPlotPLs, 'plot3k') && ishandle(hRxs)
    plotGoogleMapAfterPlot3k(gcf, 'satellite');
else
    plot_google_map('MapType', 'satellite');
end

end
% EOF