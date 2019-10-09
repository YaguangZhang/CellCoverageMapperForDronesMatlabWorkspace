function [hCurPLMap, hCurHandleTxs] ...
    = plotPathLossMap(matRxLonLatWithPathLoss, cellAntLonLats, ...
    simConfigs, flagVisible, flagZoomIn, flagCmdToPlotPLs)
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
%         Color mapping limits for plot4k. A two element vector specifying
%         the values that map to the first and last colors. This is useful
%         for generating a series of plots with identical coloring. The
%         colormap but not the colorbar) is flipped upside down if
%         'ColorRange' is given as [max min] instead of [min max].
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
%
% Outputs:
%   - hCurPLMap
%     The handle to the resultant figure.
%   - hCurHandleTxs
%     The handle to the cellular tower.
%
% Yaguang Zhang, Purdue, 10/02/2019

% We support: 'plot3k' and 'surf'(default).
if ~exist('flagCmdToPlotPLs', 'var')
    flagCmdToPlotPLs = 'surf';
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

colorRange = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB;

hCurPLMap = figure('visible', flagVisible);
hold on;

% Area of interest.
[areaOfInterestLats, areaOfInterestLons] = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
plot(polyshape(areaOfInterestLons, areaOfInterestLats), ...
    'FaceColor', 'white');

if flagZoomIn
    axis tight;
    zoomInAxisToSet = axis;
    axis auto;
end

% TX.
hCurHandleTxs = plot3(cellAntLonLats(:,1), cellAntLonLats(:,2), ...
    ones(length(cellAntLonLats(:,1)), 1)...
    .*(colorRange(2)+1), ...
    'xr', 'LineWidth', 1.5);

% Simulation results.
boolsPathLossesToShow = (matRxLonLatWithPathLoss(:,3)>= colorRange(1)) ...
    & (matRxLonLatWithPathLoss(:,3)<= colorRange(2));

hRxs = nan;
if any(boolsPathLossesToShow)
    xLabelToSet = 'Longitude (degrees)';
    yLabelToSet = 'Latitude (degrees)';
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
            sufNumPtsPerSide = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE;
            x = matRxLonLatWithPathLoss(boolsPathLossesToShow, 1);
            y = matRxLonLatWithPathLoss(boolsPathLossesToShow, 2);
            z = matRxLonLatWithPathLoss(boolsPathLossesToShow, 3);
            [xi,yi] = meshgrid( ...
                linspace(min(x), max(x), sufNumPtsPerSide), ...
                linspace(min(y), max(y), sufNumPtsPerSide));
            zi = griddata(x,y,z,xi,yi);
            
            maskMapEdge = zeros(size(zi));
            maskMapEdge(1,:) = 1; maskMapEdge(end,:) = 1; 
            maskMapEdge(:,1) = 1; maskMapEdge(:,end) = 1; 
            maskNanMapEdgePts = isnan(zi)&maskMapEdge;
            
            ziNearest = griddata(x,y,z,xi,yi, 'Nearest');
            zi(maskNanMapEdgePts) = ziNearest(maskNanMapEdgePts);
            
            hRxs = surf( xi,yi,zi, ...
                'FaceAlpha',0.5, 'EdgeColor', 'none');
            caxis(colorRange); xlabel(xLabelToSet); ylabel(yLabelToSet);
            hCB = colorbar; ylabel(hCB, cLabelToSet);
            [c,h] = contour(xi,yi,zi);
            clabel(c,h);
    end
end

legend(hCurHandleTxs, 'TXs', 'Location', 'SouthEast'); view(2);
if flagZoomIn
    axis(zoomInAxisToSet);
end

if strcmpi(flagCmdToPlotPLs, 'plot3k') && ishandle(hRxs)
    plotGoogleMapAfterPlot3k(gcf, 'satellite');
else
    plot_google_map('MapType', 'satellite');
end

end
% EOF