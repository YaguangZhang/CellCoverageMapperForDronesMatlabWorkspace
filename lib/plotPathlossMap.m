function [hCurPLMap, hCurHandleTxs] ...
    = plotPathLossMap(matRxLonLatWithPathLoss, cellAntLonLats, ...
    simConfigs, flagVisible)
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
%
% Outputs:
%   - hCurPLMap
%     The handle to the resultant figure.
%   - hCurHandleTxs
%     The handle to the cellular tower.
%
% Yaguang Zhang, Purdue, 10/02/2019

% By default, show the plot.
if ~exist('flagVisible', 'var')
    flagVisible = true;
end

colorRange = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB;

hCurPLMap = figure('visible', flagVisible);
hold on;

% Area of interest.
[areaOfInterestLats, areaOfInterestLons] = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
plot(polyshape(areaOfInterestLons, areaOfInterestLats));

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
    hRxs = plot3k(matRxLonLatWithPathLoss(boolsPathLossesToShow, :), ...
        'Labels', {'', ...
        'Longitude (degrees)', 'Latitude (degrees)', ...
        '', 'Path Loss (dB)'}, ...
        'ColorRange', colorRange);
end
legend(hCurHandleTxs, 'TXs', 'Location', 'SouthEast');
if ishandle(hRxs)
    plotGoogleMapAfterPlot3k(gcf, 'satellite');
else
    plot_google_map('MapType', 'satellite');
end
view(2);

end
% EOF