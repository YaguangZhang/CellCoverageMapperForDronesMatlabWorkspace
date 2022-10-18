%PLOTGEOFIGSFORTWC2021 Generate some overview plots of the LiDAR and ground
%elevation data for publication after the simulation is done.
%
% This helper script is for the IEEE TWC paper prepared in 2021.
%
% Yaguang Zhang, Purdue, 05/11/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath')));
addpath('.'); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Path to save the plots.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'GeoFigsForTwc2021');
% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% For the grid covering IN.
NUM_OF_PIXELS_FOR_LONGER_SIDE = 1000;

% IN boundary.
[inBoundaryLatLons, inBoundaryXYs, inBoundaryUtmZone] = loadInBoundary;

% For plotting: black - blue - light green.
customTurbo1 = turbo(256); % customTurbo1 = customTurbo1(141:-1:1, :);
% Light green - red - black.
customTurbo2 = turbo(256); % customTurbo2 = customTurbo2(142:end, :);

%% LiDAR Data Info

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading Indiana LiDAR meta data ...'])

% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(inBoundaryUtmZone);

pathToCacheGeoInfo = fullfile(pathToSaveResults, 'CachedGeoInfo.mat');
if exist(pathToCacheGeoInfo, 'file')
    load(pathToCacheGeoInfo);
else
    dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar_2019', 'IN', 'DSM');

    % Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need
    % to remove preprocessIndianaLidarDataSet from path after things are
    % done.
    addpath(fullfile(pwd, 'lib', 'lidar'));
    [lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
        = preprocessLidarDataSetDsm(dirToLidarFiles, ...
        deg2utm_speZone, utm2deg_speZone);
    rmpath(fullfile(pwd, 'lib', 'lidar'));
    lidarFileAbsDirs = cellfun(@(d) ...
        [dirToLidarFiles, strrep(d, '\', filesep)], ...
        lidarFileRelDirs, 'UniformOutput', false);

    % Extra information on the LiDAR data set.
    %   - Overall boundries for the area covered by the LiDAR data set in
    %   UTM.
    % lidarFilesXYCoveragePolyshape ...
    %     = mergePolygonsForAreaOfInterest(lidarFileXYCoveragePolyshapes,
    %     1);
    %   - Centroids for the LiDAR files in UTM.
    lidarFileXYCentroids ...
        = extractCentroidsFrom2DPolyCell(lidarFileXYCoveragePolyshapes);
    %   - The .mat copies for the LiDAR data. For the 2019 dataset, they
    %   are stored in a cache folder.
    lidarMatFileAbsDirs = lidarFileAbsDirs;
    for idxMatF = 1:length(lidarMatFileAbsDirs)
        [lidarMatFPath, lidarMatFName, ~] ...
            = fileparts(lidarMatFileAbsDirs{idxMatF});
        lidarMatFileAbsDirs{idxMatF} = fullfile(lidarMatFPath, '..', ...
            'MatlabCache', [lidarMatFName, '.mat']);
    end

    % Grid covering IN.
    [indianaGridXYPts, indianaGridResolutionInM] ...
        = buildSimGrid(inBoundaryXYs, NUM_OF_PIXELS_FOR_LONGER_SIDE);
    % Ground elevation and LiDAR z overview for shrunk IN.
    [terrainEles, lidarZs, curEleForNanPts] ...
        = generateProfileSamps( ...
        indianaGridXYPts, utm2deg_speZone, ...
        lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
        lidarMatFileAbsDirs, 'both');

    [indianaGridLats, indianaGridLons] = utm2deg_speZone( ...
        indianaGridXYPts(:,1), indianaGridXYPts(:,2));
    indianaGridLatLonPts = [indianaGridLats, indianaGridLons];

    save(pathToCacheGeoInfo, 'indianaGridLatLonPts', ...
        'indianaGridXYPts', 'indianaGridResolutionInM', ...
        'terrainEles', 'lidarZs', 'curEleForNanPts');
end

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Plots

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating overview plots ...'])

% Replace -inf LiDAR z values with NaN for plotting.
lidarZs(isinf(lidarZs)) = nan;

curFigSize = [400, 500];
curFigAxis = [-88.331812, -84.544638, 37.641598, 41.887921];
curColormap = turbo;
curColorRange = [85, 705];

% Ground elevation.
figure('Position', [0, 0, curFigSize]); colormap(curColormap);
plot3k([indianaGridLatLonPts(:,[2,1]), terrainEles], ...
    'Labels', {'', 'Longitude (degree)', 'Latitude (degree)', '', ...
    'Ground Elevation (m)'}, 'ColorRange',curColorRange);
view(2); zlim([0, max(terrainEles)]); axis(curFigAxis);
plot_google_map('MapType', 'roadmap');

saveas(gcf, fullfile(pathToSaveResults, 'IndianaGroundEle_plot3k.jpg'));

% LiDAR.
figure('Position', [0, 0, curFigSize]); colormap(curColormap);
plot3k([indianaGridLatLonPts(:,[2,1]), lidarZs], ...
    'Labels', {'', 'Longitude (degree)', 'Latitude (degree)', '', ...
    'LiDAR z (m)'}, 'ColorRange',curColorRange);
view(2); zlim([0, max(lidarZs)]); axis(curFigAxis);
plot_google_map('MapType', 'roadmap');

saveas(gcf, fullfile(pathToSaveResults, 'IndianaLidarZ_plot3k.jpg'));

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Plots for Publication

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating figures for publication ...'])

curFigSize = [270, 300];
curFigAxis = [-88.331812, -84.544638, 37.584357, 41.945162];
curColormap = customTurbo1;
curColorRange = [85, 705];
curCbTitlePos = [-9.52501, 189.15002, 0];
epsResolutionDpi = 300;

[simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(inBoundaryUtmZone);

% Ground elevation.
figure('Position', [0, 0, curFigSize]); hold on; colormap(curColormap);
xticks([]); yticks([]);
hPolyIn = plot3(inBoundaryLatLons(:, 2), inBoundaryLatLons(:, 1), ...
    ones(size(inBoundaryLatLons(:, 2))).*curColorRange(2), ...
    '-.', 'LineWidth', 2.3, 'Color', 'k');
axis(curFigAxis);
plot_google_map('MapType', 'roadmap'); axis manual;
[~, ~, hCb] = plot3k([indianaGridLatLonPts(:,[2,1]), terrainEles], ...
    'Labels', {'', '', '', '', ...
    'Ground Elevation (m)'}, 'ColorRange', curColorRange, ...
    'Marker', {'.', 1}, 'CBLabels', 6);
hCb.Title.Position = curCbTitlePos;
hLeg = legend(hPolyIn, 'Indiana', 'Location', 'southeast');
view(2); zlim([0, curColorRange(2)]); axis(curFigAxis);

curFigName = fullfile(pathToSaveResults, ...
    'IndianaGroundEle_plot3k_fixedCR');
% saveEpsFigForPaper(gcf, curFigName);
title(' '); tightfig;
curFigNameExpGraph = [curFigName, '_expGraph'];

hLeg.Position = [1.9855, 0.1606, 2.4871, 0.5027];

exportgraphics(gca, [curFigNameExpGraph, '.jpg'], ...
    'Resolution', epsResolutionDpi);
exportgraphics(gca, [curFigNameExpGraph, '.eps'], ...
    'Resolution', epsResolutionDpi);

% LiDAR.
figure('Position', [0, 0, curFigSize]); hold on; colormap(curColormap);
xticks([]); yticks([]);
plot3(inBoundaryLatLons(:, 2), inBoundaryLatLons(:, 1), ...
    ones(size(inBoundaryLatLons(:, 2))).*curColorRange(2), ...
    '-.', 'LineWidth', 2.3, 'Color', 'k');
axis(curFigAxis);
plot_google_map('MapType', 'roadmap'); axis manual;
[~, ~, hCb] = plot3k([indianaGridLatLonPts(:,[2,1]), lidarZs], ...
    'Labels', {'', '', '', '', ...
    'LiDAR z (m)'}, 'ColorRange', curColorRange, ...
    'Marker', {'.', 1}, 'CBLabels', 6);
hCb.Title.Position = curCbTitlePos;
view(2); zlim([0, curColorRange(2)]); axis(curFigAxis);

curFigName = fullfile(pathToSaveResults, 'IndianaLidarZ_plot3k_fixedCR');
% saveEpsFigForPaper(gcf, curFigName);
title(' '); tightfig;
curFigNameExpGraph = [curFigName, '_expGraph'];
exportgraphics(gca, [curFigNameExpGraph, '.jpg'], ...
    'Resolution', epsResolutionDpi);
exportgraphics(gca, [curFigNameExpGraph, '.eps'], ...
    'Resolution', epsResolutionDpi);

% Ground elevation - tighter color range.
curTightColorRange = [105, 365];
figure('Position', [0, 0, curFigSize]); hold on; colormap(curColormap);
xticks([]); yticks([]);
hPolyIn = plot3(inBoundaryLatLons(:, 2), inBoundaryLatLons(:, 1), ...
    ones(size(inBoundaryLatLons(:, 2))).*curTightColorRange(2), ...
    '-.', 'LineWidth', 2.3, 'Color', 'k');
axis(curFigAxis);
plot_google_map('MapType', 'roadmap'); axis manual;
[~, ~, hCb] = plot3k([indianaGridLatLonPts(:,[2,1]), terrainEles], ...
    'Labels', {'', '', '', '', ...
    'Ground Elevation (m)'}, 'ColorRange', curTightColorRange, ...
    'Marker', {'.', 1}, 'CBLabels', 6);
hCb.Title.Position = curCbTitlePos;
hLeg = legend(hPolyIn, 'Indiana', 'Location', 'southeast');
view(2); zlim([0, curTightColorRange(2)]); axis(curFigAxis);

curFigName = fullfile(pathToSaveResults, ...
    'IndianaGroundEle_plot3k');
% saveEpsFigForPaper(gcf, curFigName);
title(' '); tightfig;
curFigNameExpGraph = [curFigName, '_expGraph'];

hLeg.Position = [1.9855, 0.1606, 2.4871, 0.5027];

exportgraphics(gca, [curFigNameExpGraph, '.jpg'], ...
    'Resolution', epsResolutionDpi);
exportgraphics(gca, [curFigNameExpGraph, '.eps'], ...
    'Resolution', epsResolutionDpi);

strFormatter = '%.2f';
curVs = terrainEles;
curVs(isnan(curVs)) = [];
curMean = mean(curVs);
curStd = std(curVs);
text(-86.51167607294693, 37.928803183673466, 400, ...
    ['Mean: ', num2str(curMean, strFormatter), ' m']);
text(-88.22064572289156, 37.733011938775505, 400, ...
    ['Standard Deviation: ', num2str(curStd, strFormatter), ' m']);

hLeg.Position = [0.159875, 6.11372500000001, 2.4871, 0.5027];

exportgraphics(gca, [curFigNameExpGraph, '_withStats.jpg'], ...
    'Resolution', epsResolutionDpi);
exportgraphics(gca, [curFigNameExpGraph, '_withStats.eps'], ...
    'Resolution', epsResolutionDpi);

% LiDAR z - ground ele.
curColormap = customTurbo2;
figure; ecdf(lidarZs-terrainEles); grid on; grid minor;
saveas(gcf, fullfile(pathToSaveResults, ...
    'IndianaDiffLidarZMinusEle_ECDF.jpg'));

curTightestColorRange = [-5, 40];
curCbTitlePos = [-36.5250, curCbTitlePos(2:3)];
figure('Position', [0, 0, curFigSize]); hold on; colormap(curColormap);
hInB = plot3(inBoundaryLatLons(:, 2), inBoundaryLatLons(:, 1), ...
    ones(size(inBoundaryLatLons(:, 2))).*(curTightestColorRange(2)-1), ...
    '-.', 'LineWidth', 2.3, 'Color', 'k');
axis(curFigAxis);
plot_google_map('MapType', 'roadmap', 'Alpha', 0); axis manual;
[~, ~, hCb] = plot3k( ...
    [indianaGridLatLonPts(:,[2,1]), lidarZs-terrainEles], ...
    'Labels', {'', '', '', '', ...
    'LiDAR z minus Ground Elevation (m)'}, ...
    'ColorRange', curTightestColorRange, ...
    'Marker', {'.', 1}, 'CBLabels', 6);
hCb.Title.Position = curCbTitlePos;
view(2); zlim(curTightestColorRange); axis(curFigAxis);
delete(hInB); grid on; grid minor; xticklabels({}); yticklabels({});
% Add symbols ≤ and ≥ to color bar.
hCb.TickLabels = strtrim(hCb.TickLabels);
hCb.TickLabels{1} = ['≤', hCb.TickLabels{1}];
hCb.TickLabels{end} = ['≥', hCb.TickLabels{end}];
hCb.TickLabels = pad(hCb.TickLabels, 'left');
hCb.TickLabels{end} = [hCb.TickLabels{end}, ' '];

curFigName = fullfile(pathToSaveResults, ...
    'IndianaDiffLidarZMinusEle_plot3k');
% saveEpsFigForPaper(gcf, curFigName);
title(' '); tightfig;
%  curFigNameExpGraph = [curFigName, '_expGraph'];
% exportgraphics(gca, [curFigNameExpGraph, '.jpg'], ...
%     'Resolution', epsResolutionDpi);
% exportgraphics(gca, [curFigNameExpGraph, '.eps'], ...
%     'Resolution', epsResolutionDpi);
saveas(gcf, [curFigName, '.eps'], 'epsc');
saveas(gcf, [curFigName, '.jpg']);

strFormatter = '%.2f';
curVs = lidarZs-terrainEles;
curVs(isnan(curVs)) = [];
curMean = mean(curVs);
curStd = std(curVs);
text(-86.39550162650602,37.928803183673466, ...
    ['Mean: ', num2str(curMean, strFormatter), ' m']);
text(-88.22064572289156, 37.733011938775505, ...
    ['Standard Deviation: ', num2str(curStd, strFormatter), ' m']);

saveas(gcf, [curFigName, '_withStats.eps'], 'epsc');
saveas(gcf, [curFigName, '_withStats.jpg']);

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Cleanup

close all;

% EOF