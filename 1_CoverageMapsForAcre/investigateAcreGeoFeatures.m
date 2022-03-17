% INVESTIGATEACREGEOFEATURES Visualize geographical features of interest
% for ACRE.
%
% Yaguang Zhang, Purdue, 03/16/2022

clearvars -except PRESETS CARRIER_FREQUENCIES_IN_MHZ ...
    PRESET CARRIER_FREQUENCY_IN_MHZ pathToSaveSimManDiary idxFre;
clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', 'AcreGeoFeatureInvestigation');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

% We will load cached results for faster processing if available.
dirToCachedResults = fullfile(pathToSaveResults, 'cachedResults.mat');
if exist(dirToCachedResults, 'file')
    disp(' ')
    disp(['[', datestr(now, datetimeFormat), ...
        '] Loading cached geo info for ACRE ...'])

    load(dirToCachedResults);
else
    %% Load Boundary
    disp(' ')
    disp(['[', datestr(now, datetimeFormat), ...
        '] Extracting geo info for ACRE ...'])

    absDirToAcreKmz = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar', 'ACRE', 'AcreExactBoundaryRaw', 'ACRE.kmz');

    % Boundary of ACRE.
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Loading outlines of ACRE fields ...'])
    [UTM_X_Y_BOUNDARY_ACRE, UTM_ZONE_ACRE, acreKmzStruct, ...
        utmXYAcrePolygons, lonLatAcrePolygons] ...
        = extractBoundaryFromKmzFile(absDirToAcreKmz);

    % For GPS and UTM conversions.
    [deg2utm_speZone, utm2deg_speZone] ...
        = genUtmConvertersForFixedZone(UTM_ZONE_ACRE);

    UTM_LAT_LON_BOUNDARY_ACRE = nan(size(UTM_X_Y_BOUNDARY_ACRE));
    [UTM_LAT_LON_BOUNDARY_ACRE(:,1), UTM_LAT_LON_BOUNDARY_ACRE(:,2)] ...
        = utm2deg_speZone(UTM_X_Y_BOUNDARY_ACRE(:,1), ...
        UTM_X_Y_BOUNDARY_ACRE(:,2));

    %% Create a Grid

    % Corresponds to the 5 feet resolution of the LiDAR DSM.
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Creating a grid ...'])
    gridResolutionInM = 1.5;
    gridUtmXYPts = buildSimGrid(UTM_X_Y_BOUNDARY_ACRE, ...
        gridResolutionInM, true);

    %% LiDAR Data Info
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Loading Indiana LiDAR meta data ...'])

    dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar_2019', 'IN', 'DSM');

    % Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need
    % to remove preprocessIndianaLidarDataSet from path after things are
    % done.
    addpath(fullfile(pwd, 'lib', 'lidar'));
    [lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
        = preprocessIndianaLidarDataSetDsm(dirToLidarFiles, ...
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

    %% Fetch LiDAR Data
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Interpolating elevation and LiDAR data for the grid pts ...'])

    [terrainEles, lidarZs, curEleForNanPts] ...
        = generateProfileSamps( ...
        gridUtmXYPts, utm2deg_speZone, ...
        lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
        lidarMatFileAbsDirs, 'both');

    %% Fetch Cell Towers
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Fetching tower locations ...'])

    % Default to the NTIA+HIFLD cell tower locations.
    ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
        ABS_PATH_TO_SHARED_FOLDER, ...
        'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs', ...
        'NtiaLayoutMergedWithHifldCellTs_Threshold_1000m_LatLonH.csv');

    % Note: we use "height" to indicate the vertical distance from the
    % ground to the antenna; "elevation" to indicate the ground elevation;
    % and "altitude" to indicate elevation+height.
    cellAntsLatLonH = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);

    %% Cache the Results
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Caching results ...'])

    save(dirToCachedResults, ...
        ... ACRE field boudnaries:
        'UTM_X_Y_BOUNDARY_ACRE', 'UTM_ZONE_ACRE', 'acreKmzStruct', ...
        'utmXYAcrePolygons', 'lonLatAcrePolygons', ...
        ... Grid point locs:
        'gridResolutionInM', 'gridUtmXYPts', ...
        ... Ground elevation and LiDAR z values.
        'terrainEles', 'lidarZs', ...
        ... Cell tower antenna locations.
        'cellAntsLatLonH');
end

if ~exist('deg2utm_speZone', 'var')
    % For GPS and UTM conversions.
    [deg2utm_speZone, utm2deg_speZone] ...
        = genUtmConvertersForFixedZone(UTM_ZONE_ACRE);
end

gridLatLonPts = nan(size(gridUtmXYPts));
[gridLatLonPts(:,1), gridLatLonPts(:,2)] ...
    = utm2deg_speZone(gridUtmXYPts(:,1), gridUtmXYPts(:,2));

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Plots

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Generating figures for ground elevation and LiDAR z values ...'])

figureSize = [500, 500];
% Use the same color for the same height so that figures can be compared.
commonCaxisRange = [205, 400];

% For reusing path loss map plot functions.
extensionFactor= 0.2;
simConfigs.deg2utm_speZone = deg2utm_speZone;
simConfigs.utm2deg_speZone = utm2deg_speZone;

simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = UTM_X_Y_BOUNDARY_ACRE;
simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE = 100;
simConfigs.CURRENT_SIMULATION_TAG = 'AcreGeo';
simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB = [-inf, inf];

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Ground elevation point cloud ...'])

figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); colormap turbo;
plot3k([gridLatLonPts(:, [2,1]), terrainEles], ...
    'Labels', {'', '', '', '', ...
    'Ground Elevation (m)            '});
xticks([]); yticks([]); view(2); zlim([0, max(terrainEles)]);
[axisToSet, weightForWidth] ...
    = extendLonLatAxisByFactor( ...
    [min(gridLatLonPts(:,2)), max(gridLatLonPts(:,2)), ...
    min(gridLatLonPts(:,1)), max(gridLatLonPts(:,1))], ...
    extensionFactor, simConfigs);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, '1_GroundEle_GridPts.jpg'), '-m3');

figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); colormap turbo;
plot3k([gridLatLonPts(:, [2,1]), terrainEles], ...
    'Labels', {'', '', '', '', ...
    'Ground Elevation (m)            '}, 'ColorRange', commonCaxisRange);
xticks([]); yticks([]); view(2); zlim([0, max(terrainEles)]);
[axisToSet, weightForWidth] ...
    = extendLonLatAxisByFactor( ...
    [min(gridLatLonPts(:,2)), max(gridLatLonPts(:,2)), ...
    min(gridLatLonPts(:,1)), max(gridLatLonPts(:,1))], ...
    extensionFactor, simConfigs);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, ...
    '1_GroundEle_GridPts_FixedColorBar.jpg'), '-m3');

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] LiDAR z point cloud ...'])

% Need to remove the -inf values for water before plotting.
lidarZsToPlot = lidarZs;
lidarZsToPlot(isinf(lidarZs)) = nan;

figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); colormap turbo;
plot3k([gridLatLonPts(:, [2,1]), lidarZsToPlot], ...
    'Labels', {'', '', '', '', ...
    'LiDAR z (m)'});
xticks([]); yticks([]); view(2); zlim([0, max(lidarZsToPlot)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, '1_LiDARZ_DSM_GridPts.jpg'), '-m3');

figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); colormap turbo;
plot3k([gridLatLonPts(:, [2,1]), lidarZsToPlot], ...
    'Labels', {'', '', '', '', ...
    'LiDAR z (m)'}, 'ColorRange', commonCaxisRange);
xticks([]); yticks([]); view(2); zlim([0, max(lidarZsToPlot)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, ...
    '1_LiDARZ_DSM_GridPts_FixedColorBar.jpg'), '-m3');
refAxis = axis;

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] LiDAR z point - ground elevation cloud ...'])

lidarZMinusEle = lidarZsToPlot-terrainEles;

figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); colormap turbo;
plot3k([gridLatLonPts(:, [2,1]), lidarZMinusEle], ...
    'Labels', {'', '', '', '', ...
    'LiDAR z - Ground Elevation (m)                            '});
xticks([]); yticks([]); view(2); zlim([0, max(lidarZMinusEle)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
axis(refAxis);
export_fig(fullfile(pathToSaveResults, ...
    '1_LiDARZMinusEle_GridPts.jpg'), '-m3');

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Ground elevation contour plot ...'])

plotPathLossMap( ...
    [gridLatLonPts(:, [2,1]), terrainEles], [], simConfigs);
hCb = findall(gcf,'type','ColorBar');
title(hCb, 'Ground Elevation (m)            ');
set(gcf, 'Position', [0, 0, figureSize]);
xlabel(''); ylabel(''); caxis([min(terrainEles), max(terrainEles)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, '2_GroundEle_Contour.jpg'), '-m3');
caxis(commonCaxisRange);
export_fig(fullfile(pathToSaveResults, ...
    '2_GroundEle_Contour_FixedColorBar.jpg'), '-m3');

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] LiDAR z contour plot ...'])

plotPathLossMap( ...
    [gridLatLonPts(:, [2,1]), lidarZsToPlot], [], simConfigs);
hCb = findall(gcf,'type','ColorBar');
title(hCb, 'LiDAR z (m)');
set(gcf, 'Position', [0, 0, figureSize]);
xlabel(''); ylabel(''); caxis([min(lidarZsToPlot), max(lidarZsToPlot)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, '2_LiDARZ_DSM_Contour.jpg'), '-m3');
caxis(commonCaxisRange);
export_fig(fullfile(pathToSaveResults, ...
    '2_LiDARZ_DSM_Contour_FixedColorBar.jpg'), '-m3');

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] LiDAR z point - ground elevation contour plot ...'])

plotPathLossMap( ...
    [gridLatLonPts(:, [2,1]), lidarZMinusEle], [], simConfigs);
hCb = findall(gcf,'type','ColorBar');
title(hCb, ...
    'LiDAR z - Ground Elevation (m)                                 ');
set(gcf, 'Position', [0, 0, figureSize]);
xlabel(''); ylabel(''); caxis([min(lidarZMinusEle), max(lidarZMinusEle)]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
hGM = plot_google_map('MapType', 'satellite'); delete(hGM);
export_fig(fullfile(pathToSaveResults, ...
    '2_LiDARZMinusEle_Contour.jpg'), '-m3');

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Ground elevation histograms ...'])
numsOfBins = [3:7, 30, 50, 100];
for idxNoB = 1:length(numsOfBins)
    numOfBins = numsOfBins(idxNoB); figure;
    histogram(terrainEles, numOfBins, 'Normalization', 'probability');
    xlabel('Ground Elevation (m)');
    ylabel(['Relative Frequency (', num2str(numOfBins), ' Bins)']);
    grid on; grid minor;
    export_fig(fullfile(pathToSaveResults, ...
        ['3_GroundEle_Histogram_NumOfBins_', num2str(numOfBins), ...
        '.jpg']), '-m3');
end

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] LiDAR z histograms ...'])
for idxNoB = 1:length(numsOfBins)
    numOfBins = numsOfBins(idxNoB); figure;
    histogram(lidarZsToPlot, numOfBins, 'Normalization', 'probability');
    xlabel('LiDAR z (m)');
    ylabel(['Relative Frequency (', num2str(numOfBins), ' Bins)']);
    grid on; grid minor;
    export_fig(fullfile(pathToSaveResults, ...
        ['3_LiDARZ_DSM_Histogram_NumOfBins_', num2str(numOfBins), ...
        '.jpg']), '-m3');
end

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Ground elevation: group maps ...'])
numsOfGroups = numsOfBins;
for idxNoG = 1:length(numsOfGroups)
    numOfGroups = numsOfGroups(idxNoG);

    % Note: The value X(i) is in the kth bin if edges(k) ≤ X(i) <
    % edges(k+1). The last bin also includes the right bin edge, so that it
    % contains X(i) if edges(end-1) ≤ X(i) ≤ edges(end).
    [~, edges] = histcounts(terrainEles, numOfGroups);

    numOfEdges = length(edges);
    assert(numOfEdges-1==numOfGroups, ...
        'Number of edges does not match with number of groups!');

    terrainElesToShow = terrainEles;
    for idxOfG = 1:numOfGroups
        startValue = edges(idxOfG);
        endValue = edges(idxOfG+1);
        centerValue = mean([startValue, endValue]);
        terrainElesToShow((terrainEles>=startValue) ...
            & (terrainEles<endValue)) = centerValue;
    end
    % The elevation values which equal to the last edge value need to be
    % moved to the last bin center.
    terrainElesToShow(terrainEles==edges(end)) = centerValue;

    figure('Position', [0, 0, figureSize], ...
        'PaperPositionMode', 'auto'); colormap turbo;
    plot3k([gridLatLonPts(:, [2,1]), terrainElesToShow], ...
        'ColorBar', false);
    xticks([]); yticks([]); view(2); zlim([0, max(terrainEles)]);
    title([num2str(numOfGroups), ' Groups']);
    adjustFigSizeByContent(gcf, axisToSet, ...
        'height', weightForWidth*1.05);
    plot_google_map('MapType', 'satellite');
    export_fig(fullfile(pathToSaveResults, ...
        ['4_GroundEle_Grouped_NumOfGroups_', num2str(numOfGroups), ...
        '.jpg']), '-m3');
end

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Overview of cell tower locations ...'])
figure('Position', [0, 0, figureSize], ...
    'PaperPositionMode', 'auto'); hold on;
hAcre = plot(lonLatAcrePolygons);
xticks([]); yticks([]);
adjustFigSizeByContent(gcf, axisToSet, ...
    'height', weightForWidth*1.05);
plot_google_map('MapType', 'satellite');
axis manual;
hTowers = plot(cellAntsLatLonH(:,2), cellAntsLatLonH(:,1), ...
    'vr', 'LineWidth', 2);
legend(hTowers, 'Cell Towers');
export_fig(fullfile(pathToSaveResults, ...
    '5_Overview_CellTowers.jpg'), '-m3');

[axisToSetZoomedOut, weightForWidthZoomedOut] ...
    = extendLonLatAxisByFactor( ...
    [min(gridLatLonPts(:,2)), max(gridLatLonPts(:,2)), ...
    min(gridLatLonPts(:,1)), max(gridLatLonPts(:,1))], ...
    extensionFactor*30, simConfigs);
adjustFigSizeByContent(gcf, axisToSetZoomedOut, ...
    'height', weightForWidthZoomedOut*1.05);
plot_google_map('MapType', 'satellite');
export_fig(fullfile(pathToSaveResults, ...
    '5_Overview_CellTowers_ZoomedOut.jpg'), '-m3');

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

% EOF