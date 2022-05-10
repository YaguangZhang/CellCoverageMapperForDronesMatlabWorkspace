%COMPAREINSTATELIDAR Explore the difference between the new LiDAR dataset
%hosted on https://lidar.digitalforestry.org/in2019/ and the cached one.
%
% Yaguang Zhang, Purdue, 05/09/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, 'Tests');

% The absolute path to the .mat file for LiDAR meta info.
pathToLidarMetaMatFile = fullfile(pwd, '..', 'Lidar_2019', 'IN', ...
    'DSM_Frankie_PartialBackup', 'metaInfo.mat');
% The absolute path to the .mat file for a list of anomaly tiles which do
% not have imbeded projected coordinate reference system (CRS) information.
pathToAnomalyTileList = fullfile(pwd, 'lib', 'lidar', ...
    'anomalyTileNameWithProjCRS.mat');
% The absolute path to save the .kml overview file.
pathToKmlFolder = fullfile(pwd, '..', 'PostProcessingResults', ...
    'Test_LidarCoverage');
if ~exist('pathToKmlFolder', 'dir')
    mkdir(pathToKmlFolder);
end

% For plotting.
polyshapeFaceAlpha = 0.75;
polyshapeColorNormal = [0,0,1];
polyshapeColorAnomaly = [1,0,0];
deltaAltInM = 100;

%% Coverage Map for the Cached LiDAR Data

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Generating coverage map for the cached LiDAR data ...'])

% For lidarFileRelDirs, xYBoundryPolygons, and lonLatBoundryPolygons.
load(pathToLidarMetaMatFile);
% For anomalyTileNameWithProjCRS.
load(pathToAnomalyTileList);

numOfTiles = size(lidarFileRelDirs, 1);
fileNameListAllTiles = cell(numOfTiles, 1);
fileNameListAnomalyTiles = anomalyTileNameWithProjCRS(:,1);
boolsIsAnomaly = false(numOfTiles, 1);
% Save results in polyshape arrays.
polyshapesNormalTiles = repmat(polyshape, ...
    numOfTiles-length(fileNameListAnomalyTiles), 1);
polyshapesAnomalyTiles = repmat(polyshape, ...
    length(fileNameListAnomalyTiles), 1);
[normalTileCnt, anomalyTileCnt] = deal(1);
for idxTile = 1:numOfTiles
    [~, fileNameListAllTiles{idxTile}] ...
        = fileparts(lidarFileRelDirs{idxTile});
    if any(strcmp(fileNameListAllTiles{idxTile}, fileNameListAnomalyTiles))
        boolsIsAnomaly(idxTile) = true;
        polyshapesAnomalyTiles(anomalyTileCnt) ...
            = lonLatBoundryPolygons{idxTile};
        anomalyTileCnt = anomalyTileCnt+1;
    else
        polyshapesNormalTiles(normalTileCnt) ...
            = lonLatBoundryPolygons{idxTile};
        normalTileCnt = normalTileCnt+1;
    end
end

% % Merge the polygons to get the final results.
%  polyshapeNormalTiles = union(polyshapesNormalTiles);
% polyshapeAnomalyTiles = union(polyshapesAnomalyTiles);

% Double-check the results.
assert(sum(boolsIsAnomaly) == length(fileNameListAnomalyTiles), ...
    ['Anomaly records found in the LiDAR meta info', ...
    ' does not agree with the anomaly tile list!']);

% Plot the polygons.
figure; hold on;
hNormalTiles = plot(polyshapesNormalTiles, ...
    'FaceColor', polyshapeColorNormal, 'EdgeColor', 'none', ...
    'FaceAlpha', polyshapeFaceAlpha);
hAnomalyTiles = plot(polyshapesAnomalyTiles, ...
    'FaceColor', polyshapeColorAnomaly, 'EdgeColor', 'none', ...
    'FaceAlpha', polyshapeFaceAlpha);
plot_google_map('MapType', 'streetview');
legend([hNormalTiles(1), hAnomalyTiles(1)], ...
    'Valid Tiles', 'Anomaly Tiles');
xlabel('Longitude (degree)'); ylabel('Latitude (degree)');
title('Overview of Anomaly Tile Locations');

exportgraphics(gca, fullfile(pathToSaveResults, ...
    'compareInStateLidar_MapOverview.jpg'));

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% KML File

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Generating KML file for anomaly tile overview ...'])

fileNameListNormalTiles = fileNameListAllTiles(~boolsIsAnomaly);

numOfNormalTs = length(polyshapesNormalTiles);
kmlPolysForNormalTiles = cell(numOfNormalTs, 1);
for idxNormT = 1:numOfNormalTs
    kmlPolysForNormalTiles{idxNormT} = ge_poly( ...
        polyshapesNormalTiles(idxNormT).Vertices(:,1), ...
        polyshapesNormalTiles(idxNormT).Vertices(:,2), ...
        'lineWidth', 0, ...
        'polyColor', ...
        constructHexColorForKml(round(polyshapeColorNormal.*255), ...
        round(polyshapeFaceAlpha.*255)), ...
        'altitude', deltaAltInM, ...
        'altitudeMode', 'relativeToGround', ...
        'extrude', 1, ...
        'tessellate', true, ...
        'name', fileNameListNormalTiles{idxNormT});
end

numOfAnomalyTs = length(polyshapesAnomalyTiles);
kmlPolysForAnomalyTiles = cell(numOfAnomalyTs, 1);
for idxAnomT = 1:numOfAnomalyTs
    kmlPolysForAnomalyTiles{idxAnomT} = ge_poly( ...
        polyshapesAnomalyTiles(idxAnomT).Vertices(:,1), ...
        polyshapesAnomalyTiles(idxAnomT).Vertices(:,2), ...
        'lineWidth', 0, ...
        'polyColor', ...
        constructHexColorForKml(round(polyshapeColorAnomaly.*255), ...
        round(polyshapeFaceAlpha.*255)), ...
        'altitude', deltaAltInM, ...
        'altitudeMode', 'relativeToGround', ...
        'extrude', 1, ...
        'tessellate', true, ...
        'name', fileNameListAnomalyTiles{idxAnomT});
end

kmlPolyFolders = cell(2, 1);
kmlPolyFolders{1} = ge_folder('Tiles with Valid Projection Info', ...
    strcat(kmlPolysForNormalTiles{:}));
kmlPolyFolders{2} = ge_folder('Tiles without Projection Info', ...
    strcat(kmlPolysForAnomalyTiles{:}));

ge_output(fullfile(pathToKmlFolder, ...
    'compareInStateLidar_MapOverview.kml'), ...
    ge_folder('Cached LiDAR Tiles', strcat(kmlPolyFolders{:})));

% Also output a .csv file of all the file names of the anomaly tiles.
% Write the table to a CSV file
writetable(cell2table(fileNameListAnomalyTiles), ...
    fullfile(pathToKmlFolder, 'fileNameListAnomalyTiles.csv'));

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

% EOF