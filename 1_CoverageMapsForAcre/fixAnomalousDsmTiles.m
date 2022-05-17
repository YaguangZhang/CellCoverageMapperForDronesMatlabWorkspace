%FIXANOMALOUSDSMTILES Add projection information to the anomalous DSM tiles
%listed in lib/lidar/anomalyTileNameWithProjCRS.mat.
%
% Note: this script needs to be run with a complete local LiDAR data cache.
%
% Yaguang Zhang, Purdue, 05/12/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath')));
addpath('.'); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'AnomalousDsmTilesRepairment');
pathToSaveRepairedTiles = fullfile(pathToSaveResults, 'RepairedDsmTiles');
pathToSaveDebugFigs = fullfile(pathToSaveResults, 'FigsForDebugging');
% Create directories if necessary.
if exist(pathToSaveRepairedTiles, 'dir')~=7
    mkdir(pathToSaveRepairedTiles);
end
if exist(pathToSaveDebugFigs, 'dir')~=7
    mkdir(pathToSaveDebugFigs);
end

% Path to anomalous tile list.
dirToAnomalousTileList = fullfile(pwd, 'lib', 'lidar', ...
    'anomalyTileNameWithProjCRS.mat');
dirToAnomalousTileList_Debug = fullfile(pwd, 'lib', 'lidar', ...
    'anomalyTileNameWithProjCRS_Debug.mat');

% Path to DSM LiDAR tiles.
dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM');

% Turn the diary logging function on.
dirToDiaryFile = fullfile(pathToSaveResults, 'diary.txt');
% Get rid of the old diary if it exists.
if exist(dirToDiaryFile, 'file')
    delete(dirToDiaryFile);
end
diary(dirToDiaryFile);

%% Load Cached Results

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Loading cached results ...'])

disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading the list of anomalous tiles ...'])

load(dirToAnomalousTileList, 'anomalyTileNameWithProjCRS');
load(dirToAnomalousTileList_Debug, 'anomalyTileNameWithProjCRS_Debug');

disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading Indiana LiDAR meta data ...'])

% For GPS and UTM conversions.
[inBoundaryLatLons, inBoundaryXYs, inBoundaryUtmZone] = loadInBoundary;
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(inBoundaryUtmZone);

% Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need to
% remove preprocessIndianaLidarDataSet from path after things are done.
addpath(fullfile(pwd, 'lib', 'lidar'));
[lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);
rmpath(fullfile(pwd, 'lib', 'lidar'));
lidarFileAbsDirs = cellfun(@(d) ...
    [dirToLidarFiles, strrep(d, '\', filesep)], ...
    lidarFileRelDirs, 'UniformOutput', false);

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Fix the Tiles One by One

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Reparing anomalous tiles ...'])

disp(['    [', datestr(now, datetimeFormat), ...
    '] Reading and rewriting LiDAR data ...'])

numOfAnomalousTiles = size(anomalyTileNameWithProjCRS, 1);
for idxT = 1:numOfAnomalousTiles
    curAnomalousTileName = anomalyTileNameWithProjCRS{idxT, 1};
    curIdxTileDir = find(contains(lidarFileAbsDirs, curAnomalousTileName));

    assert(length(curIdxTileDir)==1, ...
        ['Exactly one raw .tif LiDAR file is expected for tile "' ...
        curAnomalousTileName, '"! ', ...
        num2str(length(curIdxTileDir)), ' found, instead.'])

    [tileContent, curProj] ...
        = readgeoraster(lidarFileAbsDirs{curIdxTileDir});
    curProjRef = anomalyTileNameWithProjCRS_Debug{idxT, 2};

    switch curProjRef.ProjectedCRS.Name
        case 'NAD83(HARN) / Indiana East (ftUS)'
            coordRefSysCode = 2967;
        case 'NAD83(HARN) / Indiana West (ftUS)'
            coordRefSysCode = 2968;
        case 'unknown'
            curProjParameters ...
                = curProjRef.ProjectedCRS.ProjectionParameters;

            if curProjParameters.LatitudeOfNaturalOrigin == 37.5 ...
                    && curProjParameters.LongitudeOfNaturalOrigin ...
                    == -85.6666666666667 ...
                    && curProjParameters.ScaleFactorAtNaturalOrigin ...
                    == 0.999966666666667 ...
                    && curProjParameters.FalseEasting == 100000 ...
                    && curProjParameters.FalseNorthing == 250000
                % 'NAD83(HARN) / Indiana East (ftUS)'
                coordRefSysCode = 2967;
            elseif curProjParameters.LatitudeOfNaturalOrigin == 37.5 ...
                    && curProjParameters.LongitudeOfNaturalOrigin ...
                    == -87.0833333333333 ...
                    && curProjParameters.ScaleFactorAtNaturalOrigin ...
                    == 0.999966666666667 ...
                    && curProjParameters.FalseEasting == 2952750 ...
                    && curProjParameters.FalseNorthing == 820208.333333333
                % 'NAD83(HARN) / Indiana West (ftUS)'
                coordRefSysCode = 2968;
            elseif curProjParameters.LatitudeOfNaturalOrigin == 37.5 ...
                    && curProjParameters.LongitudeOfNaturalOrigin ...
                    == -87.0833333333333 ...
                    && curProjParameters.ScaleFactorAtNaturalOrigin ...
                    == 0.999966666666667 ...
                    && curProjParameters.FalseEasting == 900000 ...
                    && curProjParameters.FalseNorthing == 250000
                % 'NAD83(NSRS2007) / Indiana West': 3534. However, it seems
                % the data only seem right with code 2968.
                coordRefSysCode = 2968;
            else
                error('Unexpected projection parameters!')
            end
        otherwise
            error('Unexpected projection name!')
    end
    geotiffwrite(fullfile(pathToSaveRepairedTiles, ...
        [curAnomalousTileName, '.tif']), tileContent, curProj, ...
        'CoordRefSysCode', coordRefSysCode);
end

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Debug Figs

% Any LiDAR z value too big or too small will be discarded (set to NaN).
maxAllowedAbsLidarZ = 10^38;
flagGenFigsQuietly = true;

for idxT = 1:numOfAnomalousTiles
    curAnomalousTileName = anomalyTileNameWithProjCRS{idxT, 1};
    curRefTileName = anomalyTileNameWithProjCRS_Debug{idxT, 1};

    curProjNameRef = anomalyTileNameWithProjCRS_Debug{idxT, 2} ...
        .ProjectedCRS.Name;

    % Load the repaired tile.
    [lidarDataImg, spatialRef] ...
        = readgeoraster(fullfile(pathToSaveRepairedTiles, ...
        [curAnomalousTileName, '.tif']));
    lidarDataImg(abs( ...
        lidarDataImg(:))>maxAllowedAbsLidarZ) = nan;
    % Essentailly meshgrid matrices.
    [lidarRasterXs, lidarRasterYs] = worldGrid(spatialRef);
    % Column vectors.
    [lidarLats, lidarLons] = projinv( ...
        spatialRef.ProjectedCRS, ...
        lidarRasterXs(:), lidarRasterYs(:));
    % Convert survery feet to meter.
    lidarZs = distdim(lidarDataImg(:), 'ft', 'm');

    % Plot.
    figure('Visible', ~flagGenFigsQuietly);
    plot3k([lidarLons, lidarLats, lidarZs], ...
        'Labels', {'', 'Longitude (degree)', 'Latitude (degree)', ...
        '', 'DSM LiDAR z (m)'}); view(2);
    zlim([0, max(lidarZs)]);
    title({['Repaired Anomalous Tile ', curAnomalousTileName], ...
        ['Based on Tile ', curRefTileName], ...
        ['(Ref Proj: ', char(curProjNameRef), ')']}, ...
        'Interpreter', 'none');
    plot_google_map('MapType', 'hybrid');
    saveas(gcf, fullfile(pathToSaveDebugFigs, ...
        [curAnomalousTileName, '.jpg']));
    saveas(gcf, fullfile(pathToSaveDebugFigs, ...
        [curAnomalousTileName, '.fig']));
    close(gcf);
end

%% Cleanup

close all;
diary off;

% EOF