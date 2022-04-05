% DEBUGPREPROCESSLIDARDSM For investigating tiles with unknown projection.
%
% Please make sure to run this before doing the simulation to make sure the
% environment is set up correctly.
%
% Yaguang Zhang, Purdue, 03/25/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', 'PreprocessLidarFctDebugging');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

% The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';
% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

addpath(fullfile(pwd, 'lib', 'lidar'));

% Save output into a log file.
dirToDiary = fullfile(pathToSaveResults, 'diary.log');
diary(dirToDiary);

%% Tests for Python Environment and Downloader

py_addpath(fullfile(pwd, 'lib', 'python'));

disp('Fetching by lib ''Matlab'':')
tic;
elesTestMatlab = queryElevationPointsFromUsgs(40:50, -80:-70, 'Matlab');
toc;

disp('Fetching by lib ''Python'':')
tic;
elesTestPy = queryElevationPointsFromUsgs(40:50, -80:-70);
toc;

disp('Fetching by lib ''Python'' in Chunks:')
tic;
elesTestPyChunk = queryElevationPointsFromUsgsInChunks(40:50, -80:-70);
toc;

assert(all(isnan(elesTestMatlab)==isnan(elesTestPy)) ...
    && all(elesTestMatlab(~isnan(elesTestMatlab)) ...
    ==elesTestPy(~isnan(elesTestPy))) && ...
    ...
    all(isnan(elesTestPyChunk)==isnan(elesTestPy)) ...
    && all(elesTestPyChunk(~isnan(elesTestPyChunk)) ...
    ==elesTestPy(~isnan(elesTestPy))), ...
    'Results fetched from different libraries do not agree!');

%% Debug Case: Check the Spatial Reference for IN tiles on Frankie.

absPathToSaveProjCrs = fullfile(pwd, 'lib', 'lidar', ...
    'anomalyTileNameWithProjCRS.mat');
if ~exist(absPathToSaveProjCrs, 'file')
    diary(fullfile(pathToSaveResults, 'diary.log'));

    dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar_2019', 'IN', 'DSM');
    tifFileHandles = rdir(fullfile( ...
        dirToLidarFiles, '**', '*.tif'), '', dirToLidarFiles);
    lidarFileRelDirs = {tifFileHandles(:).name}';

    % Extract number IDs and sort the file name list accordingly.
    regexpPat = ['\/(\d+)_dsm.tif|MC_(\d+)_dsm.tif', ...
        '|in\d+_(\d+)_\d+_dsm.tif|IN\d+_(\d+)_\d+_dsm.tif'];
    tileIdTokens = cellfun(@(relDir) regexp(relDir, regexpPat, 'tokens'), ...
        lidarFileRelDirs);
    % Note that some files will have same IDs.
    [allTileIdNums, indicesNewOrder] = sort(cellfun( ...
        @(token) str2double(token{1}), tileIdTokens));

    allTileAbsDirs = cellfun(@(relDir) fullfile(dirToLidarFiles, relDir), ...
        lidarFileRelDirs, 'UniformOutput', false);
    allTileAbsDirs = allTileAbsDirs(indicesNewOrder);

    anomalyTileAbsDirs = {};
    anomalyTileIdNums = [];
    for idxLidarF = 1:length(allTileAbsDirs)
        [~, R] = readgeoraster(allTileAbsDirs{idxLidarF});
        if ~isa(R.ProjectedCRS, 'projcrs')
            warning('Empty ProjectedCRS!')
            disp('    Available R:')
            disp(R)
            disp(' ')
            disp('    File path:')
            anomalyTileAbsDirs{end+1,1} ...
                = allTileAbsDirs{idxLidarF}; %#ok<SAGROW>
            anomalyTileIdNums(end+1,1) ...
                = allTileIdNums(idxLidarF); %#ok<SAGROW>
            disp(anomalyTileAbsDirs{end})
            disp(' ')
        end
    end

    % Fetch the nearest tile (based on filename) that has a valid spatial
    % reference.
    [validTileAbsDirs, indicesValidTiles] = setdiff( ...
        allTileAbsDirs, anomalyTileAbsDirs);
    validTileIdNums = allTileIdNums(indicesValidTiles);

    % For plotting LiDAR z. Any LiDAR z value too big or too small will be
    % discarded (set to NaN).
    maxAllowedAbsLidarZ = 10^38;

    numOfAnomalyTiles = length(anomalyTileAbsDirs);
    anomalyTileNameWithProjCRS = cell(numOfAnomalyTiles, 2);

    for idxAnomalyT = 1:numOfAnomalyTiles
        [curLidarImg, curR] = readgeoraster(anomalyTileAbsDirs{idxAnomalyT});
        curTileIdNum = anomalyTileIdNums(idxAnomalyT);

        [idDiff, idxValidTile] = min(abs(validTileIdNums-curTileIdNum));
        [~, tentativeR] = readgeoraster(validTileAbsDirs{idxValidTile});

        [lidarRasterXs, lidarRasterYs] = worldGrid(curR);
        [lidarLats, lidarLons] = projinv( ...
            tentativeR.ProjectedCRS, ...
            lidarRasterXs(:), lidarRasterYs(:));

        curLidarImg(abs(curLidarImg(:))>maxAllowedAbsLidarZ) = nan;
        % Do not plot anything for water.
        curLidarImg(curLidarImg(:)==0) = nan;
        lidarZs = distdim(curLidarImg(:), 'ft', 'm');

        [~, anomalyTileName] = fileparts(anomalyTileAbsDirs{idxAnomalyT});
        [~, tentativeRefTileName] = fileparts(validTileAbsDirs{idxValidTile});

        figure; plot3k([lidarLons, lidarLats, lidarZs]); view(2);
        plot_google_map('MapType', 'satellite'); zlim([0, max(lidarZs)]);
        title(['Anomaly tile ', anomalyTileName, ...
            ' projected based on ', tentativeRefTileName], ...
            'Interpreter', 'none');
        saveas(gcf, fullfile(pathToSaveResults, [anomalyTileName, '.jpg']));
        close(gcf);

        anomalyTileNameWithProjCRS{idxAnomalyT, 1} = anomalyTileName;
        anomalyTileNameWithProjCRS{idxAnomalyT, 2} = tentativeR.ProjectedCRS;
    end

    save(absPathToSaveProjCrs, 'anomalyTileNameWithProjCRS');
    diary off;
end

%% Test Case: Tiles with Unknown Projection Names

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_UnknownProjTiles');
[lidarFileRelDirs1, ~, lidarFileLonLatCoveragePolyshapes1] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone, true);

for idxF = 1:length(lidarFileRelDirs1)
    [curDir, curFileName] = fileparts(lidarFileRelDirs1{idxF});
    load(fullfile(dirToLidarFiles, 'MatlabCache', [curFileName,'.mat']));
    lidarZs(isinf(lidarZs(:))) = nan;

    curLonLats = lidarFileLonLatCoveragePolyshapes1{idxF}.Vertices;
    curLonLats(end+1,:) = curLonLats(1,:); %#ok<SAGROW>

    figure; plot3k([lidarLons, lidarLats, lidarZs]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarZs), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZs)]);

    figure; plot3k([lidarLons, lidarLats, lidarEles]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarEles), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);
end

%% Test Case: ACRE

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_InspectGroundEle_ACRE_TestNewPrep');
[ lidarFileRelDirs2, ~, lidarFileLonLatCoveragePolyshapes2 ] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone, true);

for idxF = 1:length(lidarFileRelDirs2)
    [curDir, curFileName] = fileparts(lidarFileRelDirs2{idxF});
    load(fullfile(dirToLidarFiles, 'MatlabCache', [curFileName,'.mat']));
    lidarZs(isinf(lidarZs(:))) = nan;

    curLonLats = lidarFileLonLatCoveragePolyshapes2{idxF}.Vertices;
    curLonLats(end+1,:) = curLonLats(1,:); %#ok<SAGROW>

    figure; plot3k([lidarLons, lidarLats, lidarZs]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarZs), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZs)]); title('Orignal LiDAR z');

    lidarZsInt = getLiDarZFromXYFct(lidarXs, lidarYs);
    lidarZsInt(isinf(lidarZsInt(:))) = nan;

    figure; plot3k([lidarLons, lidarLats, lidarZsInt]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarZsInt), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZsInt)]); title('LiDAR z from getLiDarZFromXYFct');

    figure; plot3k([lidarLons, lidarLats, lidarEles]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarEles), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]); title('Orignal Ground Elevation');

    lidarElesInt = getEleFromXYFct(lidarXs, lidarYs);
    figure; plot3k([lidarLons, lidarLats, lidarElesInt]); view(2); hold on;
    plot3(curLonLats(:,1), curLonLats(:,2), ...
        ones(length(curLonLats(:,1)), 1).*max(lidarElesInt), 'r-.');
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);
    title('Ground Elevation from getEleFromXYFct');
end

%% Debug Case: All IN Tiles on Frankie.
% Note 20220328: We used MATLAB R2021a on Frankie.

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM');
[lidarFileRelDirs2, lidarFileXYCoveragePolyshapes2, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);

%% Debug Case: Selected Tiles which Caused Error

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_ErrInFetchEle');
[lidarFileRelDirs3, lidarFileXYCoveragePolyshapes3, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);

%% Cleanup

diary off;

% EOF