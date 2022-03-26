% DEBUGPREPROCESSLIDARDSM For investigating tiles with unknown projection.
%
% Yaguang Zhang, Purdue, 03/25/2022

clear; clc; close all; dbstop if error;

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

% The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';
% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

addpath(fullfile(pwd, 'lib', 'lidar'));

%% Test Case: Tiles with Unknown Projection Names

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_UnknownProjTiles');
[lidarFileRelDirs1, lidarFileXYCoveragePolyshapes1, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);

for idxF = 1:length(lidarFileRelDirs1)
    [curDir, curFileName] = fileparts(lidarFileRelDirs1{idxF});
    load(fullfile(dirToLidarFiles, 'MatlabCache', [curFileName,'.mat']));
    lidarZs(isinf(lidarZs(:))) = nan; %#ok<SAGROW>

    figure; plot3k([lidarLons, lidarLats, lidarZs]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZs)]);

    figure; plot3k([lidarLons, lidarLats, lidarEles]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);
end

%% Test Case: ACRE

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_InspectGroundEle_ACRE_TestNewPrep');
[lidarFileRelDirs2, lidarFileXYCoveragePolyshapes2, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);

for idxF = 1:length(lidarFileRelDirs2)
    [curDir, curFileName] = fileparts(lidarFileRelDirs2{idxF});
    load(fullfile(dirToLidarFiles, 'MatlabCache', [curFileName,'.mat']));
    lidarZs(isinf(lidarZs(:))) = nan;

    figure; plot3k([lidarLons, lidarLats, lidarZs]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZs)]); title('Orignal LiDAR z');

    lidarZsInt = getLiDarZFromXYFct(lidarXs, lidarYs);
    lidarZsInt(isinf(lidarZsInt(:))) = nan;

    figure; plot3k([lidarLons, lidarLats, lidarZsInt]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZsInt)]); title('LiDAR z from getLiDarZFromXYFct');

    figure; plot3k([lidarLons, lidarLats, lidarEles]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]); title('Orignal Ground Elevation');

    lidarElesInt = getEleFromXYFct(lidarXs, lidarYs);
    figure; plot3k([lidarLons, lidarLats, lidarElesInt]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);
    title('Ground Elevation from getEleFromXYFct');
end

%% Debug Case: Check the Spatial Reference for IN tiles on Frankie.

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM');
tifFileHandles = rdir(fullfile( ...
    dirToLidarFiles, '**', '*.tif'), '', dirToLidarFiles);
lidarFileRelDirs = {tifFileHandles(:).name}';

anomalyTileAbsDirs = {};
for idxLidarF = 1:length(lidarFileRelDirs)
    [~, R] =readgeoraster(fullfile(dirToLidarFiles, ...
        lidarFileRelDirs{idxLidarF}));
    if ~isa(R.ProjectedCRS, 'projcrs')
        warning('Empty ProjectedCRS!')
        anomalyTileAbsDirs{end+1} = fullfile( ...
            dirToLidarFiles, lidarFileRelDirs{idxLidarF}); %#ok<SAGROW>
        disp(anomalyTileAbsDirs{end})
        disp(' ')
    end
end

%% Debug Case: All IN tiles on Frankie.

dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM');
[lidarFileRelDirs2, lidarFileXYCoveragePolyshapes2, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);

% EOF