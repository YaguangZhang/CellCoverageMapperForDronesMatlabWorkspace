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

%% Test Case with Unknown Projection Name

%   - The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';
dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM_UnknownProjTiles');

% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

addpath(fullfile(pwd, 'lib', 'lidar'));
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
    zlim([0, max(lidarZs)]);

    lidarZsInt = getLiDarZFromXYFct(lidarXs, lidarYs);
    lidarZsInt(isinf(lidarZsInt(:))) = nan;

    figure; plot3k([lidarLons, lidarLats, lidarZsInt]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarZsInt)]);

    figure; plot3k([lidarLons, lidarLats, lidarEles]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);

    lidarElesInt = getEleFromXYFct(lidarXs, lidarYs);
    figure; plot3k([lidarLons, lidarLats, lidarElesInt]); view(2);
    plot_google_map('MapType', 'satellite');
    zlim([0, max(lidarEles)]);
end

% EOF