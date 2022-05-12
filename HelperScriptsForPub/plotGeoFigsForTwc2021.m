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
NUM_OF_PIXELS_FOR_LONGER_SIDE = 100;

%% LiDAR Data Info

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading Indiana LiDAR meta data ...'])

[inBoundaryLatLons, inBoundaryXYs, inBoundaryUtmZone] = loadInBoundary;
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

%% Plots

%% Cleanup

close all;

% EOF