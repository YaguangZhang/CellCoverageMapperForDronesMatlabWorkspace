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
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
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
        'terrainEles', 'lidarZs');
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

% EOF