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

%% Load Boundary

absDirToAcreKmz = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', 'ACRE', 'AcreExactBoundaryRaw', 'ACRE.kmz');

% Boundary of ACRE.
[UTM_X_Y_BOUNDARY_ACRE, UTM_ZONE_ACRE, kmzStruct, ...
    utmXYPolygons, lonLatPolygons] ...
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
gridResolutionInM = 1.5;

gridUtmXYPts = buildSimGrid(UTM_X_Y_BOUNDARY_ACRE, ...
    gridResolutionInM, true);

% EOF