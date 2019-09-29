% SIMULATECHANNELFORINDIANA Generate the coverage and blockage maps for the
% Indiana state.
%
% If new simulation is needed, please delete previous result files to avoid
% progress recovery.
%
% Yaguang Zhang, Purdue, 09/10/2019

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Script Parameters

% The absolute path to the Lidar .las file.
LIDAR_DATA_SET_TO_USE = 'IN';

% The absolute path to the antenna infomation file.
ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'RandomizedCarrierSitesv2.csv');
% The cellular tower height value to use.
DEFAULT_TX_ANT_HEIGHT_IN_M = 50;

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '6_SimulationForIndiana');

%% Simulation Configurations

% We will organize all the necessary configurations into a structure.
%   - A string label to identify this simulation.
simConfigs.CURRENT_SIMULATION_TAG = 'IN';

%   - The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';

%   - The UTM (x, y) polygon representing the area of interest for
%   generating the coverage maps; this is default to the range covered by
%   the availabe LiDAR data set when it is empty.
simConfigs.UTM_X_Y_POLYGON_OF_INTEREST = [];

%   - We will use this number of pixels for the longer side (width/height)
%   of the map; the number of pixels for the other side will be
%   proportional to its length.
simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE = 100;

%   - The guaranteed spacial resolution for terrain profiles; a larger
%   value will decrease the simulation time but small obstacles may get
%   ingored.
simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M = 5;

%   - The guaranteed minimum number of LiDAR z (or possibly elevation)
%   elements in one terrain profile; this will ensure non-empty terrain
%   profiles.
simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE = 10;

%% Preprocessing LiDAR Data

DIR_TO_LIDAR_FILES = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', LIDAR_DATA_SET_TO_USE);

% Convert GPS degrees to UTM coordinates for the specified zone.
utmstruct_speZone = defaultm('utm');
% Remove white space in the zone label.
utmstruct_speZone.zone ...
    = simConfigs.UTM_ZONE(~isspace(simConfigs.UTM_ZONE));
utmstruct_speZone.geoid = wgs84Ellipsoid;
utmstruct_speZone = defaultm(utmstruct_speZone);

deg2utm_speZone = @(lat, lon) mfwdtran(utmstruct_speZone, lat,lon);
utm2deg_speZone = @(x, y) minvtran(utmstruct_speZone, x, y);

% Preprocess .img LiDAR data.
[lidarFileRelDirs, lidarFileXYBoundries, ~] ...
    = preprocessIndianaLidarDataSet(DIR_TO_LIDAR_FILES, ...
    deg2utm_speZone, utm2deg_speZone);
lidarFileAbsDirs = cellfun(@(d) [DIR_TO_LIDAR_FILES, d], ...
    lidarFileRelDirs, 'UniformOutput', false);

%% Load Cellular Tower Information

disp(' ')
disp('    Loading cellular antenna information ...')

% Note: we use "height" to indicate the vertical distance from the ground
% to the antenna; "elevation" to indicate the ground elevation; and
% "altitude" to indicate elevation+height.
cellAntsLatLon = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);
[numAnts, ~] = size(cellAntsLatLon);
cellAntsXYH = nan(numAnts, 3);
[cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
    = deg2utm_speZone(cellAntsLatLon(:,1), cellAntsLatLon(:,2));
cellAntsXYH(:,3) = DEFAULT_TX_ANT_HEIGHT_IN_M;

disp('    Done!')

%% Simulation



% EOF