% SIMULATECHANNELFOREXTENDEDTIPP Generate the coverage and blockage maps
% for area around Tippecanoe County, Indiana. We will consider counties in
% the Wabash Heartland Innovation Network, with an extended area to
% consider cellular towers around.
%
% If new simulation is needed, please delete previous result files to avoid
% progress recovery.
%
% Yaguang Zhang, Purdue, 08/28/2019

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Change PRESET to run the code for different areas.
SUPPORTED_PRESETS = {'ACRE', 'Tipp', 'Tipp_Extended'};
PRESET = 'ACRE';

assert(any(strcmp(SUPPORTED_PRESETS, PRESET)), ...
    ['Unsupported preset "', PRESET, '"!']);

%% Script Parameters

% The absolute path to the Lidar .las file. Currently supporting
% 'Tipp_Extended' (for ten counties in the WHIN are), 'IN' (all indiana)
% and '' (automatically pick the biggest processed set).
LIDAR_DATA_SET_TO_USE = '';

% The absolute path to the antenna infomation file.
ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'RandomizedCarrierSitesv2.csv');
% The cellular tower height value to use.
DEFAULT_TX_ANT_HEIGHT_IN_M = 50;

% The absolute path to save results.
switch PRESET
    case 'ACRE'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_1_SimulationForAcre');
    case 'Tipp'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_2_SimulationForTipp');
    case 'Tipp_Extended'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_3_SimulationForExtendedTipp');
end

%% Simulation Configurations

% We will organize all the necessary configurations into a structure.
%   - A string label to identify this simulation.
switch PRESET
    case 'ACRE'
        simConfigs.CURRENT_SIMULATION_TAG = 'ACRE';
    case 'Tipp'
        simConfigs.CURRENT_SIMULATION_TAG = 'Tipp';
    case 'Tipp_Extended'
        simConfigs.CURRENT_SIMULATION_TAG = 'ExtendedTipp';
end

%   - The UTM (x, y) polygon boundary vertices representing the area of
%   interest for generating the coverage maps; this is default to the range
%   covered by the availabe LiDAR data set when it is empty.
switch PRESET
    case 'ACRE'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.467341, -87.015762; ...
            40.501484, -86.979905]);
    case 'Tipp'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.216047, -87.093700; ...
            40.562743, -86.695913]);
    case 'Tipp_Extended'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = [];
end

%   - Carrier frequency and wavelength.
simConfigs.CARRIER_FREQUENCY_IN_MHZ = 1900;
simConfigs.CARRIER_WAVELENGTH_IN_M ...
    = physconst('LightSpeed')/simConfigs.CARRIER_FREQUENCY_IN_MHZ/1e6;

%   - The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';

%   - We will use this number of pixels for the longer side (width/height)
%   of the map; the number of pixels for the other side will be
%   proportional to its length.
simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE = 100;

%   - The guaranteed spacial resolution for terrain profiles; a larger
%   value will decrease the simulation time but small obstacles may get
%   ingored.
simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M = 50;
%   - Similarly, the guaranteed spacial resolution for LiDAR profiles.
simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M = 50;

%   - The guaranteed minimum number of LiDAR z (or possibly elevation)
%   elements in one terrain profile; this will ensure non-empty terrain
%   profiles.
simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE = 10;

%   - Rx heights for different algorithms to inspect.
simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = [1.5; 10; 30; 50; 70; 90; 100];

%   - The maximum radius that a cellular tower can cover; we will use this
%   to limit the area to consider for each cellular tower in the
%   simulation; by default (when it is not a positive scalar) it will be
%   set to be the distance that one can see at the top of the highest cell
%   tower to the highest RX above the horizon (i.e. without blockage of the
%   earth): D_BL_IN_KM ~= 3.57(sqrt(h_TX_IN_M)+sqrt(h_RX_IN_M)).
simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M = nan;

%   - The clearance ratio of the first Fresnel zone for a LoS path: at
%   least this ratio of the first Fresnel zone needs to be clear for a path
%   to be considered as "line of sight" (LoS); we expect the value to be
%   larger or equal to 50%.
simConfigs.LOS_FIRST_FRES_CLEAR_RATIO = 0.6;

%   - For the NTIA eHata library: the quantile percent not exceeded of the
%   signal. Limits: 0 < reliability < 1.
simConfigs.NTIA_EHATA_RELIABILITY = 0.95;
%   - For the NTIA eHata library: the NLCD environment code.
simConfigs.NTIA_EHATA_ENVIRO_CODE = 82; % Cultivated Crops.

%   - Functions to convert GPS degrees (lat, lon) from/to UTM (x, y). We
%   will polulate them later.

% Convert GPS degrees to UTM coordinates for the specified zone.
utmstruct_speZone = defaultm('utm');
% Remove white space in the zone label.
utmstruct_speZone.zone ...
    = simConfigs.UTM_ZONE(~isspace(simConfigs.UTM_ZONE));
utmstruct_speZone.geoid = wgs84Ellipsoid;
utmstruct_speZone = defaultm(utmstruct_speZone);

deg2utm_speZone = @(lat, lon) mfwdtran(utmstruct_speZone, lat,lon);
utm2deg_speZone = @(x, y) minvtran(utmstruct_speZone, x, y);

% Store these functions in simConfigs.
simConfigs.deg2utm_speZone = deg2utm_speZone;
simConfigs.utm2deg_speZone = utm2deg_speZone;

%% Preprocessing LiDAR Data

[LIDAR_DATA_SET_TO_USE] = verifyLidarDataSetToUse( ...
    LIDAR_DATA_SET_TO_USE, ABS_PATH_TO_SHARED_FOLDER);

DIR_TO_LIDAR_FILES = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', LIDAR_DATA_SET_TO_USE);

% Preprocess .img LiDAR data. To make Matlab R2019b work, we need to remove
% preprocessIndianaLidarDataSet from path after things are done.
simConfigs.LIDAR_DATA_SET_TO_USE = LIDAR_DATA_SET_TO_USE;
addpath(fullfile(pwd, 'lib', 'lidar'));
[lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessIndianaLidarDataSet(DIR_TO_LIDAR_FILES, ...
    deg2utm_speZone, utm2deg_speZone);
rmpath(fullfile(pwd, 'lib', 'lidar'));
lidarFileAbsDirs = cellfun(@(d) ...
    [DIR_TO_LIDAR_FILES, strrep(d, '\', filesep)], ...
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

% We need Python for concurrent HTTP requests to get elevation data from
% USGS faster.
if isempty(pyversion)
    pyversion(ABS_PATH_TO_PYTHON);
end
% Generate coverage and blockage maps.
[simState, simConfigs] ...
    = analyzeCoverageViaChannelSimulation(pathToSaveResults, ...
    lidarFileAbsDirs, lidarFileXYCoveragePolyshapes, ...
    cellAntsXYH, simConfigs);

% EOF