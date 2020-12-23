% SIMULATECHANNELFOREXTENDEDTIPP Generate the coverage and blockage maps
% for area around Tippecanoe County, Indiana. We will consider counties in
% the Wabash Heartland Innovation Network, with an extended area to
% consider cellular towers around.
%
% This simulator has configuration presets available for carrying out
% simulation with different area of interest. Please adjust the value of
% the viable PRESET below to choose the target simulation.
%
% If new simulation is needed, please delete previous result files to avoid
% progress recovery.
%
% Matlab R2019b or newer is required to utilize the custom Python module
% for downloading elevation data from USGS concurrently.
%
% Yaguang Zhang, Purdue, 08/28/2019

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Change PRESET to run the code for different areas.
%   - 'ACRE'
%     Purdue ACRE research farm with manually labeled square boundary
%     (larger than that for 'ACRE_EXACT').
%   - 'ACRE_EXACT'
%     Purdue ACRE research farm with the exact boundary imported from the
%     offical kmz boundary file.
%   - 'Tipp'
%     Tippecanoe county.
%   - 'ExtendedTipp'
%     WHIN area.
%   - 'IN'
%     Indiana state.
%   - 'Aerostat'
%     Simulation for one areostat as the cell tower case.
SUPPORTED_PRESETS = {'ACRE', 'ACRE_EXACT', ...
    'Tipp', 'ExtendedTipp', 'IN', ...
    'Aerostat'};
PRESET = 'IN';

assert(any(strcmp(SUPPORTED_PRESETS, PRESET)), ...
    ['Unsupported preset "', PRESET, '"!']);

%% Script Parameters

% The absolute path to the Lidar .las file. Currently supporting
% 'Tipp_Extended' (for ten counties in the WHIN are), 'IN' (all indiana)
% and '' (automatically pick the biggest processed set).
LIDAR_DATA_SET_TO_USE = '';

switch PRESET
    case 'Aerostat'
        % The absolute path to the antenna infomation file.
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'Aerostat', 'AerostatLatLonAlt.csv');
    otherwise
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'RandomizedCarrierSitesv2.csv');
        % The cellular tower height value to use.
        DEFAULT_TX_ANT_HEIGHT_IN_M = 50;
end

% The absolute path to save results.
switch PRESET
    case 'ACRE'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_1_SimulationForAcre');
    case 'Tipp'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_2_SimulationForTipp');
    case 'ExtendedTipp'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_3_SimulationForExtendedTipp');
    case 'IN'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_4_SimulationForIn');
    case 'ACRE_EXACT'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_5_SimulationForAcreExact');
    case 'Aerostat'
        pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', '5_6_SimulationForAerostat');
end

%% Simulation Configurations

% We will organize all the necessary configurations into a structure called
% simConfigs.
%   - A string label to identify this simulation.
simConfigs.CURRENT_SIMULATION_TAG = PRESET;

%   - The UTM (x, y) polygon boundary vertices representing the area of
%   interest for generating the coverage maps; for presets ExtendedTipp and
%   IN, this is default to the range covered by the availabe LiDAR data set
%   for the corresponding area of interest when it is empty.
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
    case 'ExtendedTipp'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = [];
    case 'IN'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = [];
    case 'ACRE_EXACT'
        % Read in the ACRE boundary.
        ABS_DIR_TO_ACRE_EXACT_BOUNDARY ...
            = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar', 'ACRE', 'AcreExactBoundaryRaw', 'ACRE.kmz');
        [simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ~] ...
            = extractBoundaryFromKmzFile(ABS_DIR_TO_ACRE_EXACT_BOUNDARY);
    case 'Aerostat'
        % We will start with the IN boundary and intersect it with a square
        % centered at the aerostat.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = [];
        squareOfInterestSideLengthInM = distdim(100, 'miles', 'm');
end

%   - Carrier frequency and wavelength. Typical values:
%       - 1900 MHz for cellular LTE
%        - 28000 MHz (28 GHz) for millimeter wave
%       - 3500 MHz (3.5 GHz) and 915 MHz for areostat applications
simConfigs.CARRIER_FREQUENCY_IN_MHZ = 28000;
simConfigs.CARRIER_WAVELENGTH_IN_M ...
    = physconst('LightSpeed')/simConfigs.CARRIER_FREQUENCY_IN_MHZ/1e6;

%   - The zone label to use in the UTM (x, y) system.
simConfigs.UTM_ZONE = '16 T';

%   - We will use this number of pixels for the longer side (width/height)
%   of the map; the number of pixels for the other side will be
%   proportional to its length.
simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE = 100; %256; % 100;

%   - The guaranteed spacial resolution for terrain profiles; a larger
%   value will decrease the simulation time but small obstacles may get
%   ingored.
simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M = 7.5; % 50;
%   - Similarly, the guaranteed spacial resolution for LiDAR profiles.
simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M = 1.5; % 50;

%   - The guaranteed minimum number of LiDAR z (or possibly elevation)
%   elements in one terrain profile; this will ensure non-empty terrain
%   profiles.
simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE = 10;

%   - Rx heights for the simulator to inspect. Note that according to FAA,
%   the maximum allowable altitude is 400 feet (~122m) above the ground.
%   Typical values:
%       - [1.5; (10:10:120)'; 125] for cellular and millimeter wave
%         inspection
%        - 1.5 for 3500 MHz (3.5 GHz) areostat application
%       - 0.1 for 915 MHz areostat application
simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
    = 1.5;

%   - The maximum radius that a cellular tower can cover; we will use this
%   to find the cellular towers that are effective and limit the area to
%   consider for each effective cellular tower during the simulation; by
%   default (when it is not a positive scalar) we will compute a value
%   based on the radio horizon.
%
%   Please see the field 'MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY' for more
%   information.
simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M = nan;

%   - Optional field (with a string value of 'LowestAntennas' or
%   'HighestAntennas'; default to 'LowestAntennas') to control how
%   simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M is computed automatically, if
%   the value of that radius is not valid. The radius will be computed
%   using the lowest (heightest) cell tower and the lowest (heightest)
%   receiver if this flag is set to 'LowestAntennas' ('HighestAntennas').
%   One can choose 'HighestAntennas' to consider more towers which may be
%   able to see the heightest receiver at some spot in the area of interest
%   for accuracy, or 'LowestAntennas' to speed up the simulation.
% simConfigs.MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY = 'HighestAntennas';

%   - The function to calculate, according to the TX and RX heights, the
%   maximum cellular tower coverage radius. For example, it can be set to
%   be the distance that one can see at the top of the highest/lowest cell
%   tower to the highest/lowest RX above the horizon (i.e. without blockage
%   of the earth):
%
%       D_BL_IN_KM ~= 3.57(sqrt(h_TX_IN_M)+sqrt(h_RX_IN_M)).
%
%   Here, we use the radio horizon instead (3.57 => 4.12). References:
%
%       Haslett, Christopher. Essentials of Radio Wave Propagation,
%       Cambridge University Press, 2007. ProQuest Ebook Central,
%       https://ebookcentral.proquest.com/lib/purdue/detail.action?docID=803105.
%
%       WikiPedia articles:
%               https://en.wikipedia.org/wiki/Horizon,
%           and
%               https://en.wikipedia.org/wiki/Line-of-sight_propagation.
%
simConfigs.getCellCoverageRadiusInM = @(TxHeightInM, RxHeightInM) ...
    4.12*(sqrt(TxHeightInM)+sqrt(RxHeightInM))*1000;

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

%   - For adjusting the feedback frequency in parfor workers.
simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT = 0.2;

%   - For plotting, we need an expected value for the maximum allowed path
%   loss. With 20 MHz bandwidth, 9 dB RX noise figure, and 100 W TX power,
%   we have a maximum path loss of ~142 dB.
simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB = [0, 150];

%   - For figures, we need to decide whether to resize them for
%   publication.
simConfigs.RESIZE_FIG_FOR_PUBLICATION = false;

% Turn the diary logging function on.
dirToSaveDiary = fullfile(pathToSaveResults, 'diary.txt');
if ~exist(dirToSaveDiary, 'file')
    if ~exist(pathToSaveResults, 'dir')
        mkdir(pathToSaveResults)
    end
    fclose(fopen(dirToSaveDiary, 'w'));
end
diary(dirToSaveDiary);

%% Preprocessing LiDAR Data

% Make sure the chosen LiDAR dataset to use in the simulation is indeed
% available, and if it is not specified (LIDAR_DATA_SET_TO_USE = ''),
% default to the bigger LiDAR dataset that has been preprocessed before for
% better coverage.
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
cellAntsLatLonAlt = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);
numAnts = size(cellAntsLatLonAlt, 1);
cellAntsXYH = nan(numAnts, 3);
[cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
    = deg2utm_speZone(cellAntsLatLonAlt(:,1), cellAntsLatLonAlt(:,2));
if size(cellAntsLatLonAlt, 2)>2
    cellAntsXYH(:,3) = cellAntsLatLonAlt(:,3);
end
boolsCellAntAltUnknown = isnan(cellAntsXYH(:,3));
if sum(boolsCellAntAltUnknown)>0
    cellAntsXYH(boolsCellAntAltUnknown,3) = DEFAULT_TX_ANT_HEIGHT_IN_M;
end

disp('    Done!')

%% Validate Area of Interest
% Generate the area of interest according to the corresponding LiDAR
% dataset for that area (indicated by PRESET), regardless of what LiDAR
% dataset to be used in the simulator for terrain/LiDAR profile generation.
% This will ensure an accurate polygon representation of the area of
% interest, while retaining the option of using other LiDAR data for the
% simulation, e.g. a big comprehensive local LiDAR data set on a powerful
% cluster machine to save time used in data fetching.

if isempty(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST)
    % So far, we have presets 'ExtendedTipp', 'IN', and 'Aerostat' that
    % need automatically generated areas of interset.
    assert(any(strcmp({'ExtendedTipp', 'IN', 'Aerostat'}, PRESET)), ...
        ['No raw Indiana LiDAR data available for ', PRESET, ...
        ' to generate the area of interest polygon!'])

    % Take care of the cases where the LiDAR dataset folder is different
    % from the simulation tag, e.g., the extended Tippecanoe area case.
    switch PRESET
        case 'ExtendedTipp'
            areaOfInterestLidarDataset = 'Tipp_Extended';
        case 'Aerostat'
            areaOfInterestLidarDataset = 'IN';
        otherwise
            areaOfInterestLidarDataset = PRESET;
    end

    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
        = generateUtmXyBoundaryOfLidarDataset( ...
        fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar', areaOfInterestLidarDataset), simConfigs);
end

% For the areostat case, we still need to refine the area of interest.
switch PRESET
    case 'Aerostat'
        % The area of interest is the intersect of the LiDAR coverage area
        % and the square centered at the aerostat.
        areaOfInterestPoly = intersect( ...
            polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
            polyshape(constructUtmRectanglePolyMatFromCenter( ...
            cellAntsLatLonAlt(1, 1:2), ...
            [1,1].*squareOfInterestSideLengthInM)));
        [areaOfInterestXs, areaOfInterestYs] ...
            = boundary(areaOfInterestPoly);
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = [areaOfInterestXs, areaOfInterestYs];
end

%% Simulation

disp(' ')
disp('    Conducting simulation ...')
disp(' ')

pathToSaveSimResults = fullfile(pathToSaveResults, 'simResults.mat');
if exist(pathToSaveSimResults, 'file')
    disp('    Loading history results ...')
    load(pathToSaveSimResults)
else
    % Computing coverage and blockage maps.
    [simState, simConfigs] ...
        = analyzeCoverageViaChannelSimulation(pathToSaveResults, ...
        lidarFileAbsDirs, lidarFileXYCoveragePolyshapes, ...
        cellAntsXYH, simConfigs);
    save(pathToSaveSimResults, 'simConfigs', 'simState', '-v7.3');
end
disp('    Done!')

% Plots.
plotCoverageSimulationResults(pathToSaveResults, simState, simConfigs);

% EOF