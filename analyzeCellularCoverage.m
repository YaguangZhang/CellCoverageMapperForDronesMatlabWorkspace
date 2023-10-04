%ANALYZECELLULARCOVERAGE Cellular coverage simulator for blockage and path
%loss maps.
%
% This script is based on:
%   5_SimulationForExtendedTipp/simulateChannelForExtendedTipp.m
%
% Notes:
%   - Please make sure to run the testing script debugPreprocessLidarDsm.m
%     under 1_CoverageMapsForAcre/ before doing the simulation to make sure
%     the environment is set up correctly.
%   - This simulator has configuration presets available for carrying out
%     simulation with different areas of interest. Please adjust the value
%     of the viable PRESET below to choose the desired simulation.
%   - If new simulation is needed, please delete previous result files to
%     avoid progress recovery.
%   - Matlab R2021a & R2021b have been used for developing and
%     testing this simulator. Release R2020b and beyond are recommended for
%     running the simulator.
%
% Yaguang Zhang, Purdue, 04/28/2022

clearvars -except PRESETS CARRIER_FREQUENCIES_IN_MHZ ...
    NUMS_OF_PIXELS_FOR_LONGER_SIDE ...
    PRESET CARRIER_FREQUENCY_IN_MHZ ...
    NUM_OF_PIXELS_FOR_LONGER_SIDE ...
    pathToSaveSimManDiary idxFre idxT idxN...
    CustomSimMetas idxSim CustomSimMeta pathToSaveCustomSimMetas ...
    simConfigFieldsToOverride;
clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Change PRESET to run the code for different areas, e.g. 'Tipp' for
% Tippecanoe county, 'ExtendedTipp' for the WHIN area, and 'IN' for Indiana
% state. Please refer to the Simulation Configurations section for a
% complete list of supported presets.
if ~exist('PRESET', 'var')
    % For testing:
    %   'ACRE_EXACT'.
    % For the journal paper:
    %   {'Tipp', 'ShrinkedIN'}.
    % For LoRaWAN on ACRE:
    %   {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R',
    %    'ACRE_LORA_HALF_MILE_R', 'ACRE_LORA_20KM_R', 'ACRE_LORAWAN_2021'}.
    % For WHIN weather stations (ExtendedTipp as area of interest with
    % weather station locs as the grid):
    %   {'WHIN_WEATHER_STATIONS'}
    % For LoRaWAN built by WHIN (ExtendedTipp as area of interest with full
    % grid user locs):
    %   {'WHIN_LORAWAN'}
    % For LoRaWAN gateways installed on the mobile trailer:
    %   {'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'}
    % Note:
    %   'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY' can be carried out by
    %   runMultipleCellCovSims.m or after idxT is manually set.
    % For custom simulations with parameters set by the variable
    % CustomSimMeta (required fields include CARRIER_FREQUENCY_IN_MHZ,
    % ABS_PATH_TO_CELL_ANTENNAS_CSV, mapGridLatLonPts,
    % UTM_X_Y_BOUNDARY_OF_INTEREST, and RX_ANT_HEIGHTS_TO_INSPECT_IN_M;
    % optional fields include FLAG_EVAL_LOSS_THROUGH_VEG and
    % pathToPostProcessingResultsFolder):
    %   {'CustomSim'}
    % For University of Colorado Boulder campus simulation for a manually
    % drawn square area.
    %   {'UniOfCoBoulderCampus'}
    % For the small grid covering part of the Norlin Quadrangle at the
    % University of Colorado Boulder. The number '#' should be 1 to 4,
    % corresponding to 4 different TX locations on the 4 sides of the grid
    % area, respectively.
    %   {'UniOfCoBoulderCampus_NorlinQuad_#'}
    % For the ACRE/Happy Hollow Park/Lindberg Village scenarios specified
    % by four/two/four tiles from Byunghyun and manually selected towers
    % based on https://www.antennasearch.com/.
    %   {'ACRE_ML', 'HappyHollowPark_ML', 'LindbergVillage_ML'}
    PRESET = 'ACRE_ML';
end

% Suppress selected warnings to reduce messages.
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'export_fig:transparency');
warning('off', 'MATLAB:loadlibrary:parsewarnings');

%% Script Parameters

% The LiDAR data set to use. Currently we only suppor the 2019 Indiana
% state-wide digital surface model (DSM) data from:
%       https://lidar.jinha.org/
% Set this to "IN_DSM_2019"/"IN_DSM_2019_DEMO" for the complete/a demo data
% set. Please refer to the Preprocessing LiDAR Data section for a complete
% list of supported LiDAR data sets.
if strcmpi('UniOfCoBoulderCampus', PRESET(1:min(20,length(PRESET))))
    LIDAR_DATA_SET_TO_USE = 'CO_DRCOG_2020';
else
    LIDAR_DATA_SET_TO_USE = 'IN_DSM_2019';
end

switch PRESET
    case {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
            'ACRE_LORA_HALF_MILE_R', 'ACRE_LORA_20KM_R', ...
            'ACRE_LORAWAN_2021'}
        % We only have one tower for LoRaWAN on ACRE.
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'PurdueAcreLoraWanTowers.csv');
    case {'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'}
        % The WHIN gateway locations.
        pathToWhinGatewayLocs = parseWhinGatewayInfo(...
            ABS_PATH_TO_SHARED_FOLDER);
        ABS_PATH_TO_CELL_ANTENNAS_CSV = pathToWhinGatewayLocs;
    case {'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'}
        % The location information for the three WHIN gateways installed on
        % the mobile trailer.
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'PurdueAcreLoraWanTowersOnTrailer.csv');
    case {'CustomSim'}
        ABS_PATH_TO_CELL_ANTENNAS_CSV ...
            = CustomSimMeta.ABS_PATH_TO_CELL_ANTENNAS_CSV;
    case {'UniOfCoBoulderCampus'}
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'UniOfCoBoulderTower.csv');
    case {'UniOfCoBoulderCampus_NorlinQuad_1', ...
            'UniOfCoBoulderCampus_NorlinQuad_2', ...
            'UniOfCoBoulderCampus_NorlinQuad_3', ...
            'UniOfCoBoulderCampus_NorlinQuad_4'}
        txSiteNumLabel = PRESET(end);
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'ColoradoSims', '20221202_NorlinQuad', ...
            ['TX_', txSiteNumLabel, '_location.csv']);
    case {'ACRE_ML'}
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', ...
            'PurdueAcreCellTowers_AntennaSearch_20230922.csv');
        % 'PurdueAcreCellTowers_AntennaSearch_20230823.csv'
    case {'HappyHollowPark_ML'}
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', ...
            'HappyHollowParkCellTowers_AntennaSearch_20230911.csv');
    case {'LindbergVillage_ML'}
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', ...
            'LindbergVillageCellTowers_AntennaSearch_20230911.csv');
    otherwise
        % Default to the NTIA+HIFLD cell tower locations.
        ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs', ...
            'NtiaLayoutMergedWithHifldCellTs_Threshold_1000m_LatLonH.csv');
end

% The cellular tower height value to use. According to the HIFLD Land
% Mobile Commercial Transmission Towers data set, the mean tower height is
% 63.5 m, while the median tower height is 47.2 m.
DEFAULT_TX_ANT_HEIGHT_IN_M = 50;
% The typical pedestrian height.
PEDESTRIAN_HEIGHT_IN_M = 1.5;

%% Simulation Configurations
% We will organize all the necessary configurations into a structure called
% simConfigs.
%   - A string label to identify this simulation.
% Note that if simConfigFieldsToOverride is defined, the fields in it will
% override the corresponding ones in simConfigs. However, one needs to make
% sure the overridden fields do not affect/depend on other simConfigs
% fields. Otherwise, the simulation may be configured inconsistently.

% ------------------------------------
%  Start of Simulation Configurations
% ------------------------------------
simConfigs.CURRENT_SIMULATION_TAG = PRESET;

%   - The zone label to use in the UTM (x, y) system.
if strcmpi('UniOfCoBoulderCampus', PRESET(1:min(20,length(PRESET))))
    simConfigs.UTM_ZONE = '13 T';
else
    simConfigs.UTM_ZONE = '16 T';
end

% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

%   - By default, do not consider the vegetation simulation (which requires
%   blockage dist results based on ground elevation).
if any(strcmpi({'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
        'ACRE_LORA_HALF_MILE_R', 'ACRE_LORA_20KM_R', ...
        'ACRE_LORAWAN_2021', ...
        'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY', ...
        'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'}, PRESET)) ...
        || strcmpi('UniOfCoBoulderCampus', ...
        PRESET(1:min(20,length(PRESET))))
    simConfigs.FLAG_EVAL_LOSS_THROUGH_VEG = true;
else
    simConfigs.FLAG_EVAL_LOSS_THROUGH_VEG = false;
end

if strcmpi(PRESET, 'CustomSim')
    if isfield(CustomSimMeta, 'FLAG_EVAL_LOSS_THROUGH_VEG')
        simConfigs.FLAG_EVAL_LOSS_THROUGH_VEG ...
            = CustomSimMeta.FLAG_EVAL_LOSS_THROUGH_VEG;
    end
end

%   - The UTM (x, y) polygon boundary vertices representing the area of
%   interest for generating the coverage maps.
pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
switch PRESET
    case 'ACRE'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.467341, -87.015762; ...
            40.501484, -86.979905]);
    case 'ACRE_EXACT'
        % Read in the ACRE boundary.
        [simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ~] ...
            = extractBoundaryFromKmzFile( ...
            fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar', 'ACRE', 'AcreExactBoundaryRaw', 'ACRE.kmz'));
    case 'ACRE_LORA_5MILE_R'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.404306, -87.086317; ...
            40.547160, -86.896734]);
    case 'ACRE_LORA_1MILE_R'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.460592, -87.010842; ...
            40.489705, -86.973423]);
    case 'ACRE_LORA_HALF_MILE_R'
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.467907, -87.001530; ...
            40.482402, -86.982589]);
    case 'ACRE_LORA_20KM_R'
        % We have calculated the area based on the gateway location.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.294772, -87.227430; ...
            40.655099, -86.755491]);
    case 'ACRE_LORAWAN_2021'
        % Based on the RX locations stored in
        % AcreRxInfo/2021/oyster_data_for_sim.csv.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.33261869, -87.2002767; ...
            40.5264963, -86.9038917]);
    case {'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'}
        % Use a square area that covers the whole WHIN region.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.46661144, -87.01630345; ...
            40.50229251, -86.97594568]);
    case 'Tipp'
        % simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
        %     = constructUtmRectanglePolyMat(...
        %      [40.216047, -87.093700; ...
        %     40.562743, -86.695913]);
        tippBoundary = load(fullfile(pathToStateInfoFolder, ...
            'Tipp', 'boundary.mat'));
        assert(strcmp(simConfigs.UTM_ZONE, ...
            tippBoundary.boundary.UTM_ZONE), 'UTM zone mismatch!');
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = tippBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    case {'ExtendedTipp', 'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'}
        % Note that for the WHIN weather station case, we will set the
        % simulation boundary, but the grid points will be set to the
        % weather station locations. Please refer to simulateCoverage.m for
        % more information.
        extTippBoundary = load(fullfile(pathToStateInfoFolder, ...
            'Tipp_Extended', 'boundary.mat'));
        assert(strcmp(simConfigs.UTM_ZONE, ...
            extTippBoundary.boundary.UTM_ZONE), 'UTM zone mismatch!');
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = extTippBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    case {'ACRE_ML'}
        % Use a square area with boundary coordinates.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.438166, -87.049971; ...
            40.480526, -86.961631]);
        % [40.4647, -87.0034; 40.4921, -86.9674]);
    case {'HappyHollowPark_ML'}
        % Use a square area with boundary coordinates.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.426118, -86.908957; ...
            40.448314, -86.859313]);
        % [40.431564, -86.904384; 40.444903, -86.896195]
    case {'LindbergVillage_ML'}
        % Use a square area with boundary coordinates.
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.438707, -86.990153;...
            40.466268, -86.951015]);
        % [40.447495, -86.974048; 40.457345, -86.960447]);
    case {'IN', 'ShrinkedIN', 'ShrinkedWHIN'}
        % For ShrinkedIN, we will start with the IN boundary and shrink it
        % after the maximum radius to consider around a cell tower is
        % specified.
        inBoundary = load(fullfile(pathToStateInfoFolder, ...
            'IN', 'boundary.mat'));
        assert(strcmp(simConfigs.UTM_ZONE, ...
            inBoundary.boundary.UTM_ZONE), 'UTM zone mismatch!');
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    case {'CustomSim'}
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = CustomSimMeta.UTM_X_Y_BOUNDARY_OF_INTEREST;
    case {'UniOfCoBoulderCampus'}
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [39.9986, -105.2801; ...
            40.0138, -105.2532], deg2utm_speZone);
    case {'UniOfCoBoulderCampus_NorlinQuad_1', ...
            'UniOfCoBoulderCampus_NorlinQuad_2', ...
            'UniOfCoBoulderCampus_NorlinQuad_3', ...
            'UniOfCoBoulderCampus_NorlinQuad_4'}
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.00828571, -105.2736924; ...
            40.00947429, -105.2712728], deg2utm_speZone);
    otherwise
        error(['Unsupported preset "', PRESET, '"!'])
end

%   - Carrier frequency and wavelength. Typical values:
%       - 915 MHz
%         For LoRaWAN.
%       - 1900 MHz
%         For cellular 4G LTE
%       - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%         For cellular 5G sub 6G
%       - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%         For broadband wireless backhaul.
%       - mmWave 28000 MHz (28 GHz)
%         For cellular 5G millimeter wave
if ~exist('CARRIER_FREQUENCY_IN_MHZ', 'var')
    if ismember(PRESET, {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
            'ACRE_LORA_HALF_MILE_R', 'ACRE_LORA_20KM_R', ...
            'ACRE_LORAWAN_2021', ...
            'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN', ...
            'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'})
        % LoRa simulations.
        CARRIER_FREQUENCY_IN_MHZ = 915;
    elseif strcmpi(PRESET, 'CustomSim')
        CARRIER_FREQUENCY_IN_MHZ = CustomSimMeta.CARRIER_FREQUENCY_IN_MHZ;
    elseif strcmpi('UniOfCoBoulderCampus_NorlinQuad_', ...
            PRESET(1:min(32,length(PRESET))))
        CARRIER_FREQUENCY_IN_MHZ = 2000;
    else
        % The cellular carrier frequency can be specified here.
        CARRIER_FREQUENCY_IN_MHZ = 1900;
    end
end
simConfigs.CARRIER_FREQUENCY_IN_MHZ = CARRIER_FREQUENCY_IN_MHZ;
% Clear the user-set parameter CARRIER_FREQUENCY_IN_MHZ to make sure it is
% properly loaded for each simulation.
clearvars CARRIER_FREQUENCY_IN_MHZ;

simConfigs.CARRIER_WAVELENGTH_IN_M ...
    = physconst('LightSpeed')/simConfigs.CARRIER_FREQUENCY_IN_MHZ/1e6;

%   - We will use this number of pixels for the longer side (width/height)
%   of the map; the number of pixels for the other side will be
%   proportional to its length. Typical values:
%       - 100
%         For normal simulations.
%       - 256
%         For high-resolution simulations.
if exist('NUM_OF_PIXELS_FOR_LONGER_SIDE', 'var')
    simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE ...
        = NUM_OF_PIXELS_FOR_LONGER_SIDE;
else
    simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE = 100;
end

%   - The guaranteed spacial resolution for terrain profiles; a larger
%   value will decrease the simulation time but small obstacles may be
%   ingored. For the USGS 1/3 arc-second terrain elevation data set, we
%   have a grid resolution of approximately 10 meters.
simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M = 1.5;
%   - Similarly, the guaranteed spacial resolution for LiDAR profiles. For
%   the 2019 DSM LiDAR data, we have a grid resolution of 5 feet (~1.52 m).
%   Note: set the two profile resolutions to be the same will help speed up
%   the simulation.
simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M = 1.5;

%   - The guaranteed minimum number of LiDAR z (or possibly elevation)
%   elements in one terrain profile; this will ensure non-empty terrain
%   profiles.
simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE = 10;

%   - Rx heights for the simulator to inspect. Note that according to FAA,
%   the maximum allowable altitude is 400 feet (~122m) above the ground.
%   Typical values:
%       - [1.5; 3; 5; 7.5; (10:10:120)'; 125] for cellular and millimeter
%         wave inspection
%        - 1.5 for 3500 MHz (3.5 GHz) areostat application
%       - 0.1 for 915 MHz areostat application
switch PRESET
    case {'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'}
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
            = [1.5; 2.5];
    case {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
            'ACRE_LORA_HALF_MILE_R', 'ACRE_LORA_20KM_R', ...
            'ACRE_LORAWAN_2021'}
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
            = 0:0.25:5;
    case {'ACRE_LORA_TRAILER', 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'}
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = 1.5;
    case {'CustomSim'}
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
            = CustomSimMeta.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;
    case {'UniOfCoBoulderCampus', ...
            'UniOfCoBoulderCampus_NorlinQuad_1', ...
            'UniOfCoBoulderCampus_NorlinQuad_2', ...
            'UniOfCoBoulderCampus_NorlinQuad_3', ...
            'UniOfCoBoulderCampus_NorlinQuad_4', ...
            'ACRE_ML', 'HappyHollowPark_ML', 'LindbergVillage_ML'}
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = 1.5;
    otherwise
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
            = [1.5; 3; 5; 7.5; (10:10:120)'; 125];
end

%   - The function to calculate, according to the TX and RX heights, the
%   maximum cellular tower coverage radius. For example, it can be set to
%   be the distance that one can see at the top of the cell tower to the RX
%   above the horizon (i.e. without blockage of the earth):
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

%   - The maximum radius that a cellular tower can cover; we will use this
%   to find the cellular towers that are effective and limit the area to
%   consider for each effective cellular tower during the simulation; by
%   default (when it is not a positive scalar) we will compute a value
%   based on the radio horizon.
%
%   Please see the field 'MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY' for more
%   information.
simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M ...
    = simConfigs.getCellCoverageRadiusInM( ...
    DEFAULT_TX_ANT_HEIGHT_IN_M, PEDESTRIAN_HEIGHT_IN_M);

% Shrink the outline for PRESET 'ShrinkedIN'.
if strcmpi(PRESET, 'ShrinkedIN') || strcmpi(PRESET, 'ShrinkedWHIN')
    utmXYPolyAreaOfInt = polybuffer( ...
        polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
        -simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = utmXYPolyAreaOfInt.Vertices;
end

% Refine the outline for PRESET 'ShrinkedWHIN'.
if strcmpi(PRESET, 'ShrinkedWHIN')
    utmXYPolyShrinkedIN = ...
        polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST);

    % We specify the UTM zone to use to make sure the WHIN boundary is
    % expressed in the same zone.
    [~, whinBoundaryXYs, ~] ...
        = loadWhinBoundary(simConfigs.UTM_ZONE);
    utmXYPolyExtendedTipp = polyshape(whinBoundaryXYs);

    utmXYPolyAreaOfInt ...
        = intersect(utmXYPolyShrinkedIN, utmXYPolyExtendedTipp);
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST = utmXYPolyAreaOfInt.Vertices;
end

%   - The clearance ratio of the first Fresnel zone for a LoS path: at
%   least this ratio of the first Fresnel zone needs to be clear for a path
%   to be considered as "line of sight" (LoS); we expect the value to be
%   larger or equal to 50%.
simConfigs.LOS_FIRST_FRES_CLEAR_RATIO = 0.6;

%   - For the NTIA eHata library: the quantile percent not exceeded of the
%   signal. Limits: 0 < reliability < 1.
simConfigs.NTIA_EHATA_RELIABILITY = 0.95; % 0.5;
%   - For the NTIA eHata library: the NLCD environment code.
simConfigs.NTIA_EHATA_ENVIRO_CODE = 82; % Cultivated Crops.

%   - Functions to convert GPS degrees (lat, lon) from/to UTM (x, y). We
%   will populate them later.
simConfigs.deg2utm_speZone = deg2utm_speZone;
simConfigs.utm2deg_speZone = utm2deg_speZone;

%   - For adjusting the feedback frequency in parfor workers.
simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT = 0.25;

%   - For plotting, we need an expected value for the maximum allowed path
%   loss. With 20 MHz bandwidth, 9 dB RX noise figure, and 100 W TX power,
%   we have a maximum path loss of ~142 dB.
simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB = [0, 150];

%   - For figures, we need to decide whether to resize them for
%   publication. If this is set to be true, figures will be generated with
%   a smaller size for papers.
simConfigs.RESIZE_FIG_FOR_PUBLICATION = false;

% Enable ITM point-to-point path loss evaluation. If the field
% itmParameters is not define, the simulator will skip the ITM computation.
%
% Default values:
%   climate                  5            [Continental Temperate]
%    N_0                     301.00       (N-Units)
%   pol                      1            [Vertical]
%    epsilon                 15           [Average Ground]
%   sigma                    0.005        [Average Ground]
%    mdvar                   2            [Mobile Mode]
%   time                     95
%    location                95
%   situation                95
% Other parameters needed for the NTIA point-to-point ITM function includ:
% h_tx__meter, h_rx__meter, f__mhz, and terrian profile pfl.
simConfigs.itmParameters = struct( ...
    'climate', 5, 'n_0', 301.00, 'pol', 1, ...
    'epsilon', 15, 'sigma', 0.005, ...
    'mdvar', 2, ... 'time', 50, 'location', 50, 'situation', 50
    'time', 95, 'location', 95, 'situation', 95);
% Set this to be 'DEM' to use bare-earth ground elevation data, or 'DSM' to
% use LiDAR z data, in the terrain profile generation for the ITM model.
if ismember(PRESET, {'UniOfCoBoulderCampus'})
    % 'DSM' for Colorado campus case to consider buildings.
    simConfigs.itmParameters.terrainProfileSource = 'DSM';
else
    % For rural cases, use 'DEM' instead because most of obstacles will be
    % vegetation.
    simConfigs.itmParameters.terrainProfileSource = 'DEM';
end

% Override simConfigs fields defined in simConfigFieldsToOverride.
if exist('simConfigFieldsToOverride', 'var')
    fieldsToOverride = fieldnames(simConfigFieldsToOverride);
    numOfFieldsToOverride = length(fieldsToOverride);
    for idxFToOverride = 1:numOfFieldsToOverride
        curFieldName = fieldsToOverride{idxFToOverride};
        simConfigs.(curFieldName) ...
            = simConfigFieldsToOverride.(curFieldName);
    end
end
% ----------------------------------
%  End of Simulation Configurations
% ----------------------------------

% The absolute path to the folder for saving the results.
if exist('NUMS_OF_PIXELS_FOR_LONGER_SIDE', 'var') ...
        && (length(NUMS_OF_PIXELS_FOR_LONGER_SIDE)>1)
    simResultsFolderName = ['Simulation_', PRESET, ...
        '_Carrier_', num2str(simConfigs.CARRIER_FREQUENCY_IN_MHZ), ...
        'MHz_NGrid_', num2str(simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE), ...
        '_LiDAR_', LIDAR_DATA_SET_TO_USE];
else
    simResultsFolderName = ['Simulation_', PRESET, ...
        '_Carrier_', num2str(simConfigs.CARRIER_FREQUENCY_IN_MHZ), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE];
end

% By default, save results to folder PostProcessingResults.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', simResultsFolderName);

if strcmpi(PRESET, 'CustomSim')
    % Overwrite the dafault pathToSaveResults if necessary.
    if isfield(CustomSimMeta, 'pathToPostProcessingResultsFolder')
        pathToSaveResults = fullfile( ...
            CustomSimMeta.pathToPostProcessingResultsFolder, ...
            simResultsFolderName);
    end

    % Extract the polygon to ingore LiDAR z values (e.g., those for power
    % lines).
    if isfield(CustomSimMeta, 'polyLatLonsToUseGroundEle')
        simConfigs.polyLatLonsToUseGroundEle = ...
            CustomSimMeta.polyLatLonsToUseGroundEle;

        % Convert (lat, lon) to UTM.
        [polyXsToUseGroundEle, polyYsToUseGroundEle] = ...
            simConfigs.deg2utm_speZone( ...
            simConfigs.polyLatLonsToUseGroundEle(:, 1), ...
            simConfigs.polyLatLonsToUseGroundEle(:, 2));
        simConfigs.polyXYsToUseGroundEle ...
            = [polyXsToUseGroundEle, polyYsToUseGroundEle];
    end

    % Append simulation index if it shows up.
    if exist('idxSim', 'var')
        pathToSaveResults = [pathToSaveResults, ...
            '_idxSim_', num2str(idxSim)];
    end
elseif strcmpi(PRESET, 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY') ...
        && exist('idxT', 'var')
    pathToSaveResults = [pathToSaveResults, '_idxT_', num2str(idxT)];
end

% Turn the diary logging function on.
dirToSaveDiary = fullfile(pathToSaveResults, 'diary.txt');
if ~exist(dirToSaveDiary, 'file')
    if ~exist(pathToSaveResults, 'dir')
        mkdir(pathToSaveResults)
    end
    fclose(fopen(dirToSaveDiary, 'w'));
end
diary(dirToSaveDiary);

disp(' ')
if strcmpi(PRESET, 'CustomSim')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Configuring the simulation for PRESET ', PRESET, ...
        ' (idxSim = ', num2str(idxSim), ')...'])
else
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Configuring the simulation for PRESET ', PRESET, ...
        ' (', num2str(simConfigs.CARRIER_FREQUENCY_IN_MHZ), ' MHz)...'])
end

% Save simConfigs if it is not yet done.
dirToSaveSimConfigs = fullfile(pathToSaveResults, 'simConfigs_Raw.mat');
if exist(dirToSaveSimConfigs, 'file')
    histSimConfigs = load(dirToSaveSimConfigs);
    [~, er1, er2] = comp_struct( ...
        histSimConfigs.simConfigs, simConfigs, 0);
    % There are three function handles that are expected to be unequal.
    numOfDiffFieldsExpected = 3;
    if length(fieldnames(er1)) ~= numOfDiffFieldsExpected ...
            || length(fieldnames(er2)) ~= numOfDiffFieldsExpected
        error(['        [', datestr(now, datetimeFormat), ...
            '] The settings for this PRESET have changed!']);
    end
else
    % Note that the simConfigs saved now only contains user-specified
    % parameters.
    save(dirToSaveSimConfigs, 'simConfigs', '-v7.3');
end

disp(['    [', datestr(now, datetimeFormat), '] Done!'])

%% Preprocessing LiDAR Data
% Note: the first time of this may take a long time, depending on the size
% of the LiDAR data set, but (1) it supports recovery from interruptions,
% and (2) once we have gone through all the data once, loading the
% information would be very fast.

simConfigs.LIDAR_DATA_SET_TO_USE = LIDAR_DATA_SET_TO_USE;

% Set the dir to find the LiDAR data set.
switch simConfigs.LIDAR_DATA_SET_TO_USE
    case 'IN_DSM_2019_DEMO'
        dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar_2019', 'IN', 'DSM_Demo');
    case 'IN_DSM_2019_DEMO_ACRE'
        dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar_2019', 'IN', 'DSM_Demo_ACRE');
    case 'IN_DSM_2019'
        dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar_2019', 'IN', 'DSM');
    case 'CO_DRCOG_2020'
        dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, 'Lidar', ...
            'CO', 'USGS_CO DRCOG 1 2020', 'TilesTif');
    otherwise
        error(['Unkown LiDAR data set ', ...
            simConfigs.LIDAR_DATA_SET_TO_USE, '!'])
end

% Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need to
% remove preprocessIndianaLidarDataSet from path after things are done.
addpath(fullfile(pwd, 'lib', 'lidar'));
[lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone);
rmpath(fullfile(pwd, 'lib', 'lidar'));
lidarFileAbsDirs = cellfun(@(d) ...
    [dirToLidarFiles, strrep(d, '\', filesep)], ...
    lidarFileRelDirs, 'UniformOutput', false);

% Save simConfigs again if it is not yet done.
dirToSaveSimConfigs = fullfile(pathToSaveResults, 'simConfigs.mat');
if ~exist(dirToSaveSimConfigs, 'file')
    save(dirToSaveSimConfigs, 'simConfigs', '-v7.3');
end

%% Load Cellular Tower Information

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading cellular antenna information ...'])

% Note: we use "height" to indicate the vertical distance from the ground
% to the antenna; "elevation" to indicate the ground elevation; and
% "altitude" to indicate elevation+height.
cellAntsLatLonH = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1); %#ok<CSVRD>
if strcmpi(PRESET, 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY')
    cellAntsLatLonH = cellAntsLatLonH(idxT, :);
end
numAnts = size(cellAntsLatLonH, 1);
cellAntsXYH = nan(numAnts, 3);
[cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
    = deg2utm_speZone(cellAntsLatLonH(:,1), cellAntsLatLonH(:,2));
if size(cellAntsLatLonH, 2)>2
    cellAntsXYH(:,3) = cellAntsLatLonH(:,3);
end
boolsCellAntAltUnknown = isnan(cellAntsXYH(:,3));
if any(boolsCellAntAltUnknown)
    cellAntsXYH(boolsCellAntAltUnknown,3) = DEFAULT_TX_ANT_HEIGHT_IN_M;
end

disp(['    [', datestr(now, datetimeFormat), '] Done!'])

%% Simulation

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Conducting simulation ...'])

pathToSaveSimResults = fullfile(pathToSaveResults, 'simState.mat');
if exist(pathToSaveSimResults, 'file')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Loading history results ...'])
    load(pathToSaveSimResults)
else
    % Computing coverage and blockage maps.
    simState = simulateCoverage(pathToSaveResults, ...
        dirToLidarFiles, lidarFileAbsDirs, ...
        lidarFileXYCoveragePolyshapes, cellAntsXYH, simConfigs);
    save(pathToSaveSimResults, 'simState', '-v7.3');
end
disp(['    [', datestr(now, datetimeFormat), '] Done!'])

% Plot results.
if ~isunix % Adoid generating plots on the headless Linux clusters.
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Generating plots ...'])

    plotSimulationResults(pathToSaveResults, simState, simConfigs);

    if strcmpi(PRESET, 'WHIN_WEATHER_STATIONS')
        addpath('10_NetworkPlanner');
        genExtraPlotsForWhinWeatherStations;
    end

    % Update the cached simState.
    save(pathToSaveSimResults, 'simState', '-v7.3');

    disp(['    [', datestr(now, datetimeFormat), '] Done!'])
end

%% Export Results to .csv Files
% For easier cooperation,  export into .csv files the raw data of path loss
% maps for selected scenarios, e.g., ACRE LoRaWAN.
if startsWith(PRESET, {'ACRE_LORA_', 'UniOfCoBoulderCampus'})
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Exporting results to .csv files ...'])

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
        [~, lidarMatFName, ~] ...
            = fileparts(lidarMatFileAbsDirs{idxMatF});
        lidarMatFileAbsDirs{idxMatF} = fullfile(dirToLidarFiles, ...
            'MatlabCache', [lidarMatFName, '.mat']);

        assert(exist(lidarMatFileAbsDirs{idxMatF}, 'file'));
    end

    % Also, append the terrain elevation and LiDAR z values for the grid.
    [ele_usgs_m, lidar_z_m] = generateProfileSamps( ...
        simState.mapGridXYPts, utm2deg_speZone, ...
        lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
        lidarMatFileAbsDirs, 'both');
    simState.mapGridGroundEleInM = ele_usgs_m;
    simState.mapGridLiarZInM = lidar_z_m;

    % Update the cached simState.
    save(pathToSaveSimResults, 'simState', '-v7.3');

    exportMaps(pathToSaveResults, simConfigs, simState);

    disp(['    [', datestr(now, datetimeFormat), '] Done!'])
end

% Close the parallel pool if necessary.
delete(gcp('nocreate'));

diary off;

% EOF