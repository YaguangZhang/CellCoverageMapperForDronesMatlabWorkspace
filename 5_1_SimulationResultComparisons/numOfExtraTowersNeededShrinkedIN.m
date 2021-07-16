% NUMOFEXTRATOWERSNEEDEDSHRINKEDIN Estimate the number of towers needed to
% cover the remaining part of the area of interest (for preset ShrinkedIN),
% with the users located at a 1.5 m height.
%
% Yaguang Zhang, Purdue, 07/15/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Parameters

pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
areaInSqKmIN = 94320;

ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs', ...
    'NtiaLayoutMergedWithHifldCellTs_Threshold_1000m_LatLonH.csv');

% Presets.
PRESET = 'ShrinkedIN';
% Carrier frequencies.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%     For broadband wireless backhaul.
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700, 7000, 13000, 28000};
% Lidar data set used in the simulations.
LIDAR_DATA_SET_USED = 'IN_DSM_2019';

% Max allowed path loss.
maxAllowedPathLossInDb = 150;

% Expected user height.
expectedRxHInM = 1.5;

% 'FSPL' or 'eHata'
cellSizeEstiMethod = 'eHata';

% The absolute paths to simulation results.
numOfCarriers = length(CARRIER_FREQUENCIES_IN_MHZ);
absPathsToSimResults = cell(numOfCarriers, 1);

for idxCarrier = 1:numOfCarriers
    absPathsToSimResults{idxCarrier} ...
        = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', PRESET, ...
        '_Carrier_', num2str(CARRIER_FREQUENCIES_IN_MHZ{idxCarrier}), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_USED]);
end

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'Simulation_Comparisons');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end
dirToSaveSummaryTable = fullfile(pathToSaveResults, ...
    ['numOfExTs_', PRESET, ...
    '_Lidar_', LIDAR_DATA_SET_USED, ...
    '_maxPL_', num2str(maxAllowedPathLossInDb), ...
    'dB_CellRangeBy_', cellSizeEstiMethod, '.csv']);

%% Get Functions for Pathloss Value Evaluation
% This is for estimating cell size.
switch cellSizeEstiMethod
    case 'FSPL'
        % Free Space Path Loss (FSPL) inverse.
        fsplInDb2distInM = @(fInMz, lossInDb) 10.^(lossInDb./20) ...
            .*(physconst('LightSpeed')./(fInMz.*(10^6))) ./ 4 ./ pi;
    case 'eHata'
        % Load the NTIA eHata library first, if necessary, to avoid the
        % "unable to find ehata" error. Note that a compiler is needed:
        %    https://www.mathworks.com/support/requirements/supported-compilers.html
        if ~libisloaded('ehata')
            loadlibrary('ehata');
        end
        
        % The distance values to inspect for creating the cell radius vs
        % path loss data set.
        distsInMToInspect = 1:20000;
        
        expectedTowerHInM = 50;
    otherwise
        error(['Unsupported cell size estimation method ', ...
            cellSizeEstiMethod, '!']);
end

%% Load Data

% GPS boundary for IN.
inBoundary = load(fullfile(pathToStateInfoFolder, ...
    'IN', 'boundary.mat'));
estAreaInSqKmIN = polyarea( ...
    inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2))./(10^6);

% Cell tower locs.
cellAntsLatLonH = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);
% Preassign memory for the UTM version of cell tower records.
numAnts = size(cellAntsLatLonH, 1);
cellAntsXYH = nan(numAnts, 3);

%% Compute Uncovered Area for Each Simulation

presets = cell(numOfCarriers, 1);
[areasInSqKm, coveredAreasInSqKm, ...
    coveredRadiiPerCellInKm, coveredAreasPerCellInSqKm, ...
    numsOfExistingTs, numsOfExtraTs] ...
    = deal(nan(numOfCarriers, 1));
maxAllowedPLsInDb = ones(numOfCarriers, 1).*maxAllowedPathLossInDb;
for idxCarrier = 1:numOfCarriers
    presets{idxCarrier} = PRESET;
    curAbsPathToSimResMat = fullfile(absPathsToSimResults{idxCarrier}, ...
        'simState.mat');
    curAbsPathToSimConfMat = fullfile(absPathsToSimResults{idxCarrier}, ...
        'simConfigs_Raw.mat');
    clearvars simState simConfig;
    load(curAbsPathToSimResMat);
    load(curAbsPathToSimConfMat);
    
    % For GPS and UTM conversions.
    [deg2utm_speZone, utm2deg_speZone] ...
        = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);
    
    % Compute the area for the region of interest considering the error
    % caused by reference system conversion.
    areasInSqKm(idxCarrier) = polyarea( ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2))./(10^6) ...
        ./estAreaInSqKmIN.*areaInSqKmIN;
    
    assert( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(1)==expectedRxHInM, ...
        ['The first RX ant height expected (', ...
        num2str(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(1)), ...
        ' m) is not ', num2str(expectedRxHInM), ' m!']);
    coveredAreasInSqKm(idxCarrier) = ...
        sum(simState.coverageMaps{1}<=maxAllowedPathLossInDb) ...
        ./length(simState.coverageMaps{1}).*areasInSqKm(idxCarrier);
    
    switch cellSizeEstiMethod
        case 'FSPL'
            coveredRadiiPerCellInKm(idxCarrier) ...
                = fsplInDb2distInM(CARRIER_FREQUENCIES_IN_MHZ{idxCarrier}, ...
                maxAllowedPathLossInDb)/(10^3);
        case 'eHata'
            numOfDistToInspect = length(distsInMToInspect);
            preCompPathLossRecords = nan(numOfDistToInspect, 1);
            for idxD = 1:numOfDistToInspect
                curDInM = distsInMToInspect(idxD);
                preCompPathLossRecords(idxD) = ...
                    computeCoveragePL([0,0,expectedTowerHInM], ...
                    [curDInM,0,expectedRxHInM], ...
                    [9, curDInM, zeros(1,10)], simConfigs, zeros(1,3));
            end
            coveredRadiiPerCellInKm(idxCarrier) = ...
                interp1(preCompPathLossRecords, ...
                distsInMToInspect, maxAllowedPathLossInDb)/1000;
    end
    
    coveredAreasPerCellInSqKm(idxCarrier) ...
        = pi*((coveredRadiiPerCellInKm(idxCarrier))^2);
    
    [cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
        = deg2utm_speZone(cellAntsLatLonH(:,1), cellAntsLatLonH(:,2));
    
    numsOfExistingTs(idxCarrier) = sum( inpolygon( ...
        cellAntsXYH(:,1), cellAntsXYH(:,2), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2)) );
    
    numsOfExtraTs(idxCarrier) = ...
        (areasInSqKm(idxCarrier) ...
        - coveredAreasInSqKm(idxCarrier)) ...
        /coveredAreasPerCellInSqKm(idxCarrier);
end

%% Put Results in a Table

carrierFreqInMHz = CARRIER_FREQUENCIES_IN_MHZ';
summaryTable = table(presets, areasInSqKm, carrierFreqInMHz, ...
    maxAllowedPLsInDb, coveredAreasInSqKm, ...
    coveredRadiiPerCellInKm, coveredAreasPerCellInSqKm, ...
    numsOfExtraTs, numsOfExistingTs);
disp(summaryTable);

% Export the results to a file.
writetable(summaryTable, dirToSaveSummaryTable);

% EOF