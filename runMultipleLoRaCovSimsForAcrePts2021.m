%RUNMULTIPLELORACOVSIMSFORACREPTS2021 Run multiple analyzeCellularCoverage
%for point-by-point simulation results for the 2021 ACRE LoRaWAN dataset.
%
% This script can also be used to update the figures for completed
% simulations.
%
% Yaguang Zhang, Purdue, 06/20/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Initialization

% Predefined simulation group label.
%   - 'acreLoRaWan2021'
%     The ACRE LoRaWAN results collected in 2021: one tower gateway +
%     multiple vechiles.
SIM_GROUP_PRESET = 'AcreLoRaWan2021';

pathToPostProcessingResultsFolder ...
    = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    ['PostProcessingResults_', SIM_GROUP_PRESET]);
if ~exist(pathToPostProcessingResultsFolder, 'dir')
    mkdir(pathToPostProcessingResultsFolder);
end
pathToSaveSimManDiary = fullfile( ...
    pathToPostProcessingResultsFolder, 'simManDiary.txt');
diary(pathToSaveSimManDiary);

% Preset of interest.
PRESET = 'CustomSim';

%------------------------
% Required for CustomSim
%------------------------
% Carrier frequencies.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%     For broadband wireless backhaul.
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CustomSimMetaDefault.CARRIER_FREQUENCY_IN_MHZ = 915;
% Tower location.
CustomSimMetaDefault.ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
    ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'PurdueAcreLoraWanTowers.csv');
% Area of interest is to be extended Tipp.
pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
extTippBoundary = load(fullfile(pathToStateInfoFolder, ...
    'Tipp_Extended', 'boundary.mat'));
CustomSimMetaDefault.UTM_X_Y_BOUNDARY_OF_INTEREST ...
    = extTippBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;

% Rx location information TBD.
CustomSimMetaDefault.mapGridLatLonPts = nan;
CustomSimMetaDefault.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = nan;

%------------------------
% Optional for CustomSim
%------------------------
% Enable simulations for vegetation.
CustomSimMetaDefault.FLAG_EVAL_LOSS_THROUGH_VEG = true;
% Dir to save results.
CustomSimMetaDefault.pathToPostProcessingResultsFolder ...
    = pathToPostProcessingResultsFolder;

%% Create CustomSimMetas

ABS_PATH_TO_RX_LOCS = fullfile( ...
    ABS_PATH_TO_SHARED_FOLDER, ...
    'AcreRxInfo', '2021', 'oyster_data_for_sim.csv');
rxGidLonLatHInMs = readmatrix(ABS_PATH_TO_RX_LOCS);

rxHInMsUnique = unique(rxGidLonLatHInMs(:,4));
numOfSims = length(rxHInMsUnique);
CustomSimMetas = cell(numOfSims, 1);

for idxSim = 1:length(CustomSimMetas)
    curRxHInM = rxHInMsUnique(idxSim);

    CustomSimMetas{idxSim} = CustomSimMetaDefault;
    CustomSimMetas{idxSim}.mapGridLatLonPts = rxGidLonLatHInMs( ...
        rxGidLonLatHInMs(:,4)==curRxHInM, [3,2]);
    CustomSimMetas{idxSim}.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = curRxHInM;
end

%% Run Sims

for idxSim = 1:length(CustomSimMetas)
    disp(' ')
    disp(['[', datestr(now, datetimeFormat), ...
        '] Running sim for ', SIM_GROUP_PRESET, ...
        ' with PRESET ', PRESET, ' (', ...
        num2str(idxSim), ') ...'])

    CustomSimMeta = CustomSimMetas{idxSim};
    try
        diary off;
        analyzeCellularCoverage;
        diary(pathToSaveSimManDiary);
    catch err
        diary(pathToSaveSimManDiary);
        disp(getReport(err))
        rethrow(err);
    end

    disp(['[', datestr(now, datetimeFormat), ...
        '] Done!'])
end

diary off;

% EOF