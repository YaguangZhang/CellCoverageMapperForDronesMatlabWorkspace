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

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Initializing ...'])

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

pathToSaveCustomSimMetas = fullfile( ...
    pathToPostProcessingResultsFolder, 'CustomSimMetas.mat');

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

% For feedback.
CustomSimMetaDefault.SIM_GROUP_PRESET = SIM_GROUP_PRESET;

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Create CustomSimMetas

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Setting up simulations ...'])

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

save(pathToSaveCustomSimMetas, ...
    'SIM_GROUP_PRESET', 'PRESET', 'ABS_PATH_TO_RX_LOCS', ...
    'CustomSimMetas', 'rxGidLonLatHInMs', 'rxHInMsUnique', 'numOfSims');

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Run Sims

for idxSim = 1:length(CustomSimMetas)
    CustomSimMeta = CustomSimMetas{idxSim};

    disp(' ')
    disp(['[', datestr(now, datetimeFormat), ...
        '] Running sim for ', CustomSimMeta.SIM_GROUP_PRESET, ...
        ' with PRESET ', PRESET, ' (', ...
        num2str(idxSim), ') ...'])

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

%% Aggregate Results

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Aggregating simulation results ...'])

% The simulation would clear some of the parameters we need.
load(pathToSaveCustomSimMetas);

% Retrieve simulation results.
[rxEHataPLsInDb, rxAccBlkDistsInM, ...
    rxAccBlkDistsInM_Terrain, rxAccBlkDistsInM_Veg] ...
    = deal(cell(numOfSims, 1));
for idxSim = 1:length(CustomSimMetas)
    CustomSimMeta = CustomSimMetas{idxSim};
    curRxHInM = rxHInMsUnique(idxSim);

    % Locate the simulation result folder based on the name pattern.
    curSimResultsFolderName = ['Simulation_', PRESET, ...
        '_Carrier_', num2str(CustomSimMeta.CARRIER_FREQUENCY_IN_MHZ), ...
        'MHz_LiDAR_*', '_idxSim_', num2str(idxSim)];

    if ~exist('pathToPostProcessingResultsFolder', 'var')
        pathToPostProcessingResultsFolder ...
            = CustomSimMeta.pathToPostProcessingResultsFolder;
    else
        assert(strcmp(pathToPostProcessingResultsFolder, ...
            CustomSimMeta.pathToPostProcessingResultsFolder), ...
            'Results are expected to be saved in the same folder!');
    end

    curSimResultFolder = rdir( ...
        fullfile(pathToPostProcessingResultsFolder, ...
        curSimResultsFolderName), 'isdir');
    assert(length(curSimResultFolder)==1, ...
        'Exactly one sim result folder is expected!');

    dirToSimState = fullfile( ...
        curSimResultFolder.name, 'simState.mat');
    curSimState = load(dirToSimState, 'simState');

    assert(length(curSimState.simState.coverageMaps) == 1, ...
        'Expecting one and only one map!');

    rxEHataPLsInDb{idxSim} = curSimState.simState.coverageMaps{1};
    rxAccBlkDistsInM{idxSim} = curSimState.simState.blockageDistMaps{1};
    rxAccBlkDistsInM_Terrain{idxSim} ...
        = curSimState.simState.blockageByTerrainDistMaps{1};
    rxAccBlkDistsInM_Veg{idxSim} ...
        = curSimState.simState.blockageByVegDistMaps{1};

    if length(rxEHataPLsInDb{idxSim}) ...
            ~= sum(rxGidLonLatHInMs(:,4)==curRxHInM)
        disp(['Num of pts in current sim results: ', ...
            num2str(length(rxEHataPLsInDb{idxSim}))]);
        disp(['Num of pts expected based on RX height:', ...
            num2str(sum(rxGidLonLatHInMs(:,4)==curRxHInM))])
        error('Sim results have unexpected number of points!');
    end
end

% Arrange results into a table for easy exportation.
gid = rxGidLonLatHInMs(:,1);
lon = rxGidLonLatHInMs(:,2);
lat = rxGidLonLatHInMs(:,3);
install_height_m = rxGidLonLatHInMs(:,4);
ehata_path_loss_db = vertcat(rxEHataPLsInDb{:});
accumulate_blockage_dist_m = vertcat(rxAccBlkDistsInM{:});
accumulate_blockage_dist_m_terrain = vertcat(rxAccBlkDistsInM_Terrain{:});
accumulate_blockage_dist_m_vegetation = vertcat(rxAccBlkDistsInM_Veg{:});

simResultsTable = table(gid, lon, lat, install_height_m, ...
    ehata_path_loss_db, accumulate_blockage_dist_m, ...
    accumulate_blockage_dist_m_terrain, ...
    accumulate_blockage_dist_m_vegetation);

writetable(simResultsTable, ...
    fullfile(pathToPostProcessingResultsFolder, ...
    ['simResults_', SIM_GROUP_PRESET, '.csv']));

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Debugging Plots

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Aggregating simulation results ...'])

% Path loss.
curData = ehata_path_loss_db;
colorRangeToSet = [60, 150];

curDataRange = [min(curData), max(curData)];
disp(['    Current data value range: ', ...
    num2str(curDataRange(1)), ', ', num2str(curDataRange(2))]);

assert(colorRangeToSet(1)<=curDataRange(1), ...
    colorRangeToSet(2)>=curDataRange(2), ...
    'colorRangeToSet does not cover curDataRange!');

figure;
plot3k([lon, lat, curData], 'ColorRange', colorRangeToSet, ...
    'Labels', {'eHata Path Loss Overview', ...
    'Longitude (degree)', 'Latitude (degree)', '', 'Path Loss (dB)'});
plot_google_map('MapType', 'hybrid');
zlim([0, colorRangeToSet(2)]); view(2);

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_eHata');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);

% LoS blockage.
curData = accumulate_blockage_dist_m;
colorRangeToSet = [0, 8000];

curDataRange = [min(curData), max(curData)];
disp(['    Current data value range: ', ...
    num2str(curDataRange(1)), ', ', num2str(curDataRange(2))]);

assert(colorRangeToSet(1)<=curDataRange(1), ...
    colorRangeToSet(2)>=curDataRange(2), ...
    'colorRangeToSet does not cover curDataRange!');

figure;
plot3k([lon, lat, curData], 'ColorRange', colorRangeToSet, ...
    'Labels', {'Accumulate LoS Blockage Distance Overview', ...
    'Longitude (degree)', 'Latitude (degree)', '', ...
    'Blockage Distance (m)'}); view(2);
plot_google_map('MapType', 'hybrid');
zlim([0, colorRangeToSet(2)]);

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_LoSBlockageDist');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

diary off;

% EOF