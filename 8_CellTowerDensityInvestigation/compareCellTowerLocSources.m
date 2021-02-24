%COMPARECELLTOWERLOCSOURCES Visualize and compare the cellular towers from
%different sources.
%
% Yaguang Zhang, Purdue, 02/19/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Parameters
% The directories to load cellular tower location information. We have:
%   - OpenCelliD
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'OpenCelliDUsa_20210219', ...
%            'Cellular_Towers_LatLonHR.csv');
%   - NTIA randomized U.S. cellular laydown
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'RandomizedCarrierSitesv2.csv')
%   - Homeland Infrastructure Foundation-Level Data (HIFLD) Cellular Towers
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'HIFLD', ...
%            'CellTowers', 'Cellular_Towers_LatLon.csv');
%   - HIFLD Land Mobile Commercial Transmission Towers
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'HIFLD', ...
%            'LandMobileCommercialTxTowers', ...
%           'Land_Mobile_Commercial_Transmission_Towers_LatLonH.csv');
ABS_PATHS_TO_CELL_ANTENNAS_CSV = { ...
    ...
    fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'OpenCelliDUsa_20210219', ...
    'Cellular_Towers_LatLonHR.csv');
    ...
    fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'RandomizedCarrierSitesv2.csv'); ...
    ...
    fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'HIFLD', ...
    'CellTowers', 'Cellular_Towers_LatLonH.csv'); ...
    ...
    fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'HIFLD', ...
    'LandMobileCommercialTxTowers', ...
    'Land_Mobile_Commercial_Transmission_Towers_LatLonH.csv')};
LABELS_SRC = {'OpenCelliD', 'NTIA', ...
    'HIFLDCellularTs', 'HIFLDLandMobCommTxTs'};

% The directory to load the simulation results for Indiana State.
pathToLoadSimResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResultsBackups', 'In_Backup_Res_50_LowestAnts', ...
    'simResults.mat');

% The directory to load the state boundary for Colorado.
pathToColoradoBound = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'Colorado_State_Boundary', ...
    'Colorado_State_Boundary.shp');

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '8_CellTowerDensityInvestigation');

%% Load Data

disp(' ')
disp('    Loading cellular tower information ...')

numOfSources = length(ABS_PATHS_TO_CELL_ANTENNAS_CSV);
assert(numOfSources ==length(LABELS_SRC), ...
    'Unexpected number of data set labels!')
% Note: we use "height" to indicate the vertical distance from the ground
% to the antenna; "elevation" to indicate the ground elevation; and
% "altitude" to indicate elevation+height.
towerLatLonHsCell = cell(numOfSources, 1);
for idxSrc = 1:numOfSources
    towerLatLonHsCell{idxSrc} = csvread( ...
        ABS_PATHS_TO_CELL_ANTENNAS_CSV{idxSrc}, 1, 1);
end

disp('    Loading the boundary for Colorado ...')
coBound = shaperead(pathToColoradoBound);

disp('    Done!')

%% Initialization

disp(' ')
disp('    Loading history simulation results ...')

load(pathToLoadSimResults);

% Functions to convert GPS degrees (lat, lon) from/to UTM (x, y).
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

% GPS boundary for the area of interest (here we have IN and CO).
labelsAreaOfInterest = {'IN', 'CO'};
numOfAreasOfInt = length(labelsAreaOfInterest);
[boundOfInterestLatsCell, boundOfInterestLonsCell] = deal(cell(2,1));
[boundOfInterestLatsCell{1}, boundOfInterestLonsCell{1}] ...
    = utm2deg_speZone(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
boundOfInterestLatsCell{2} = coBound.Y;
boundOfInterestLonsCell{2} = coBound.X;

disp('    Done!')

%% Plots

for idxA = 1:numOfAreasOfInt
    labelAreaOfInterest = labelsAreaOfInterest{idxA};
    boundOfInterestLats = boundOfInterestLatsCell{idxA};
    boundOfInterestLons = boundOfInterestLonsCell{idxA};
    
    disp(' ')
    disp(['    Generating plots for ', labelAreaOfInterest, ' ...'])
    
    for idxSrc = 1:numOfSources
        % Cellular tower on road map for US.
        hCellTowerOnRoadMap = figure('Visible', true, ...
            'Unit', 'pixel', 'Position', [0,0,1600,1200]);
        hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
        hPolyIn = plot(boundOfInterestLons, boundOfInterestLats, ...
            'r-', 'LineWidth', 7);
        newAxis = extendAxisByFactor(axis, 0.1);
        axis(newAxis);
        tightfig;
        plot_google_map('MapType', 'roadmap');
        hCellTowers = plot(towerLatLonHsCell{idxSrc}(:,2), ...
            towerLatLonHsCell{idxSrc}(:,1), ...
            'b.', 'MarkerSize', 5);
        uistack(hPolyIn,'top'); box on;
        legend([hCellTowers, hPolyIn], ...
            'Cell Tower', labelAreaOfInterest, ...
            'FontSize', 25, 'Location', 'southeast', 'LineWidth', 3);
        % Move IN state to the center.
        moveIndianaStateToCenter;
        saveas(hCellTowerOnRoadMap, ...
            fullfile(pathToSaveResults, ...
            ['cellTowerOnRoadMap_', labelAreaOfInterest, '_', ...
            LABELS_SRC{idxSrc}, '.png']));
    end
    
    % Cellular tower on road map for all sources.
    hCellTowerOnRoadMap = figure('Visible', true, ...
        'Unit', 'pixel', 'Position', [0,0,1600,1200]);
    hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
    hPolyIn = plot(boundOfInterestLons, boundOfInterestLats, ...
        'r-', 'LineWidth', 7);
    newAxis = extendAxisByFactor(axis, 0.1);
    axis(newAxis);
    tightfig;
    plot_google_map('MapType', 'roadmap');
    hCellTowersCell = cell(numOfSources, 1);
    srcOrderForPlotting = [1 2 4 3];
    mksToUse = {'b.', 'rx', 'go', 'k*'};
    for idxSrc = 1:numOfSources
        curSrcNum = srcOrderForPlotting(idxSrc);
        hCellTowersCell{idxSrc} = plot( ...
            towerLatLonHsCell{curSrcNum}(:,2), ...
            towerLatLonHsCell{curSrcNum}(:,1), ...
            mksToUse{idxSrc}, 'MarkerSize', 5, 'LineWidth', 1);
    end
    uistack(hPolyIn,'top'); box on;
    legend([hCellTowersCell{1:end}, hPolyIn], ...
        [LABELS_SRC(srcOrderForPlotting), {labelAreaOfInterest}], ...
        'FontSize', 25, 'Location', 'southeast', 'LineWidth', 3);
    % Move IN state to the center.
    moveIndianaStateToCenter;
    saveas(hCellTowerOnRoadMap, ...
        fullfile(pathToSaveResults, ...
        ['cellTowerOnRoadMap_', labelAreaOfInterest, 'All.png']));
    
    disp('    Done!')
end

% EOF