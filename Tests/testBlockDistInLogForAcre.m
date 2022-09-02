% TESTBLOCKDISTINLOGFORACRE A snippet to generate some maps with blockage
% distance colored in log scale.
%
% Yaguang Zhang, Purdue, 08/30/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Load the Simulation Results

dirToSimResultFolder = fullfile('..', 'PostProcessingResults', ...
    'Simulation_ACRE_LORA_20KM_R_Carrier_915MHz_LiDAR_IN_DSM_2019');
assert(exist(dirToSimResultFolder, 'dir'), ...
    'Simulation results not found!');

load(fullfile(dirToSimResultFolder, 'simConfigs.mat'));
load(fullfile(dirToSimResultFolder, 'simState.mat'));

dirToSaveFigs = fullfile(dirToSimResultFolder, 'blkDistMapsInLogScale');
if ~exist(dirToSaveFigs, 'dir')
    mkdir(dirToSaveFigs);
end

%% Extract and Plot the Blockage Dist

numOfRxHs = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);
for idxH = 1:numOfRxHs
    curRxHInM = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
    curDstBlkInM = simState.blockageDistMaps{idxH};

    hFig = figure; colormap turbo;
    curDstBlkInM(curDstBlkInM==0) = nan;
    plot3k([simState.mapGridLatLonPts(:,2:-1:1), log10(curDstBlkInM)], ...
        'Labels', {' ', 'Longitude (degree)', 'Latitude (degree)', '', ...
        'Blockage Distance in Bel Meters                              '}, ...
        'ColorRange', [0.1, 4.4]);
    view(2);
    plot_google_map('MapType', 'hybrid');
    tightfig;

    curFigFileName = ['BlkDistMapInLogScale_RxH_', ...
        num2str(curRxHInM), '_m'];
    saveas(hFig, fullfile(dirToSaveFigs, [curFigFileName, '.fig']));
    saveas(hFig, fullfile(dirToSaveFigs, [curFigFileName, '.jpg']));
    saveas(hFig, fullfile(dirToSaveFigs, [curFigFileName, '.png']));
end

% EOF