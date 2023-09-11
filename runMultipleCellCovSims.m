%RUNMULTIPLECELLCOVSIMS Run multiple analyzeCellularCoverage with different
%areas/carrier frequencies of interest.
%
% This script can also be used to update the figures for completed
% simulations.
%
% Yaguang Zhang, Purdue, 05/13/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

pathToPostProcessingResultsFolder ...
    = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults');
if ~exist(pathToPostProcessingResultsFolder, 'dir')
    mkdir(pathToPostProcessingResultsFolder);
end
pathToSaveSimManDiary = fullfile( ...
    pathToPostProcessingResultsFolder, 'simManDiary.txt');
diary(pathToSaveSimManDiary);

% Shortcuts for predefined simulation groups.
%   - 'cellularCov'
%     Large-scale cellular coverage simulations.
%   - 'cellularCovExt'
%     Extensive large-scale cellular coverage simulations covering
%     different carrier frequencies with different channel models.
%   - 'whinLoRaWan'
%     The WHIN LoRaWAN + WHIN weather stations.
%   - 'acreLoRaWan'
%     The ACRE LoRaWAN.
%   - 'acreLoRaTrailer'
%     The ACRE LoRa gateways installed on the mobile trailer, simulated (1)
%     individually and (2) together.
%   - 'uniOfCoBoulder' and 'uniOfCoBoulderMedianPL'
%     The scenario with one cell tower on University of Colorado, Boulder
%     campus, with reliabilty being 0.95 and 0.5, respectively.
%   - 'diffResTipp'
%     Simulate the Tipp scenario with different resolutions of the user
%     location grid.
%   - 'diffFreqAcreML', 'diffFreqHappyHollowParkML',
%      'diffFreqLindbergVillageML'
%     Simulate the ACRE/Happy Hollow Park/Lindberg Village scenarios
%     specified by four/two/four tiles from Byunghyun with high user
%     location resolution and different carrier frequencies.
SIM_GROUP_PRESET = 'diffFreqHappyHollowParkML';

switch SIM_GROUP_PRESET
    case 'cellularCov'
        % Presets of interest.
        PRESETS = {'Tipp', 'ShrinkedWHIN', 'ShrinkedIN'};

        % Carrier frequency:
        %	- 1900 MHz
        %     For cellular 4G LTE
        CARRIER_FREQUENCIES_IN_MHZ = {1900};
    case 'cellularCovExt'
        % Presets of interest.
        PRESETS = {'Tipp', 'ShrinkedWHIN', 'ShrinkedIN'};

        % Carrier frequencies.
        %	- 1900 MHz
        %     For cellular 4G LTE.
        %   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
        %     For cellular 5G sub 6G.
        %   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
        %     For broadband wireless backhaul.
        %   - 20000 MHz (20 GHz)
        %     Max frequency supported by ITM.
        %   - mmWave 28000 MHz (28 GHz)
        %     For cellular 5G millimeter wave.
        CARRIER_FREQUENCIES_IN_MHZ ...
            = {1900, 3700, 4700, 7000, 13000, 20000, 28000};

        % Enable vegetation blockage evaluation.
        simConfigFieldsToOverride.FLAG_EVAL_LOSS_THROUGH_VEG = true;
        % % Simulate for the median path loss.
        % simConfigs.NTIA_EHATA_RELIABILITY = 50;
    case 'whinLoRaWan'
        PRESETS = {'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'};
        CARRIER_FREQUENCIES_IN_MHZ = {915};
    case 'acreLoRaWan'
        PRESETS = {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
            'ACRE_LORA_HALF_MILE_R'};
        CARRIER_FREQUENCIES_IN_MHZ = {915};
    case 'acreLoRaTrailer'
        PRESETS = {'ACRE_LORA_TRAILER', ...
            'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY'};
        CARRIER_FREQUENCIES_IN_MHZ = {915};
    case 'uniOfCoBoulder'
        PRESETS = {'UniOfCoBoulderCampus'};
        CARRIER_FREQUENCIES_IN_MHZ = {700, 1500, 3500};

        % Override the mininum distance to consider weighing FSPL.
        simConfigFieldsToOverride ...
            .MIN_TX_TO_RX_DIST_FOR_EHATA_AND_ITM_IN_M = 10;
    case 'uniOfCoBoulderMedianPL'
        PRESETS = {'UniOfCoBoulderCampus'};
        CARRIER_FREQUENCIES_IN_MHZ = {700, 1500, 3500};

        % Override the mininum distance to consider weighing FSPL.
        simConfigFieldsToOverride ...
            .MIN_TX_TO_RX_DIST_FOR_EHATA_AND_ITM_IN_M = 10;
        % Simulate for the median path loss.
        simConfigs.NTIA_EHATA_RELIABILITY = 50;
        simConfigs.itmParameters = struct( ...
            'climate', 5, 'n_0', 301.00, 'pol', 1, ...
            'epsilon', 15, 'sigma', 0.005, ...
            'mdvar', 2, ... 'time', 95, 'location', 95, 'situation', 95
            'time', 50, 'location', 50, 'situation', 50, ...
            'terrainProfileSource', 'DSM');
    case 'diffResTipp'
        % Presets of interest.
        PRESETS = {'Tipp'};

        % Carrier frequency:
        %	- 1900 MHz
        %     For cellular 4G LTE.
        CARRIER_FREQUENCIES_IN_MHZ = {1900};

        % Different numbers of pixels for the longer side of the area of
        % interest to simulate.
        NUMS_OF_PIXELS_FOR_LONGER_SIDE = 2.^(1:8);
    case 'diffFreqAcreML'
        % Presets of interest.
        PRESETS = {'ACRE_ML'};
        % Carrier frequencies from Ahmed.
        CARRIER_FREQUENCIES_IN_MHZ = ...
            {1900, ...
            624.5, 632.45, 731.5, 739, 742.5, 751, 763, 885, ...
            1940, 1957.5, 1970, 1977.5, 1982.6, 1990, ...
            2115.1, 2125, 2137.5, 2147.5, 2150, 2152.5, 2155, 2175, ...
            2506, 2525.8, 3605.3, 3620.2, 3640};

        % Number of pixels for the longer side of the area of interest to
        % simulate.
        NUMS_OF_PIXELS_FOR_LONGER_SIDE = 256;
    case 'diffFreqHappyHollowParkML'
        % Presets of interest.
        PRESETS = {'HappyHollowPark_ML'};
        % Carrier frequencies from Ahmed.
        CARRIER_FREQUENCIES_IN_MHZ = ...
            {1900, ...
            731.5, 742.5, 763, ...
            1955, 1970.1, 1982.5, 1985, ...
            2115, 2125, 2132.5, 2140, 2157.5, 2160, 2167.5, ...
            2355, 2506, 2525.8, 2538.2};

        % Number of pixels for the longer side of the area of interest to
        % simulate.
        NUMS_OF_PIXELS_FOR_LONGER_SIDE = 256;
    case 'diffFreqLindbergVillageML'
        % Presets of interest.
        PRESETS = {'LindbergVillage_ML'};
        % Carrier frequencies from Ahmed.
        CARRIER_FREQUENCIES_IN_MHZ = ...
            {1900, ...
            731.5, 742.5, 763, ...
            1955, 1970.1, 1982.5, 1985, ...
            2115, 2125, 2132.5, 2140, 2157.5, 2160, 2167.5, ...
            2355, 2506, 2525.8, 2538.2};

        % Number of pixels for the longer side of the area of interest to
        % simulate.
        NUMS_OF_PIXELS_FOR_LONGER_SIDE = 256;
    otherwise
        error(['Unknown simulation group: ', SIM_GROUP_PRESET, '!']);
end

% Configuration variable MIN_TX_TO_RX_DIST_FOR_EHATA_AND_ITM_IN_M is
% optional.
if exist('NUMS_OF_PIXELS_FOR_LONGER_SIDE', 'var')
    numOfNGrid = length(NUMS_OF_PIXELS_FOR_LONGER_SIDE);
else
    numOfNGrid = 1;
    clear NUM_OF_PIXELS_FOR_LONGER_SIDE;
end

for idxN = 1:numOfNGrid
    if exist('NUMS_OF_PIXELS_FOR_LONGER_SIDE', 'var')
        NUM_OF_PIXELS_FOR_LONGER_SIDE ...
            = NUMS_OF_PIXELS_FOR_LONGER_SIDE(idxN);
    end

    for idxFre = 1:length(CARRIER_FREQUENCIES_IN_MHZ)
        for idxPreset = 1:length(PRESETS)
            PRESET = PRESETS{idxPreset};
            CARRIER_FREQUENCY_IN_MHZ = CARRIER_FREQUENCIES_IN_MHZ{idxFre};
            disp(' ')
            disp(['[', datestr(now, datetimeFormat), ...
                '] Running sim for ', PRESET, ' (', ...
                num2str(CARRIER_FREQUENCY_IN_MHZ), ' MHz) ...'])

            try
                diary off;
                if strcmpi(PRESET, 'ACRE_LORA_TRAILER_INDIVIDUAL_GATEWAY')
                    for idxT = 1:3
                        analyzeCellularCoverage;
                    end
                else
                    analyzeCellularCoverage;
                end
                diary(pathToSaveSimManDiary);
            catch err
                diary(pathToSaveSimManDiary);
                disp(getReport(err))
                rethrow(err);
            end
            disp(['[', datestr(now, datetimeFormat), ...
                '] Done!'])
        end
    end
end

diary off;

% EOF