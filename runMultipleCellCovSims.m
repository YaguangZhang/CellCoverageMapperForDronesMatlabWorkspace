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
%   - 'whinLoRaWan'
%     The WHIN LoRaWAN + WHIN weather stations.
%   - 'acreLoRaWan'
%     The ACRE LoRaWAN.
SIM_GROUP_PRESET = 'cellularCov';

switch SIM_GROUP_PRESET
    case 'cellularCov'
        % Presets of interest.
        PRESETS = {'Tipp', 'ShrinkedWHIN', 'ShrinkedIN'};

        % Carrier frequencies.
        %	- 1900 MHz
        %     For cellular 4G LTE
        %   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
        %     For cellular 5G sub 6G
        %   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
        %     For broadband wireless backhaul.
        %   - mmWave 28000 MHz (28 GHz)
        %     For cellular 5G millimeter wave
        CARRIER_FREQUENCIES_IN_MHZ ...
            = {1900}; % {1900, 3700, 4700, 7000, 13000, 28000};
    case 'whinLoRaWan'
        PRESETS = {'WHIN_WEATHER_STATIONS', 'WHIN_LORAWAN'};
        CARRIER_FREQUENCIES_IN_MHZ = {915};
    case 'acreLoRaWan'
        PRESETS = {'ACRE_LORA_5MILE_R', 'ACRE_LORA_1MILE_R', ...
            'ACRE_LORA_HALF_MILE_R'};
        CARRIER_FREQUENCIES_IN_MHZ = {915};
    otherwise
        error(['Unknown simulation group: ', SIM_GROUP_PRESET, '!']);
end

for idxFre = 1:length(CARRIER_FREQUENCIES_IN_MHZ)
    for idxPreset = 1:length(PRESETS)
        PRESET = PRESETS{idxPreset};
        CARRIER_FREQUENCY_IN_MHZ = CARRIER_FREQUENCIES_IN_MHZ{idxFre};
        disp(' ')
        disp(['    [', datestr(now, datetimeFormat), ...
            '] Running sim for ', PRESET, ' (', ...
            num2str(CARRIER_FREQUENCY_IN_MHZ), ' MHz) ...'])

        try
            diary off;
            analyzeCellularCoverage;
            diary(pathToSaveSimManDiary);
        catch err
            diary(pathToSaveSimManDiary);
            disp(getReport(err))
            rethrow(err);
        end
        disp(['    [', datestr(now, datetimeFormat), ...
            '] Done!'])
    end
end

diary off;

% EOF