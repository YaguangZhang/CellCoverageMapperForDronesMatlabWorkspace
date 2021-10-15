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

pathToSaveSimManDiary = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'simManDiary.txt');
diary(pathToSaveSimManDiary);

% Presets.
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
CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700, 7000, 13000, 28000};

for idxFre = 1:length(CARRIER_FREQUENCIES_IN_MHZ)
    for idxPreset = 1:length(PRESETS)
        PRESET = PRESETS{idxPreset};
        CARRIER_FREQUENCY_IN_MHZ = CARRIER_FREQUENCIES_IN_MHZ{idxFre};

        try
            diary off;
            analyzeCellularCoverage;
            diary(pathToSaveSimManDiary);
            delete(gcp('nocreate'));
        catch err
            diary(pathToSaveSimManDiary);
            disp(getReport(err))
        end
    end
end

diary off;

% EOF