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

% Presets.
PRESETS = {'ACRE_EXACT', 'Tipp', 'ShrinkedIN'};
% Carrier frequencies.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700};

for idxFre = 1:length(CARRIER_FREQUENCIES_IN_MHZ)
    for idxPreset = 1:length(PRESETS)
        PRESET = PRESETS{idxPreset};
        CARRIER_FREQUENCY_IN_MHZ = CARRIER_FREQUENCIES_IN_MHZ{idxFre};
        
        analyzeCellularCoverage;
    end
end

% EOF