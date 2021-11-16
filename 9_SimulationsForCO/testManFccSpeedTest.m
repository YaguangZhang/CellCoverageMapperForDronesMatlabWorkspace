%TESTMANFCCSPEEDTEST A manager to run testFccSpeedTestAppResults multiple
%times for specified PRESETs.
%
% Yaguang Zhang, Purdue, 11/16/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

%% Call the Test Script

testManPresets = {'FCC_SpeedTest_LongmontCampaign_20211105', ...
    'FCC_SpeedTest_LongmontCampaign_20211105_USNA', ...
    'FCC_SpeedTest_LongmontCampaign_20211105_Longmont', ...
    'FCC_SpeedTest_LongmontCampaign_20211105_LongmontExt'};
for presetCell = testManPresets
    PRESET = presetCell{1};
    testFccSpeedTestAppResults;
end

% EOF