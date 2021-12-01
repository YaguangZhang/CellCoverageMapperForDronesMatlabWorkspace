%TESTMANSIGCAP A manager to run testSigCapAppResults multiple
%times for specified PRESETs.
%
% Yaguang Zhang, Purdue, 11/29/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

%% Call the Test Script

testManPresets = {'SigCap_LongmontCampaign_20211105', ...
    'SigCap_LongmontCampaign_20211105_Longmont', ...
    'SigCap_LongmontCampaign_20211105_LongmontExt'};
for presetCell = testManPresets
    PRESET = presetCell{1};
    testSigCapAppResults;
end

% EOF