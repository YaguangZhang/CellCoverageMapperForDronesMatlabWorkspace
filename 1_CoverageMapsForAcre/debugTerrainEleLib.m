% DEBUGTERRAINELELIB For fixing the bug (shifted tiles) in Terrain
% Elevation Library under Linux.
%
% Yaguang Zhang, Purdue, 03/25/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', 'AcreGeoFeatureInvestigation');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

%% Test Case

lowResGridLatRange = [40.467071825508839; 40.501647105371035];
lowResGridLonRange = [-87.015786097515942; -86.976458345550711];
regionRef = fetchregion(lowResGridLatRange, lowResGridLonRange);
rawRefElevData = regionRef.readelevation( ...
    lowResGridLatRange, lowResGridLonRange, 'sampleFactor', 1);

dispelev(rawRefElevData, 'mode', 'latlong'); plot_google_map;

% EOF