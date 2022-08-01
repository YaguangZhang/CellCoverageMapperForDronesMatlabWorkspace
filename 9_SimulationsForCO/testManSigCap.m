%TESTMANSIGCAP A manager to run testSigCapAppResults multiple times for
%specified PRESETs.
%
% Yaguang Zhang, Purdue, 11/29/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

% Set this to be true to generate an overview .kml file for all data
% processed (set by testManPresets).
flagGenKmlOverview = true;

%% Call the Test Script

% Please refer to testSigCapAppResults.m for a list of supported data sets.
testManPresets = {'SigCapCSV20210406', 'SigCapCSV20210928_Anderson', ...
    'SigCap_LongmontCampaign_20211105_Longmont', ...
    'SigCap_LongmontCampaign_20211105_LongmontExt', ...
    'SigCap_LongmontCampaign_AdditionalData_20220708'};

for presetCell = testManPresets
    PRESET = presetCell{1};
    testSigCapAppResults;
end

%% Generate Overview .kml File

if flagGenKmlOverview
    pathToSaveOverviewKml = fileparts(pathToSaveResults);

    % For Google Earth markers.
    iconStrBase = 'http://maps.google.com/mapfiles/kml/pushpin/';
    iconUrlAtt = [iconStrBase, 'red-pushpin.png'];
    iconUrlVer = [iconStrBase, 'blue-pushpin.png'];
    iconUrlTM = [iconStrBase, 'grn-pushpin.png'];

    % Reload cache files and show results on an overview .kml map for
    % Google Earth.
    numOfPresets = length(testManPresets);
    cachedResultsCell = cell(numOfPresets, 1);

    for idxPreset = 1:numOfPresets
        curPreset = testManPresets{idxPreset};

        cachedResultsCell{idxPreset} = load(fullfile( ...
            pathToSaveOverviewKml, curPreset, 'cache.mat'));
    end

    % Generate .kml string for measurement locations with Google Earth
    % markers.
    kmlPresetFolders = cell(numOfPresets, 1);
    for idxPreset = 1:numOfPresets
        kmlPtGroupsCell = {};
        allLats = cachedResultsCell{idxPreset}.allLats;
        allLons = cachedResultsCell{idxPreset}.allLons;

        boolsIsAtt = cachedResultsCell{idxPreset}.boolsIsAtt;
        boolsIsVer = cachedResultsCell{idxPreset}.boolsIsVer;
        boolsIsTM = cachedResultsCell{idxPreset}.boolsIsTM;

        if any(boolsIsAtt)
            kmlPtGroupsCell{end+1} = ge_folder('ATT', ...
                ge_point_new( ...
                allLons(boolsIsAtt), ...
                allLats(boolsIsAtt), 0,...
                'iconURL', iconUrlAtt)); %#ok<SAGROW>
        end

        if any(boolsIsVer)
            kmlPtGroupsCell{end+1} = ge_folder('Verizon', ...
                ge_point_new( ...
                allLons(boolsIsVer), ...
                allLats(boolsIsVer), 0,...
                'iconURL', iconUrlVer)); %#ok<SAGROW>
        end

        if any(boolsIsTM)
            kmlPtGroupsCell{end+1} = ge_folder('T-Mobile', ...
                ge_point_new( ...
                allLons(boolsIsTM), ...
                allLats(boolsIsTM), 0,...
                'iconURL', iconUrlTM)); %#ok<SAGROW>
        end

        kmlPresetFolders{idxPreset} ...
            = ge_folder(testManPresets{idxPreset}, ...
            strcat(kmlPtGroupsCell{:}));
    end

    kmlFileName = fullfile(pathToSaveOverviewKml, ...
        'Overview_MeasLocs.kml');

    ge_output(kmlFileName, ge_folder('Datasets', ...
        strcat(kmlPresetFolders{:}) ...
        ));
end

% EOF