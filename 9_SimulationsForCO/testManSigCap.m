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
    'SigCap_LongmontCampaign_AdditionalData_20220708', ...
    'SigCap_LongmontCampaign_AdditionalData_20220812'};

for presetCell = testManPresets
    PRESET = presetCell{1};
    testSigCapAppResults;
end

%% Generate Overview .kml File

if flagGenKmlOverview
    pathToSaveOverviewKml = fileparts(pathToSaveResults);

    % For Google Earth markers.
    %   - 'CostumIcon'
    %     Icon URLs can be used to set the icon to use. Note: this will
    %     generate big .kml files which Google Earth may not well support.
    %   - 'DefaultIcon'
    %     Instead of icon URL, we will set icon color. TODO: ge_plot_new
    %     (and the function it uses, ge_plot) does not support iconColor
    %     without iconUrl.
    %   - 'Line'
    %     Instead of an icon for each point, a line will be drawn for a
    %     group of points.
    MARKER_STYLE = 'line';
    % Custom online icons can be set here.
    iconStrBase = 'http://maps.google.com/mapfiles/kml/pushpin/';
    iconUrlAtt = [iconStrBase, 'red-pushpin.png'];
    iconUrlVer = [iconStrBase, 'blue-pushpin.png'];
    iconUrlTM = [iconStrBase, 'grn-pushpin.png'];
    % Or we will use the default icon with specified colors.
    iconColorAtt = ge_color('r', 1); % Red 'FF0000FF'.
    iconColorVer = ge_color('b', 1); % Blue 'FFFF0000'.
    iconColorTM  = ge_color('g', 1); % Green 'FF00FF00'.
    % Line width. We will adjust line width and extrude amount to make sure
    % lines are not completely blocking each other visually.
    %   - TODO: somehow the colors follow the traditional AARRGGBB format
    %   instead of the KML AABBGGRR format.
    lineWidth = 6;
    lineWidthMinus = lineWidth*0.5;
    lineWidthPlus = lineWidth*1.5;
    lineHInM = 2;
    lineHInMMinus = lineHInM-1;
    lineHInMPlus = lineHInM+1;
    lineExtrude = 0;
    lineColorAtt = '80FF0000';
    lineColorVer = '800000FF';
    lineColorTM = '8000FF00';

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

        % allTimestampsDatetimeUtc ...
        %   = cachedResultsCell{idxPreset}.allTimestampsDatetimeUtc;
        numbersOfPtsForAllTracks ...
            = cachedResultsCell{idxPreset}.numbersOfPtsForAllTracks;

        switch lower(MARKER_STYLE)
            case lower('CostumIcon')
                if any(boolsIsAtt)
                    kmlPtGroupsCell{end+1} = ge_folder('ATT', ...
                        ge_point_new( ...
                        allLons(boolsIsAtt), ...
                        allLats(boolsIsAtt), 0, ...
                        'iconURL', iconUrlAtt)); %#ok<SAGROW>
                end
                if any(boolsIsVer)
                    kmlPtGroupsCell{end+1} = ge_folder('Verizon', ...
                        ge_point_new( ...
                        allLons(boolsIsVer), ...
                        allLats(boolsIsVer), 0, ...
                        'iconURL', iconUrlVer)); %#ok<SAGROW>
                end
                if any(boolsIsTM)
                    kmlPtGroupsCell{end+1} = ge_folder('T-Mobile', ...
                        ge_point_new( ...
                        allLons(boolsIsTM), ...
                        allLats(boolsIsTM), 0, ...
                        'iconURL', iconUrlTM)); %#ok<SAGROW>
                end
            case lower('DefaultIcon')
                if any(boolsIsAtt)
                    kmlPtGroupsCell{end+1} = ge_folder('ATT', ...
                        ge_point_new( ...
                        allLons(boolsIsAtt), ...
                        allLats(boolsIsAtt), 0, ...
                        'iconColor', iconColorAtt)); %#ok<SAGROW>
                end
                if any(boolsIsVer)
                    kmlPtGroupsCell{end+1} = ge_folder('Verizon', ...
                        ge_point_new( ...
                        allLons(boolsIsVer), ...
                        allLats(boolsIsVer), 0, ...
                        'iconColor', iconColorVer)); %#ok<SAGROW>
                end
                if any(boolsIsTM)
                    kmlPtGroupsCell{end+1} = ge_folder('T-Mobile', ...
                        ge_point_new( ...
                        allLons(boolsIsTM), ...
                        allLats(boolsIsTM), 0, ...
                        'iconColor', iconColorTM)); %#ok<SAGROW>
                end
            case lower('Line')
                % Special attention is needed to avoid connecting different
                % GPS tracks.
                maxTimeGapAllowedInOneSegInS = 60;

                curNumOfTracks = length(numbersOfPtsForAllTracks);
                [gePlotsAtt, gePlotsVer, gePlotsTM] = deal([]);
                ptCnt = 0;

                for idxT = 1:curNumOfTracks
                    curNumOfPts = numbersOfPtsForAllTracks(idxT);

                    curPtIndexRange = ptCnt+(1:curNumOfPts);
                    ptCnt = ptCnt + curNumOfPts;

                    curBoolsIsAtt = boolsIsAtt(curPtIndexRange);
                    curBoolsIsVer = boolsIsVer(curPtIndexRange);
                    curBoolsIsTM = boolsIsTM(curPtIndexRange);

                    if any(curBoolsIsAtt)
                        gePlotsAtt = strcat(gePlotsAtt, ge_plot3( ...
                            allLons(curBoolsIsAtt), ...
                            allLats(curBoolsIsAtt), ...
                            ones(size(allLons( ...
                            curBoolsIsAtt))).*lineHInMMinus, ...
                            'name', ['Track ', num2str(idxT)], ...
                            'lineColor', lineColorAtt, ...
                            'lineWidth', lineWidthPlus, ...
                            'extrude', lineExtrude));
                    end
                    if any(curBoolsIsVer)
                        gePlotsVer = strcat(gePlotsVer, ge_plot3( ...
                            allLons(curBoolsIsVer), ...
                            allLats(curBoolsIsVer), ...
                            ones(size(allLons( ...
                            curBoolsIsVer))).*lineHInM, ...
                            'name', ['Track ', num2str(idxT)], ...
                            'lineColor', lineColorVer, ...
                            'lineWidth', lineWidth, ...
                            'extrude', lineExtrude));
                    end
                    if any(curBoolsIsTM)
                        gePlotsTM = strcat(gePlotsTM, ge_plot3( ...
                            allLons(curBoolsIsTM), ...
                            allLats(curBoolsIsTM), ...
                            ones(size(allLons( ...
                            curBoolsIsTM))).*lineHInMPlus, ...
                            'name', ['Track ', num2str(idxT)], ...
                            'lineColor', lineColorTM, ...
                            'lineWidth', lineWidthMinus, ...
                            'extrude', lineExtrude));
                    end
                end
                assert(ptCnt == sum(numbersOfPtsForAllTracks), ...
                    'Unexpected number of processed points!');

                if any(boolsIsAtt)
                    kmlPtGroupsCell{end+1} = ge_folder('ATT', ...
                        gePlotsAtt); %#ok<SAGROW>
                end
                if any(boolsIsVer)
                    kmlPtGroupsCell{end+1} = ge_folder('Verizon', ...
                        gePlotsVer); %#ok<SAGROW>
                end
                if any(boolsIsTM)
                    kmlPtGroupsCell{end+1} = ge_folder('T-Mobile', ...
                        gePlotsTM); %#ok<SAGROW>
                end
        end
        kmlPresetFolders{idxPreset} ...
            = ge_folder(testManPresets{idxPreset}, ...
            strcat(kmlPtGroupsCell{:}));
    end

    kmlFileName = fullfile(pathToSaveOverviewKml, ...
        ['Overview_MeasLocs_', ...
        char(datetime('today', 'Format', 'yyyyMMdd')), '.kml']);

    ge_output(kmlFileName, ge_folder('Datasets', ...
        strcat(kmlPresetFolders{:}) ...
        ));
end

% EOF