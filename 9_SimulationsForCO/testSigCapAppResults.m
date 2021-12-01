%TESTSIGCAPAPPRESULTS A script to read in and visualize results from the
%SigCap app.
%
% Yaguang Zhang, Purdue, 11/29/2021
if exist('testManPresets', 'var')
    clearvars -except PRESET testManPresets;
else
    clear;
end
clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

%% Script Parameters

% The absolute path to the folder with the FCC speed test app results.
%   - SigCap_LongmontCampaign_20211105
%     Results from the 2021 Longmont measurement campaign.
if ~exist('PRESET', 'var')
    PRESET = 'SigCap_LongmontCampaign_20211105';
end
if ismember(PRESET, ...
        {'SigCap_LongmontCampaign_20211105_Longmont', ...
        'SigCap_LongmontCampaign_20211105_LongmontExt'})
    folderNameForTestResults = 'SigCap_LongmontCampaign_20211105';
else
    folderNameForTestResults = PRESET;
end
pathToSigCapAppResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'FCC Rural Broadband', folderNameForTestResults);

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '9_SimulationsForCO', ...
    'SigCapAppTests', PRESET);
% Create directories if necessary.
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults);
end

% When specified, "latLonsForAreaOfInterest" will be used to filter out
% tests based on their location information.
switch PRESET
    case 'SigCap_LongmontCampaign_20211105_Longmont'
        latLonsForAreaOfInterest = ptPairToRect([ ...
            40.062060015391786, -105.25808569083334, ...
            40.277115808065005, -104.94420519323936]);
    case 'SigCap_LongmontCampaign_20211105_LongmontExt'
        latLonsForAreaOfInterest = ptPairToRect([ ...
            39.80668743376393, -105.38602021079713, ...
            40.40409728561462, -104.68631223445952]);
end

%% Read All the Data

logsCsvs = rdir(fullfile(pathToSigCapAppResults, '**', '*.csv'));
assert(~isempty(logsCsvs), 'Error! No CSV log files found!')

logsStruct = struct('carrier', {}, ...
    'latitude', {}, 'longitude', {}, 'lte_primary_rsrp', {});
numOfLogs = length(logsCsvs);
for idxLog = 1:numOfLogs
    % Read logs and extract the information of interest.
    logTable = readtable(logsCsvs(idxLog).name, ...
        'VariableNamingRule', 'preserve');
    for idxRecord = 1:size(logTable,1)
        logsStruct(end+1) = struct( ...
            'carrier', logTable.carrier{idxRecord}, ...
            'latitude', logTable.latitude(idxRecord), ...
            'longitude', logTable.longitude(idxRecord), ...
            'lte_primary_rsrp', logTable.lte_primary_rsrp(idxRecord) ...
            ); %#ok<SAGROW>
    end
end

allLats = vertcat(logsStruct.latitude);
allLons = vertcat(logsStruct.longitude);
allPriRsrps = vertcat(logsStruct.lte_primary_rsrp);
allCarriers = {logsStruct.carrier}';

% If necessary, remove results out of the area of interest.
if exist('latLonsForAreaOfInterest', 'var')
    boolsToKeep = inpolygon(allLats, allLons, ...
        latLonsForAreaOfInterest(:,1), latLonsForAreaOfInterest(:,2));
    logsStruct = logsStruct(boolsToKeep);
    clearvars latLonsForAreaOfInterest;

    allLats = vertcat(logsStruct.latitude);
    allLons = vertcat(logsStruct.longitude);
    allPriRsrps = vertcat(logsStruct.lte_primary_rsrp);
    allCarriers = {logsStruct.carrier}';
end

%% Plots

% RSRP on map.
hRsrp = figure; hold on;
plot3k([allLons, allLats, allPriRsrps], ...
    'Labels', {'LTE Primary RSRP', '', '', '', ...
    'RSRP (dBm)'}, ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hRsrp, fullfile(pathToSaveResults, 'primRsrp_3D.jpg'));
view(2);
saveas(hRsrp, fullfile(pathToSaveResults, 'primRsrp_2D.jpg'));
saveas(hRsrp, fullfile(pathToSaveResults, 'primRsrp_2D.fig'));

% RSRP by carrier.
boolsIsAtt = cellfun(@(v) strcmp(v(1:3), 'AT&'), allCarriers);
boolsIsVer = cellfun(@(v) strcmp(v(1:3), 'Ver'), allCarriers);
boolsIsTM = cellfun(@(v) strcmp(v(1:3), 'T-M'), allCarriers);
assert(all(boolsIsAtt|boolsIsVer|boolsIsTM), ...
    'Unknown carrier label(s)!');

markerForAtt = ':r.';
markerForVer = ':b.';
markerForTM = ':g.';

hRsrpByCar = figure; hold on;
hMAtt = stem3(allLons(boolsIsAtt), allLats(boolsIsAtt), ...
    -allPriRsrps(boolsIsAtt), markerForAtt);
hMVer = stem3(allLons(boolsIsVer), allLats(boolsIsVer), ...
    -allPriRsrps(boolsIsVer), markerForVer);
hMTM = stem3(allLons(boolsIsTM), allLats(boolsIsTM), ...
    -allPriRsrps(boolsIsTM), markerForTM);
plot_google_map('MapType', 'Hybrid');
title('Amplitude of LTE Primary RSRP by Carriers');
zlabel('-RSRP (dBm)');
xticklabels([]); yticklabels([]); view(3);
legend([hMAtt, hMVer, hMTM], 'AT&T', 'Verizon', 'T-Mobile', ...
    'Location','best');

saveas(hRsrpByCar, ...
    fullfile(pathToSaveResults, 'primRsrpByCar_3D.jpg'));
view(2);
legend([hMAtt, hMVer, hMTM], 'AT&T', 'Verizon', 'T-Mobile', ...
    'No Service', 'Location','best');
saveas(hRsrpByCar, ...
    fullfile(pathToSaveResults, 'primRsrpByCar_2D.jpg'));
saveas(hRsrpByCar, ...
    fullfile(pathToSaveResults, 'primRsrpByCar_2D.fig'));

hRsrpAtt = figure; hold on;
plot3k([allLons(boolsIsAtt), allLats(boolsIsAtt), ...
    allPriRsrps(boolsIsAtt)], ...
    'Labels', {'LTE Primary RSRP for ATT', '', '', '', ...
    'RSRP (dBm)'}, 'ColorRange', [min(allPriRsrps), max(allPriRsrps)], ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hRsrpAtt, fullfile(pathToSaveResults, 'primRsrp_att_3D.jpg'));
view(2);
saveas(hRsrpAtt, fullfile(pathToSaveResults, 'primRsrp_att_2D.jpg'));
saveas(hRsrpAtt, fullfile(pathToSaveResults, 'primRsrp_att_2D.fig'));

hRsrpVer = figure; hold on;
plot3k([allLons(boolsIsVer), allLats(boolsIsVer), ...
    allPriRsrps(boolsIsVer)], ...
    'Labels', {'LTE Primary RSRP for Verizon', '', '', '', ...
    'RSRP (dBm)'}, 'ColorRange', [min(allPriRsrps), max(allPriRsrps)], ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hRsrpVer, fullfile(pathToSaveResults, 'primRsrp_ver_3D.jpg'));
view(2);
saveas(hRsrpVer, fullfile(pathToSaveResults, 'primRsrp_ver_2D.jpg'));
saveas(hRsrpVer, fullfile(pathToSaveResults, 'primRsrp_ver_2D.fig'));

hRsrpTM = figure; hold on;
plot3k([allLons(boolsIsTM), allLats(boolsIsTM), ...
    allPriRsrps(boolsIsTM)], ...
    'Labels', {'LTE Primary RSRP for T-Mobile', '', '', '', ...
    'RSRP (dBm)'}, 'ColorRange', [min(allPriRsrps), max(allPriRsrps)], ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hRsrpTM, fullfile(pathToSaveResults, 'primRsrp_tm_3D.jpg'));
view(2);
saveas(hRsrpTM, fullfile(pathToSaveResults, 'primRsrp_tm_2D.jpg'));
saveas(hRsrpTM, fullfile(pathToSaveResults, 'primRsrp_tm_2D.fig'));
% EOF