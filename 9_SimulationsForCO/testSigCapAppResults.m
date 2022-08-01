%TESTSIGCAPAPPRESULTS A script to read in and visualize results from the
%SigCap app.
%
% Yaguang Zhang, Purdue, 11/29/2021
if exist('testManPresets', 'var')
    clearvars -except PRESET testManPresets flagGenKmlOverview;
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

% Supported PRESET values:
%   - SigCapCSV20210406
%     FCC preliminary measurements in Longmont.
%   - SigCapCSV20210928_Anderson
%     Test measurements from Professor Anderson, conducted in USNA.
%   - SigCap_LongmontCampaign_20211105
%     Results from the 2021 Longmont measurement campaign by Professor
%     Anderson and his students. This dataset is a copy of the SigCap
%     output files under folder "Longmont Lab_Raw".
%       - SigCap_LongmontCampaign_20211105_Longmont and
%       SigCap_LongmontCampaign_20211105_LongmontExt
%         Filtered version of SigCap_LongmontCampaign_20211105 to only show
%         results in/out of Longmont.
%   - SigCap_LongmontCampaign_AdditionalData_20220708
%     Extra results from the 2022 Longmont measurement campaign by
%     Professor Anderson.
%
% Note:
%   If the variable PRESET is already defined out of this script, we will
%   use that value directly. This provides a way to set PRESET to any
%   needed value.
if ~exist('PRESET', 'var')
    PRESET = 'SigCapCSV20210928_Anderson';
end
if ismember(PRESET, ...
        {'SigCap_LongmontCampaign_20211105_Longmont', ...
        'SigCap_LongmontCampaign_20211105_LongmontExt'})
    folderNameForTestResults = 'SigCap_LongmontCampaign_20211105';
else
    folderNameForTestResults = PRESET;
end
% The absolute path to the folder with the FCC speed test app results.
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
    if ismember('lte_primary_rsrp', logTable.Properties.VariableNames)
        % Old version of SigCap has the field lte_primary_rsrp.
        for idxRecord = 1:size(logTable,1)
            logsStruct(end+1) = struct( ...
                'carrier', logTable.carrier{idxRecord}, ...
                'latitude', logTable.latitude(idxRecord), ...
                'longitude', logTable.longitude(idxRecord), ...
                'lte_primary_rsrp', logTable.lte_primary_rsrp(idxRecord) ...
                ); %#ok<SAGROW>
        end
    elseif ismember('primary/other*', logTable.Properties.VariableNames)
        % New version of SigCap has, instead, the field primary/other*.
        for idxRecord = 1:size(logTable,1)
            if strcmpi(logTable.("primary/other*"){idxRecord}, 'primary')
                logsStruct(end+1) = struct( ...
                    'carrier', logTable.carrier{idxRecord}, ...
                    'latitude', logTable.latitude(idxRecord), ...
                    'longitude', logTable.longitude(idxRecord), ...
                    'lte_primary_rsrp', logTable.rsrp(idxRecord) ...
                    ); %#ok<SAGROW>
            end
        end
    else
        error('Unsupported SigCap log file!');
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

if any(boolsIsAtt | boolsIsVer | boolsIsTM)
    hRsrpByCar = figure; hold on;
    if any(boolsIsAtt)
        hMAtt = stem3(allLons(boolsIsAtt), allLats(boolsIsAtt), ...
            -allPriRsrps(boolsIsAtt), markerForAtt);
    end
    if any(boolsIsVer)
        hMVer = stem3(allLons(boolsIsVer), allLats(boolsIsVer), ...
            -allPriRsrps(boolsIsVer), markerForVer);
    end
    if any(boolsIsTM)
        hMTM = stem3(allLons(boolsIsTM), allLats(boolsIsTM), ...
            -allPriRsrps(boolsIsTM), markerForTM);
    end

    plot_google_map('MapType', 'Hybrid');
    axis manual;

    if ~any(boolsIsAtt)
        hMAtt = stem3([], [], [], markerForAtt);
    end
    if ~any(boolsIsVer)
        hMVer = stem3([], [], [], markerForVer);
    end
    if ~any(boolsIsTM)
        hMTM = stem3([], [], [], markerForTM);
    end

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
end

if any(boolsIsAtt)
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
end

if any(boolsIsVer)
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
end

if any(boolsIsTM)
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
end

%% Cache Results

save(fullfile(pathToSaveResults, 'cache.mat'), ...
    'allLats', 'allLons', 'allPriRsrps', 'allCarriers', ...
    'boolsIsAtt', 'boolsIsVer', 'boolsIsTM');

% EOF