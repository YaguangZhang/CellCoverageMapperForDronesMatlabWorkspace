%TESTFCCSPEEDTESTAPPRESULTS A script to read in and visualize results from
%the FCC speed test app.
%
% Yaguang Zhang, Purdue, 10/14/2021
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
%   - FCC_SpeedTest_Anderson_20211007
%     Four tests from Prof. Anderson.
%   - FCC_SpeedTest_Yaguang_20211014
%     Two tests from Yaguang.
%   - FCC_SpeedTest_Anderson_20211014
%     Raw .zip files from four students.
%   - FCC_SpeedTest_LongmontCampaign_20211105
%     Copies of the Longmont measurement campaign results.
%   - FCC_SpeedTest_LongmontCampaign_20211105_USNA
%     Only show tests conducted in USNA.
%   - FCC_SpeedTest_LongmontCampaign_20211105_Longmont
%     Only show tests conducted in Longmont.
%   - FCC_SpeedTest_LongmontCampaign_20211105_LongmontExt
%     Only show tests conducted in and around Longmont.
if ~exist('PRESET', 'var')
    PRESET = 'FCC_SpeedTest_LongmontCampaign_20211105_LongmontExt';
end
if ismember(PRESET, ...
        {'FCC_SpeedTest_LongmontCampaign_20211105_Longmont', ...
        'FCC_SpeedTest_LongmontCampaign_20211105_LongmontExt', ...
        'FCC_SpeedTest_LongmontCampaign_20211105_USNA'})
    folderNameForTestResults = 'FCC_SpeedTest_LongmontCampaign_20211105';
else
    folderNameForTestResults = PRESET;
end
pathToFccSpeedTestAppResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'FCC Rural Broadband', folderNameForTestResults);

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '9_SimulationsForCO', ...
    'FccSpeedTestAppTests', PRESET);
% Create directories if necessary.
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults);
end

% When specified, "latLonsForAreaOfInterest" will be used to filter out
% tests based on their location information.
switch PRESET
    case 'FCC_SpeedTest_LongmontCampaign_20211105_Longmont'
        latLonsForAreaOfInterest = ptPairToRect([ ...
            40.062060015391786, -105.25808569083334, ...
            40.277115808065005, -104.94420519323936]);
    case 'FCC_SpeedTest_LongmontCampaign_20211105_LongmontExt'
        latLonsForAreaOfInterest = ptPairToRect([ ...
            39.80668743376393, -105.38602021079713, ...
            40.40409728561462, -104.68631223445952]);
    case 'FCC_SpeedTest_LongmontCampaign_20211105_USNA'
        latLonsForAreaOfInterest = ptPairToRect([ ...
            38.96013809274135, -76.51826130542474, ...
            38.99854785776838, -76.4573199249457]);
end

%% Read All the Data

logsJson = rdir(fullfile(pathToFccSpeedTestAppResults, '**', '*.json'));
logsZip = rdir(fullfile(pathToFccSpeedTestAppResults, '**', '*.zip'));
if isempty(logsJson) || ...
        ((~isempty(logsJson)) && (length(logsJson)<length(logsZip)))
    % Unzip the raw .zip files first.
    for idxLogZip = 1:length(logsZip)
        unzip(logsZip(idxLogZip).name, logsZip(idxLogZip).name(1:(end-4)));
    end

    % We need to get rid of the logs under the __MACOSX folders.
    invalidLogs = rdir(fullfile(pathToFccSpeedTestAppResults, ...
        '**', '__MACOSX', '**', '*.json'));
    for idxInvalidLog = 1:length(invalidLogs)
        delete(invalidLogs(idxInvalidLog).name);
    end

    % Reload the files.
    logsJson = rdir(fullfile(pathToFccSpeedTestAppResults, '**', '*.json'));

    % This check does not always apply. For example, when the app exports
    % multiple tests in one attemp.
    %
    % logsZip = rdir(fullfile(pathToFccSpeedTestAppResults, ...
    %     '**', '*.zip'));
    %
    % assert(length(logsJson)==length(logsZip), ...
    %     ['Error! The number of raw .zip files (', ...
    %      num2str(length(logsZip)), ...
    %     ') does not agree with the number of JSON files (', ...
    %      num2str(length(logsJson)), ')!']);
end
assert(~isempty(logsJson), 'Error! No JSON log files found!')

logsStructAndroid = struct('version', {}, ...
    'carrier_name', {}, 'iso_country_code', {}, 'manufacturer', {}, ...
    'mobile_country_code', {}, 'mobile_network_code', {}, 'model', {}, ...
    'operating_system_version', {}, 'tests', {}, 'submission_type', {});
logsStructIos = struct('device_environment', {}, ...
    'tests', {}, 'metadata', {});
numOfLogs = length(logsJson);
[indicesAndroidLogs, indicesIosLogs] = deal([]);
for idxLog = 1:numOfLogs
    % Read files as text.
    logStr = fileread(logsJson(idxLog).name);
    % Examine the JSON object and adjust the structure if necessary.
    logNew = jsondecode(logStr);
    try
        logsStructAndroid(end+1) = logNew; %#ok<SAGROW>
        indicesAndroidLogs = [indicesAndroidLogs; idxLog]; %#ok<AGROW>
        assert(strcmp(logNew.operating_system_version(1:7), 'Android'), ...
            'Error! Input log is generated by ', ...
            logNew.operating_system_version, ' instead of Android!')
    catch
        logsStructIos(end+1) = logNew; %#ok<SAGROW>
        indicesIosLogs = [indicesIosLogs; idxLog]; %#ok<AGROW>
        osVersion = logNew.device_environment.operating_system_version;
        assert(strcmp(osVersion(1:3), 'iOS'), ...
            'Error! Input log is generated by ', osVersion, ...
            ' instead of iSO!')
    end
end

%% Extract the Information of Interest
% Note that we will reorder the logs by the OS name & version.

[osVersions, carriers, testIds] = deal(cell(numOfLogs, 1));
[downSpeedsBps, upSpeedsBps] = deal(nan(numOfLogs, 1));
[downEndLatLons, upEndLatLons] = deal(nan(numOfLogs, 2));
flagsIsCell = false(numOfLogs, 1);

numOfAndroidLogs = length(logsStructAndroid);
% indicesIgnoredLog = [];
for idxLog = 1:numOfLogs
    try
        if idxLog<=numOfAndroidLogs
            [osVersions{idxLog}, carriers{idxLog}, ...
                downSpeedsBps(idxLog), upSpeedsBps(idxLog), ...
                downEndLatLons(idxLog, :), upEndLatLons(idxLog, :), ...
                flagsIsCell(idxLog), testIds{idxLog}] ...
                = parseFccSpeedTestAppLog(logsStructAndroid(idxLog));
        else
            [osVersions{idxLog}, carriers{idxLog}, ...
                downSpeedsBps(idxLog), upSpeedsBps(idxLog), ...
                downEndLatLons(idxLog, :), upEndLatLons(idxLog, :), ...
                flagsIsCell(idxLog), testIds{idxLog}] ...
                = parseFccSpeedTestAppLog(logsStructIos(...
                idxLog-numOfAndroidLogs));
        end

        % If we failed to fetch the test ID from the log file, the file
        % name will be used, instead.
        if isempty(testIds{idxLog})
            if idxLog<=numOfAndroidLogs
                curLogDir = logsJson(indicesAndroidLogs(idxLog));
            else
                curLogDir = logsJson(indicesIosLogs( ...
                    idxLog-numOfAndroidLogs));
            end
            [~, testIds{idxLog}] = fileparts(curLogDir.name);
        end
    catch err
        % indicesIgnoredLog(end+1) = idxLog; %#ok<SAGROW>
        warning(['Ignoring invalid log #', num2str(idxLog), '!']);
        if idxLog<=numOfAndroidLogs
            disp(logsJson(indicesAndroidLogs(idxLog)));
        else
            disp(logsJson(indicesIosLogs(idxLog-numOfAndroidLogs)));
        end
        disp(err);
    end
end

% Remove invalid cellular speed test results.

% indicesLogsWithUnknownOs = find(cellfun(@isempty, osVersions))';
%  for idxLogWithUnkownOs = indicesLogsWithUnknownOs
%     osVersions{idxLogWithUnkownOs} = 'Unkown';
%  end

% If necessary, remove tests out of the area of interest.
if exist('latLonsForAreaOfInterest', 'var')
    boolsToKeep = inpolygon(downEndLatLons(:,1), downEndLatLons(:,2), ...
        latLonsForAreaOfInterest(:,1), latLonsForAreaOfInterest(:,2)) ...
        & inpolygon(upEndLatLons(:,1), upEndLatLons(:,2), ...
        latLonsForAreaOfInterest(:,1), latLonsForAreaOfInterest(:,2));
    flagsEleToRemove = ~boolsToKeep; %#ok<NASGU>
    cleanFccLogVars;

    flagsIsCell = flagsIsCell(boolsToKeep);
    clearvars latLonsForAreaOfInterest;
end

% Output raw results to a .csv file.
tableHeader = {'osVersion','carrier', ...
    'downEndLatLons','downSpeedsBps', ...
    'upEndLatLons','upSpeedsBps', ...
    'testId'};
rawTestResults = cell2table([osVersions, carriers, ...
    num2cell(downEndLatLons, 2), num2cell(downSpeedsBps), ...
    num2cell(upEndLatLons, 2), num2cell(upSpeedsBps), testIds]);
rawTestResults.Properties.VariableNames = tableHeader;
writetable(rawTestResults, ...
    fullfile(pathToSaveResults, 'validTestResults.csv'));

flagsEleToRemove = ~flagsIsCell;
cleanFccLogVars;

% Output cellular results to a .csv file.
rawTestResults = cell2table([osVersions, carriers, ...
    num2cell(downEndLatLons, 2), num2cell(downSpeedsBps), ...
    num2cell(upEndLatLons, 2), num2cell(upSpeedsBps), testIds]);
rawTestResults.Properties.VariableNames = tableHeader;
writetable(rawTestResults, ...
    fullfile(pathToSaveResults, 'validCellTestResults.csv'));

%% Plots

% Down/up speed on map.
hDownSpeed = figure; hold on;
plot3k([downEndLatLons(:, [2, 1]), downSpeedsBps./(10^6)], ...
    'Labels', {'Download Speed on Map', '', '', '', ...
    'Data Rate (Mbps)'}, ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hDownSpeed, ...
    fullfile(pathToSaveResults, 'downloadSpeed_3D.jpg'));
view(2);
saveas(hDownSpeed, ...
    fullfile(pathToSaveResults, 'downloadSpeed_2D.jpg'));
saveas(hDownSpeed, ...
    fullfile(pathToSaveResults, 'downloadSpeed_2D.fig'));

hUpSpeed = figure; hold on;
plot3k([downEndLatLons(:, [2, 1]), upSpeedsBps./(10^6)], ...
    'Labels', {'Upload Speed on Map', '', '', '', ...
    'Data Rate (Mbps)'}, ...
    'Plottype','stem', 'Marker', {'.', 12});
xticklabels([]); yticklabels([]); view(3);
plot_google_map('MapType', 'Hybrid');

saveas(hUpSpeed, ...
    fullfile(pathToSaveResults, 'upSpeed_3D.jpg'));
view(2);
saveas(hUpSpeed, ...
    fullfile(pathToSaveResults, 'upSpeed_2D.jpg'));
saveas(hUpSpeed, ...
    fullfile(pathToSaveResults, 'upSpeed_2D.fig'));

% Down/up speed by OS.
boolsIsAndroid = cellfun(@(v) strcmp(v(1:3), 'And'), osVersions);
markerForAndroid = ':.r';
markerForIos = ':.b';

hDownSpeedByOs = figure; hold on;
hMAnd = stem3(downEndLatLons(boolsIsAndroid, 2), ...
    downEndLatLons(boolsIsAndroid, 1), ...
    downSpeedsBps(boolsIsAndroid)./(10^6), markerForAndroid);
hMIos = stem3(downEndLatLons(~boolsIsAndroid, 2), ...
    downEndLatLons(~boolsIsAndroid, 1), ...
    downSpeedsBps(~boolsIsAndroid)./(10^6), markerForIos);
title('Download Speed by Operating System');
zlabel('Data Rate (Mbps)');
xticklabels([]); yticklabels([]); view(3);
legend([hMAnd, hMIos], 'Android', 'iOS', 'Location','best');
plot_google_map('MapType', 'Hybrid');

saveas(hDownSpeedByOs, ...
    fullfile(pathToSaveResults, 'downloadSpeedByOs_3D.jpg'));
view(2); legend([hMAnd, hMIos], 'Android', 'iOS', 'Location','best');
saveas(hDownSpeedByOs, ...
    fullfile(pathToSaveResults, 'downloadSpeedByOs_2D.jpg'));
saveas(hDownSpeedByOs, ...
    fullfile(pathToSaveResults, 'downloadSpeedByOs_2D.fig'));

hUpSpeedByOs = figure; hold on;
hMAnd = stem3(upEndLatLons(boolsIsAndroid, 2), ...
    upEndLatLons(boolsIsAndroid, 1), ...
    upSpeedsBps(boolsIsAndroid)./(10^6), markerForAndroid);
hMIos = stem3(upEndLatLons(~boolsIsAndroid, 2), ...
    upEndLatLons(~boolsIsAndroid, 1), ...
    upSpeedsBps(~boolsIsAndroid)./(10^6), markerForIos);
title('Upload Speed by Operating System');
zlabel('Data Rate (Mbps)');
xticklabels([]); yticklabels([]); view(3);
legend([hMAnd, hMIos], 'Android', 'iOS', 'Location','best');
plot_google_map('MapType', 'Hybrid');

saveas(hUpSpeedByOs, ...
    fullfile(pathToSaveResults, 'uploadSpeedByOs_3D.jpg'));
view(2); legend([hMAnd, hMIos], 'Android', 'iOS', 'Location','best');
saveas(hUpSpeedByOs, ...
    fullfile(pathToSaveResults, 'uploadSpeedByOs_2D.jpg'));
saveas(hUpSpeedByOs, ...
    fullfile(pathToSaveResults, 'uploadSpeedByOs_2D.fig'));

% Down/up speed by carrier.
boolsIsAtt = cellfun(@(v) strcmp(v(1:3), 'AT&'), carriers);
boolsIsVer = cellfun(@(v) strcmp(v(1:3), 'Ver'), carriers);
boolsIsTM = cellfun(@(v) strcmp(v(1:3), 'T-M'), carriers);
boolsNoSig = cellfun(@(v) strcmp(v(1:3), 'no_'), carriers);
assert(all(boolsIsAtt|boolsIsVer|boolsIsTM|boolsNoSig), ...
    'Unknown carrier label(s)!');

markerForAtt = ':r.';
markerForVer = ':b.';
markerForTM = ':g.';
markerForNoSig = ':xk';

hDownSpeedByCar = figure; hold on;
hMAtt = stem3(downEndLatLons(boolsIsAtt, 2), ...
    downEndLatLons(boolsIsAtt, 1), ...
    downSpeedsBps(boolsIsAtt)./(10^6), markerForAtt);
hMVer = stem3(downEndLatLons(boolsIsVer, 2), ...
    downEndLatLons(boolsIsVer, 1), ...
    downSpeedsBps(boolsIsVer)./(10^6), markerForVer);
hMTM = stem3(downEndLatLons(boolsIsTM, 2), ...
    downEndLatLons(boolsIsTM, 1), ...
    downSpeedsBps(boolsIsTM)./(10^6), markerForTM);
plot_google_map('MapType', 'Hybrid');
hMNoSig = stem3(downEndLatLons(boolsNoSig, 2), ...
    downEndLatLons(boolsNoSig, 1), ...
    downSpeedsBps(boolsNoSig)./(10^6), markerForNoSig);
title('Download Speed by Carriers');
zlabel('Data Rate (Mbps)');
xticklabels([]); yticklabels([]); view(3);
legend([hMAtt, hMVer, hMTM, hMNoSig], 'AT&T', 'Verizon', 'T-Mobile', ...
    'No Service', 'Location','best');

saveas(hDownSpeedByCar, ...
    fullfile(pathToSaveResults, 'downloadSpeedByCar_3D.jpg'));
view(2);
legend([hMAtt, hMVer, hMTM, hMNoSig], 'AT&T', 'Verizon', 'T-Mobile', ...
    'No Service', 'Location','best');
saveas(hDownSpeedByCar, ...
    fullfile(pathToSaveResults, 'downloadSpeedByCar_2D.jpg'));
saveas(hDownSpeedByCar, ...
    fullfile(pathToSaveResults, 'downloadSpeedByCar_2D.fig'));

hUpSpeedByCar = figure; hold on;
hMAtt = stem3(upEndLatLons(boolsIsAtt, 2), ...
    upEndLatLons(boolsIsAtt, 1), ...
    upSpeedsBps(boolsIsAtt)./(10^6), markerForAtt);
hMVer = stem3(upEndLatLons(boolsIsVer, 2), ...
    upEndLatLons(boolsIsVer, 1), ...
    upSpeedsBps(boolsIsVer)./(10^6), markerForVer);
hMTM = stem3(upEndLatLons(boolsIsTM, 2), ...
    upEndLatLons(boolsIsTM, 1), ...
    upSpeedsBps(boolsIsTM)./(10^6), markerForTM);
plot_google_map('MapType', 'Hybrid');
hMNoSig = stem3(upEndLatLons(boolsNoSig, 2), ...
    upEndLatLons(boolsNoSig, 1), ...
    upSpeedsBps(boolsNoSig)./(10^6), markerForNoSig);
title('Upload Speed by Carriers');
zlabel('Data Rate (Mbps)');
xticklabels([]); yticklabels([]); view(3);
legend([hMAtt, hMVer, hMTM, hMNoSig], 'AT&T', 'Verizon', 'T-Mobile', ...
    'No Service', 'Location','best');

saveas(hUpSpeedByCar, ...
    fullfile(pathToSaveResults, 'uploadSpeedByCar_3D.jpg'));
view(2);
legend([hMAtt, hMVer, hMTM, hMNoSig], 'AT&T', 'Verizon', 'T-Mobile', ...
    'No Service', 'Location','best');
saveas(hUpSpeedByCar, ...
    fullfile(pathToSaveResults, 'uploadSpeedByCar_2D.jpg'));
saveas(hUpSpeedByCar, ...
    fullfile(pathToSaveResults, 'uploadSpeedByCar_2D.fig'));
% EOF