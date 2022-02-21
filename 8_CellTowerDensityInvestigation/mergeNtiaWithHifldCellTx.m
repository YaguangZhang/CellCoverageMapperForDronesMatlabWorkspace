% MERGENTIAWITHHIFLDCELLTX Merge the NTIA randomized U.S. laydown with the
% HIFLD Cell Tower locations.
%
% Yaguang Zhang, Purdue, 03/08/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Parameters
% The directories to load cellular tower location information. We have:
%   - National Telecommunications and Information Administration (NTIA)
%   randomized U.S. cellular laydown
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'RandomizedCarrierSitesv2.csv')
%   - Homeland Infrastructure Foundation-Level Data (HIFLD) Cellular Towers
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'HIFLD', ...
%            'CellTowers', 'Cellular_Towers_LatLonH.csv');
absPathToNtiaCsv = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'RandomizedCarrierSitesv2.csv');
absPathToHifldCsv = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'HIFLD', ...
    'CellTowers', 'Cellular_Towers_LatLonH.csv');

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs');

% Create directories if necessary.
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults);
end

% The directory to load the simulation results for Indiana State.
pathToLoadSimResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResultsBackups', 'In_Backup_Res_50_LowestAnts', ...
    'simResults.mat');

%% Load Data

disp(' ')
disp('    Loading cellular tower and area of interest information ...')

towerLatLonHsNtia = csvread(absPathToNtiaCsv, 1, 1); %#ok<*CSVRD>
towerLatLonHsHifld = csvread(absPathToHifldCsv, 1, 1);

load(pathToLoadSimResults);

% Functions to convert GPS degrees (lat, lon) from/to UTM (x, y).
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

% GPS boundary for the area of interest (here we have IN).
labelForAreaOfInterest = 'IN';
[boundOfInterestLats, boundOfInterestLons] ...
    = utm2deg_speZone(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

disp('    Done!')

%% Compute Distances

disp(' ')
disp('    Computing distMinInMToExistingT ...')

pathToSaveDistMin ...
    = fullfile(pathToSaveResults, 'distMinInMToExistingT.mat');

numOfExistingTs = size(towerLatLonHsNtia, 1);
numOfNewTs = size(towerLatLonHsHifld, 1);

areaOfIntLonLatPolyshape = ...
    polyshape(boundOfInterestLons, boundOfInterestLats);
boolsExistingTxInAreaOfInt = isinterior(areaOfIntLonLatPolyshape, ...
    towerLatLonHsNtia(:,2), towerLatLonHsNtia(:,1));
boolsNewTxInAreaOfInt = isinterior(areaOfIntLonLatPolyshape, ...
    towerLatLonHsHifld(:,2), towerLatLonHsHifld(:,1));
numOfExistingTsIndiana = sum(boolsExistingTxInAreaOfInt);
numOfNewTsIndiana = sum(boolsNewTxInAreaOfInt);

if exist(pathToSaveDistMin, 'file')
    disp('        Loading history results ...')
    load(pathToSaveDistMin);
else
    [distMinInMToExistingT, indicesToNearestExT] ...
        = deal(nan(numOfNewTs, 1));

    towerLatLonsNtia = towerLatLonHsNtia(:, 1:2);
    for idxNewT = 1:numOfNewTs
        curNewTLatLon = towerLatLonHsHifld(idxNewT, 1:2);

        curDistsInM = nan(numOfExistingTs, 1);

        parfor idxExistingT = 1:numOfExistingTs
            curDistsInM(idxExistingT) = lldistkm(curNewTLatLon, ...
                towerLatLonsNtia(idxExistingT, :)).*1000;
        end
        assert(~any(isnan(curDistsInM)), 'NaN distance value found!');

        [distMinInMToExistingT(idxNewT), indicesToNearestExT(idxNewT)] ...
            = min(curDistsInM);
    end

    save(pathToSaveDistMin, ...
        'distMinInMToExistingT', 'indicesToNearestExT', '-v7.3');
end
disp('    Done!')

%% Expirical CDF for distMinInMToExistingT

figResolutionFactor = 0.8;
hFigDistMinECDF = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,1600,1200].*figResolutionFactor);
hold on;
ecdf(distMinInMToExistingT);
[f, x] = ecdf(distMinInMToExistingT);
grid on; grid minor;
set(gca, 'XScale', 'log')

saveas(hFigDistMinECDF, ...
    fullfile(pathToSaveResults, 'distMinToExistingTEmpCdf.png'));
saveas(hFigDistMinECDF, ...
    fullfile(pathToSaveResults, 'distMinToExistingTEmpCdf.fig'));

%% Plot the Low Dist Locs on Map

distsToInspectInM = [50, 100, 250, 500, 1000, 1500, 1600];
numOfDistsToIns = length(distsToInspectInM);

figResolutionFactor = 1.25;
curAlpha = 0.5;
curDotSize = 10;

for idxDistToIns = 1:numOfDistsToIns
    curDistToInsInM = distsToInspectInM(idxDistToIns);
    boolsDistLessThanThreshold = distMinInMToExistingT<curDistToInsInM;

    curPathToSaveFig = fullfile(pathToSaveResults, ...
        ['lowDistNewLocs_NoTsWithDistsLessThan_', ...
        num2str(curDistToInsInM), 'm']);

    hFigLowDistNewLocs = figure('Visible', true, ...
        'Unit', 'pixel', 'Position', [0,0,1600,1200].*figResolutionFactor);
    hold on;

    hPolyIn = plot(boundOfInterestLons, boundOfInterestLats, ...
        'r-', 'LineWidth', 7);
    axisToSetIn = axis;
    set(hPolyIn, 'Visible', false);

    hNtia = plot(towerLatLonHsNtia(:,2), towerLatLonHsNtia(:,1), ...
        'b.', 'MarkerSize', curDotSize);
    axisToSet = axis;

    hNewLocs = scatter( ...
        towerLatLonHsHifld(~boolsDistLessThanThreshold,2), ...
        towerLatLonHsHifld(~boolsDistLessThanThreshold,1), ...
        curDotSize.*ones(sum(~boolsDistLessThanThreshold),1), ...
        'filled', 'MarkerFaceColor', 'g', 'Marker', 'o', ...
        'MarkerFaceAlpha', curAlpha);
    hLowDistNewLocs = plot( ...
        towerLatLonHsHifld(boolsDistLessThanThreshold,2), ...
        towerLatLonHsHifld(boolsDistLessThanThreshold,1), ...
        'c*', 'MarkerSize', curDotSize);
    hReplacedNtia = plot(towerLatLonHsNtia( ...
        indicesToNearestExT(boolsDistLessThanThreshold),2), ...
        towerLatLonHsNtia( ...
        indicesToNearestExT(boolsDistLessThanThreshold),1), ...
        'r.', 'MarkerSize', curDotSize);

    curNumOfLocsDisc = sum(boolsDistLessThanThreshold);
    title(['Threshold = ', num2str(curDistToInsInM), ' m (', ...
        num2str(curNumOfLocsDisc), '/', num2str(numOfNewTs), ...
        '~=', num2str(curNumOfLocsDisc/numOfNewTs*100), ...
        '% Locs Discarded)']);
    legend([hNtia, hReplacedNtia, hLowDistNewLocs, hNewLocs], ...
        'NTIA', 'Old Locs To Update', 'Updated New Locs', ...
        'Added New Locs');
    xticks([]); yticks([]); box on;

    plot_google_map('MapType', 'roadmap');

    % Overview plot for all data plotted.
    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_All.png']);
    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_All.fig']);

    % Overview plot over the area covered by the NTIA data set.
    axis(axisToSet);
    plot_google_map('MapType', 'roadmap');

    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_All.png']);
    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_All.fig']);

    % Overview plot for Indiana.
    curNumOfLocsDiscIndiana = ...
        sum(boolsNewTxInAreaOfInt&boolsDistLessThanThreshold);
    title(['Threshold = ', num2str(curDistToInsInM), ' m (', ...
        num2str(curNumOfLocsDiscIndiana), ...
        '/', num2str(numOfNewTsIndiana), ...
        '~=', num2str(curNumOfLocsDiscIndiana/numOfNewTsIndiana*100), ...
        '% IN Locs Discarded)']);
    legend([hPolyIn, hNtia, hReplacedNtia, hLowDistNewLocs, hNewLocs], ...
        'Indiana', 'NTIA', 'Old Locs To Update', 'Updated New Locs', ...
        'Added New Locs');

    set(hPolyIn, 'Visible', true); uistack(hPolyIn,'top');
    axis(axisToSetIn);
    plot_google_map('MapType', 'roadmap');

    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_IN.png']);
    saveas(hFigLowDistNewLocs, [curPathToSaveFig, '_IN.fig']);

    close all;
end

%% Export the Merged Results as a .csv File

THRESHOLD_DIST_IN_M = 1000;

dirCsvFileOut = fullfile(pathToSaveResults, ...
    ['NtiaLayoutMergedWithHifldCellTs_Threshold_', ...
    num2str(THRESHOLD_DIST_IN_M), 'm_LatLonH.csv']);

boolsDistLessThanThreshold = distMinInMToExistingT<THRESHOLD_DIST_IN_M;

lats = towerLatLonHsHifld(:,1);
lons = towerLatLonHsHifld(:,2);
hInM = towerLatLonHsHifld(:,3);

% Note: the nearest old loc could be the same for multiple new locs, so we
% may replace less old locs than expected.
towerLatLonHsNtiaNew = towerLatLonHsNtia;
towerLatLonHsNtiaNew( ...
    indicesToNearestExT(boolsDistLessThanThreshold), :) = [];
lats = [lats; towerLatLonHsNtiaNew(:,1)];
lons = [lons; towerLatLonHsNtiaNew(:,2)];
% Tower heights are not available in the NTIA data set.
hInM = [hInM; nan(length(towerLatLonHsNtiaNew(:,1)),1)];

% Sort by lat and lon.
latLonHInMs = sortrows([lats, lons, hInM],[1 2 3]);
numOfTowers = size(latLonHInMs, 1);

dataOut = array2table([(1:numOfTowers)', latLonHInMs]);
dataOut.Properties.VariableNames(1:4) ...
    = {'Site', 'Lat', 'Long', 'Height (m)'};
writetable(dataOut, dirCsvFileOut);

%% Overview Figure for Publication
% Cellular towers on road map.

hCellTowerOnRoadMap = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,1600,1200].*0.8);
hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
hPolyIn = plot(boundOfInterestLons, boundOfInterestLats, ...
    '-.', 'LineWidth', 5, 'Color', 'k'); % ones(1,3).*0.3
newAxis = extendAxisByFactor(axis, 0.1);
axis(newAxis);
tightfig;
plot_google_map('MapType', 'roadmap', 'Alpha', 0.5);
hNtiaNewRecords ...
    = plot(towerLatLonHsNtiaNew(:,2), towerLatLonHsNtiaNew(:,1), ...
    'o', 'Color', [0.225 0.675 1], 'MarkerSize', 5, 'LineWidth', 2);
hHifldCellTowers ...
    = plot(towerLatLonHsHifld(:,2), towerLatLonHsHifld(:,1), ...
    '.', 'Color', 'b', 'MarkerSize', 16);
plot( ...
    towerLatLonHsNtia(indicesToNearestExT( ...
    boolsDistLessThanThreshold),2), ...
    towerLatLonHsNtia(indicesToNearestExT( ...
    boolsDistLessThanThreshold),1), ...
    'wx', 'MarkerSize', 14, 'LineWidth', 3);
hNtiaIgnoredRecords ...
    = plot( ...
    towerLatLonHsNtia(indicesToNearestExT( ...
    boolsDistLessThanThreshold),2), ...
    towerLatLonHsNtia(indicesToNearestExT( ...
    boolsDistLessThanThreshold),1), ...
    'rx', 'MarkerSize', 14, 'LineWidth', 2);
uistack(hPolyIn,'top'); box on;
hLeg = legend([hHifldCellTowers, hNtiaNewRecords, hNtiaIgnoredRecords, ...
    hPolyIn], ...
    'HIFLD cell tower', 'New tower from NTIA', 'Ignored NTIA record', ...
    'Indiana boundary', ...
    'FontSize', 25, 'Location', 'southeast', 'LineWidth', 3);
% Move IN state to the center.
moveIndianaStateToCenter;

% Further tighten the figure.
%   tightfig(hCellTowerOnRoadMap); % No Effect...

% Manually adjust the legend position.
set(hLeg, 'Position', [0.6235, 0.0109, 0.3689, 0.2185]);

saveas(hCellTowerOnRoadMap, ...
    fullfile(pathToSaveResults, 'cellTowerOnRoadMap.png'));
saveas(hCellTowerOnRoadMap, ...
    fullfile(pathToSaveResults, 'cellTowerOnRoadMap.fig'));
saveas(hCellTowerOnRoadMap, ...
    fullfile(pathToSaveResults, 'cellTowerOnRoadMap.eps'), 'epsc');

% EOF