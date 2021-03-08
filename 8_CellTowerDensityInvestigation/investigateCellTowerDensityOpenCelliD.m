%INVESTIGATECELLTOWERDENSITY Visualize the cellular towers and explore
%their geographic density.
%
% Yaguang Zhang, Purdue, 01/06/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Script Parameters

% The directories to load cellular tower location information. We have:
%   - OpenCelliD
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'OpenCelliDUsa_20210219', ...
%            'Cellular_Towers_LatLonHR.csv');
%   - NTIA randomized U.S. cellular laydown
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'RandomizedCarrierSitesv2.csv')
%   - Homeland Infrastructure Foundation-Level Data (HIFLD) Cellular Towers
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'HIFLD', ...
%            'CellTowers', 'Cellular_Towers_LatLon.csv');
%   - HIFLD Land Mobile Commercial Transmission Towers
%       fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
%           'CellTowerInfo', 'HIFLD', ...
%            'LandMobileCommercialTxTowers', ...
%           'Land_Mobile_Commercial_Transmission_Towers_LatLonH.csv');
ABS_PATH_TO_CELL_ANTENNAS_CSV =  fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'OpenCelliDUsa_20210219', ...
    'Cellular_Towers_LatLonHR.csv');
% The values to show contour lines. For the NTIA randomized cell tower
% layout, it is recommended to use:
%       contourLevelVs = [1, 3, 5, 10, 15, 30, 50, 80, 100, 150];
% For the HIFLD Land Mobile Commercial Transmission Towers:
%       contourLevelVs = [2, 3, 5, 7, 9, 15, 20, 50, 80, 110];
contourLevelVs = [150, 200, 300, 500, 800, 1300, 2000, 3000, 5000];

% Used to adjust the exported resolution via the figure position attribute.
figResolutionFactor = 0.8;
% The position for the legend in the tower density plots. Recommendations:
%   - [0.6,0.01,0.37,0.1] for figResolutionFactor = 1;
%    - [0.53, 0.01, 0.46, 0.1] for figResolutionFactor = 0.8;
towerDensityLegendPos = [0.53, 0.01, 0.46, 0.1];

% The radius to inspect around a point of interest for estimating its
% cellular tower density. Based on:
%   https://en.wikipedia.org/wiki/Cell_site
% for 3G/4G/5G (FR1) Mobile base station tower: it is technically possible
% to cover up to 50 km-150 km (Macrocell). However, according to:
%   https://amphenolwireless.com/macro/
% a typical macrocell covers a distance from 1 km to 30 km. We will count
% the number of towers in the range (distance <= radius to inspect) for the
% cellular tower density computation.
RADIUS_TO_INSPECT_IN_M = 30000;

% The directory to load the simulation results for Indiana State.
pathToLoadSimResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResultsBackups', 'In_Backup_Res_50_LowestAnts', ...
    'simResults.mat');

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '8_CellTowerDensityInvestigation_OpenCelliD');

%% Initialization

disp(' ')
disp('    Loading history simulation results ...')

load(pathToLoadSimResults);

% Functions to convert GPS degrees (lat, lon) from/to UTM (x, y). We will
% populate them later.

% Convert GPS degrees to UTM coordinates for the specified zone.
utmstruct_speZone = defaultm('utm');
% Remove white space in the zone label.
utmstruct_speZone.zone ...
    = simConfigs.UTM_ZONE(~isspace(simConfigs.UTM_ZONE));
utmstruct_speZone.geoid = wgs84Ellipsoid;
utmstruct_speZone = defaultm(utmstruct_speZone);

deg2utm_speZone = @(lat, lon) mfwdtran(utmstruct_speZone, lat,lon);
utm2deg_speZone = @(x, y) minvtran(utmstruct_speZone, x, y);

% GPS boundary for the area of interest (here we have IN).
[boundOfInterestLats, boundOfInterestLons] ...
    = utm2deg_speZone(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

disp('    Done!')

%% Load Cellular Tower Information

disp(' ')
disp('    Loading cellular antenna information ...')

% Note: we use "height" to indicate the vertical distance from the ground
% to the antenna; "elevation" to indicate the ground elevation; and
% "altitude" to indicate elevation+height.
cellAntsLatLonAlt = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);
numAnts = size(cellAntsLatLonAlt, 1);
cellAntsXYH = nan(numAnts, 3);
[cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
    = deg2utm_speZone(cellAntsLatLonAlt(:,1), cellAntsLatLonAlt(:,2));
if size(cellAntsLatLonAlt, 2)>2
    cellAntsXYH(:,3) = cellAntsLatLonAlt(:,3);
end

disp('    Done!')

%% Estimate Cellular Tower Density for IN

disp(' ')
disp('    Estimating cellular tower density for IN ...')

numGridPts = size(simState.mapGridXYPts, 1);
cellTowerDensitiesNumPer1000SqKms = nan(numGridPts,1);
areaToInspectInSqKm = pi*(RADIUS_TO_INSPECT_IN_M/1000)^2;
for idxGridPt = 1:numGridPts
    curGridXY = simState.mapGridXYPts(idxGridPt,:);
    cellTowerDists = pdist2(curGridXY, cellAntsXYH(:,1:2), 'euclidean');
    cellTowerDensitiesNumPer1000SqKms(idxGridPt) = ...
        sum(cellTowerDists<=RADIUS_TO_INSPECT_IN_M) ...
        /(areaToInspectInSqKm/1000);
end

disp('    Done!')

%% Plots

% Cellular tower on road map for US.
hCellTowerOnRoadMap = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,1600,1200].*figResolutionFactor);
hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
hPolyIn = plot(boundOfInterestLons, boundOfInterestLats, ...
    'r-', 'LineWidth', 7);
newAxis = extendAxisByFactor(axis, 0.1);
axis(newAxis);
tightfig;
plot_google_map('MapType', 'roadmap');
hCellTowers = plot(cellAntsLatLonAlt(:,2), cellAntsLatLonAlt(:,1), ...
    'b.', 'MarkerSize', 5);
uistack(hPolyIn,'top'); box on;
legend([hCellTowers, hPolyIn], 'Cell Tower', 'Indiana', ...
    'FontSize', 25, 'Location', 'southeast', 'LineWidth', 3);
% Move IN state to the center.
moveIndianaStateToCenter;
saveas(hCellTowerOnRoadMap, ...
    fullfile(pathToSaveResults, 'cellTowerOnRoadMap.png'));
saveas(hCellTowerOnRoadMap, ...
    fullfile(pathToSaveResults, 'cellTowerOnRoadMap.fig'));

% Cellular tower density for IN.
hCellTowerDensity = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,900,1200].*figResolutionFactor);
hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
plot3(boundOfInterestLons, boundOfInterestLats, ...
    ones(size(boundOfInterestLons)) ...
    .*(max(cellTowerDensitiesNumPer1000SqKms)+2), ...
    'r-', 'LineWidth', 7);
newAxis = extendAxisByFactor(axis, -0.1);
axis(newAxis);
tightfig;
plot_google_map('MapType', 'roadmap');
moveIndianaStateToCenter;
% plot3k([simState.mapGridLatLonPts(:, 2:-1:1) ...
%     cellTowerDensitiesNumPerSqKm]);
hCellTowers = plot(cellAntsLatLonAlt(:,2), cellAntsLatLonAlt(:,1), ...
    'b.', 'MarkerSize', 5);
surfOnSimGrid([simState.mapGridLatLonPts(:, 2:-1:1) ...
    cellTowerDensitiesNumPer1000SqKms], simConfigs, ...
    contourLevelVs, 0.75, 22, 3, 'k');
colormap(colormap(flipud(hot))); set(gca, 'ColorScale', 'log'); box on;
hAT = autoText('Cell Tower Density (# per 1000 km^2)', ...
    'Location','southeast', 'FontSize', 22, 'LineWidth', 3);
set(hAT, 'FitBoxToText', 'on', 'LineStyle', '-', 'FontSize', 22, ...
    'LineWidth', 3, 'BackGroundColor', 'w', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Position', towerDensityLegendPos);
% Move IN state to the center.
moveIndianaStateToCenter;
saveas(hCellTowerDensity, ...
    fullfile(pathToSaveResults, 'cellTowerDensity.png'));
saveas(hCellTowerDensity, ...
    fullfile(pathToSaveResults, 'cellTowerDensity.fig'));

% Cellular tower density (no tower locs) for IN.
hCellTowerDensityOnly = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,900,1200].*figResolutionFactor);
hold on; xticks([]); xticklabels({}); yticks([]); yticklabels({});
hPolyIn = plot3(boundOfInterestLons, boundOfInterestLats, ...
    ones(size(boundOfInterestLons)) ...
    .*(max(cellTowerDensitiesNumPer1000SqKms)+2), ...
    'r-', 'LineWidth', 7);
newAxis = extendAxisByFactor(axis, -0.1);
axis(newAxis);
tightfig;
plot_google_map('MapType', 'roadmap');
moveIndianaStateToCenter;
% plot3k([simState.mapGridLatLonPts(:, 2:-1:1) ...
%     cellTowerDensitiesNumPerSqKm]);
% hCellTowers = plot(cellAntsLatLonAlt(:,2), cellAntsLatLonAlt(:,1), ...
%     'b.', 'MarkerSize', 5);
hSurf = surfOnSimGrid([simState.mapGridLatLonPts(:, 2:-1:1) ...
    cellTowerDensitiesNumPer1000SqKms], simConfigs, ...
    contourLevelVs, 0.75, 22, 3, 'k');
colormap(colormap(flipud(hot))); set(gca, 'ColorScale', 'log'); box on;
hAT = autoText('Cell Tower Density (# per 1000 km^2)', ...
    'Location','southeast', 'FontSize', 22, 'LineWidth', 3);
set(hAT, 'FitBoxToText', 'on', 'LineStyle', '-', 'FontSize', 22, ...
    'LineWidth', 3, 'BackGroundColor', 'w', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'Position', towerDensityLegendPos);
saveas(hCellTowerDensityOnly, ...
    fullfile(pathToSaveResults, 'cellTowerDensityOnly.png'));
saveas(hCellTowerDensityOnly, ...
    fullfile(pathToSaveResults, 'cellTowerDensityOnly.fig'));

%% Empirical CDF for cellular tower density of IN

hCellTowerDensityECDF = figure('Visible', true, ...
    'Unit', 'pixel', 'Position', [0,0,1600,1200].*figResolutionFactor);
hold on;
ecdf(cellTowerDensitiesNumPer1000SqKms);
[f, x] = ecdf(cellTowerDensitiesNumPer1000SqKms);
contourLevelVs = interp1(f(2:end), x(2:end), 0:0.1:1);
F_300 = interp1(x(2:end), f(2:end), 300);
F_600 = interp1(x(2:end), f(2:end), 600);
title({['F(300) = ', num2str(F_300), '; F(600) = ', num2str(F_600)]; ...
    ['Ref contourLevelVs = ', num2str(contourLevelVs)]});
grid on; grid minor;
saveas(hCellTowerDensityECDF, ...
    fullfile(pathToSaveResults, 'cellTowerDensityECDF.png'));
saveas(hCellTowerDensityECDF, ...
    fullfile(pathToSaveResults, 'cellTowerDensityECDF.fig'));

% EOF