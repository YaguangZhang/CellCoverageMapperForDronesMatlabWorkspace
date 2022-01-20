%PLOTFIGSFORTWC2021 Generate some plots for publication after the
%simulation is done.
%
% This helper script is for the IEEE TWC paper prepared in 2021.
%
% Yaguang Zhang, Purdue, 10/29/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% Presets of interest.
PRESETS = {'Tipp', 'ShrinkedWHIN', 'ShrinkedIN'};
% Carrier frequencies of interest.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%     For broadband wireless backhaul.
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700, 7000, 13000, 28000};

% The LiDAR data set used in the simualtions.
LIDAR_DATA_SET_TO_USE = 'IN_DSM_2019';

% For plotting.
areaOfInterestColor = [0.9290 0.6940 0.1250];
lightBlue = [0.3010 0.7450 0.9330];
darkBlue = [0 0.4470 0.7410];
colorEffectiveTowers = 'b';
markerEffectiveTowers = '.';
markerSizeEffectiveTowers = 10;
colorIneffectiveTowers = 'r';
markerIneffectiveTowers = 'x';
lineWidthIneffectiveTowers = 1.5;

% The figure size was gotten by:
%   defaultFigPos = get(0,'defaultfigureposition');
%    customFigSize = defaultFigPos(3:4);
customFigSize = [560, 420];

% Path to save the plots.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'FigsForTwc2021');
% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% Turn the diary logging function on.
dirToDiaryFile = fullfile(pathToSaveResults, 'diary.txt');
% Get rid of the old diary if it exists.
if exist(dirToDiaryFile, 'file')
    delete(dirToDiaryFile);
end
diary(dirToDiaryFile);

% Reference boundaries in UTM.
pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
inBoundary = load(fullfile(pathToStateInfoFolder, ...
    'IN', 'boundary.mat'));
inBoundaryXYs = inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
UTM_ZONE = inBoundary.boundary.UTM_ZONE;
[~, whinBoundaryXYs, ~] ...
    = loadWhinBoundary(UTM_ZONE);

% For GPS and UTM conversions.
[deg2utm_speZone, utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(UTM_ZONE);

% Reference boundaries in (lat, lon).
[whinBoundaryLats, whinBoundaryLons] ...
    = utm2deg_speZone(whinBoundaryXYs(:,1), whinBoundaryXYs(:,2));
[inBoundaryLats, inBoundaryLons] ...
    = utm2deg_speZone(inBoundaryXYs(:,1), inBoundaryXYs(:,2));

%% Cell Towers to Consider on Roadmaps + User Location Grid

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating simulation overview plots ...'])

% We will use the 1900 MHz results.
numOfPresets = length(PRESETS);
freqInMhz = 1900;
for idxPreset = 1:numOfPresets
    preset = PRESETS{idxPreset};

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(freqInMhz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    load(fullfile(pathToReadResults, 'simConfigs.mat'));
    load(fullfile(pathToReadResults, 'simState.mat'));

    %% Simulation Area Overview
    % A plot for cell towers to consider in GPS (lon, lat).
    [effeCellAntsLats, effeCellAntsLons] ...
        = simConfigs.utm2deg_speZone( ...
        simState.CellAntsXyhEffective(:,1), ...
        simState.CellAntsXyhEffective(:,2));

    utmXYBoundaryPolyToKeepCellTowers ...
        = polybuffer(polyshape( ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
        simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);
    utmXYBoundaryToKeepCellTowers ...
        = utmXYBoundaryPolyToKeepCellTowers.Vertices;

    [gpsLatsBoundaryToKeepCellTowers, ...
        gpsLonsBoundaryToKeepCellTowers] ...
        = simConfigs.utm2deg_speZone( ...
        utmXYBoundaryToKeepCellTowers(:,1), ...
        utmXYBoundaryToKeepCellTowers(:,2));
    [gpsLatsBoundaryOfInterest, gpsLonsBoundaryOfInterest] ...
        = simConfigs.utm2deg_speZone( ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

    effeCellAntsBools = false(size(simState.CellAntsXyhAll, 1), 1);
    effeCellAntsBools(simState.CellAntsEffectiveIds) = true;
    inEffeCellAntsXYH = simState.CellAntsXyhAll(~effeCellAntsBools, :);
    [inEffeCellAntsLats, inEffeCellAntsLons] ...
        = simConfigs.utm2deg_speZone(inEffeCellAntsXYH(:,1), ...
        inEffeCellAntsXYH(:,2));

    % Resize the figure for this type of plots.
    curCustomFigSize = customFigSize.*0.9;

    hFigCellOverview = figure('Position', [0,0,curCustomFigSize]);
    hold on; set(gca, 'fontWeight', 'bold');
    hIneffeCells = plot(inEffeCellAntsLons, ...
        inEffeCellAntsLats, markerIneffectiveTowers, ...
        'Color', colorIneffectiveTowers, ...
        'LineWidth', lineWidthIneffectiveTowers);
    hExtendedArea = plot(polyshape([gpsLonsBoundaryToKeepCellTowers, ...
        gpsLatsBoundaryToKeepCellTowers; ...
        nan nan; ...
        gpsLonsBoundaryOfInterest, gpsLatsBoundaryOfInterest]));
    hAreaOfInterest = plot( ...
        polyshape([gpsLonsBoundaryOfInterest, gpsLatsBoundaryOfInterest]), ...
        'FaceColor', areaOfInterestColor);
    hEffeCells = plot(effeCellAntsLons, effeCellAntsLats, ...
        markerEffectiveTowers, ...
        'MarkerSize', markerSizeEffectiveTowers, ...
        'Color', colorEffectiveTowers);
    % Extend the content by a constant factor in the UTM system.
    extensionFactor = 0.2;
    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'ExtendedTipp')
        extensionFactor = 0.25;
    end
    [axisLonLatToSet, weightForWidth] ...
        = extendLonLatAxisByFactor( ...
        [min(gpsLonsBoundaryToKeepCellTowers), ...
        max(gpsLonsBoundaryToKeepCellTowers), ...
        min(gpsLatsBoundaryToKeepCellTowers), ...
        max(gpsLatsBoundaryToKeepCellTowers)], extensionFactor, simConfigs);
    adjustFigSizeByContent(hFigCellOverview, axisLonLatToSet, ...
        'height', weightForWidth.*0.9);
    view(2); plot_google_map;
    % grid on; grid minor;
    xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
    box on; xtickangle(0); xticks('manual');
    % Manually adjust the figure for publication.
    switch simConfigs.CURRENT_SIMULATION_TAG
        case 'Tipp'
            makescale('sw', 'units', 'si');
            hLeg = legend( ...
                [hAreaOfInterest, hExtendedArea, ...
                hEffeCells, hIneffeCells], ...
                'Area of interest', 'Extended area', ...
                'Cell tower to consider', ...
                'Ineffective cell tower', ...
                'defaultLegendAutoUpdate','off');
            set(hLeg, 'Position', [0.1375, 0.7417, 0.4718, 0.1839]);
        case 'ShrinkedWHIN'
            h = makescale(3.7, 'se', 'units', 'si');
            hPolyWhin = plot(whinBoundaryLons, whinBoundaryLats, ...
                ':', 'LineWidth', 2.3, 'Color', 'k');
            hLeg = legend(hPolyWhin, 'WHIN Boundary');
            set(hLeg, 'Position', [0.1686, 0.8792, 0.4669, 0.0476]);
        case 'ShrinkedIN'
            h = makescale(3.15, 'se', 'units', 'si');
            hPolyIn = plot(inBoundaryLons, inBoundaryLats, ...
                '-.', 'LineWidth', 2.3, 'Color', 'k');
            hLeg = legend(hPolyIn, 'Indiana Boundary');
            set(hLeg, 'Position', [0.1969, 0.8792, 0.5934, 0.0476]);
    end
    curDirToSave = fullfile(pathToSaveResults, ...
        ['Overview_CellularTowersToConsider_RoadMap-', preset]);
    saveEpsFigForPaper(hFigCellOverview, curDirToSave, false);

    %% User Location Grid

    % Resize the figure for this type of plots.
    curCustomFigSize = customFigSize.*0.9;

    hFigAreaOfInterest = figure('Position', [0,0,curCustomFigSize]);
    hCurAxis = gca;
    hold on; set(hCurAxis, 'fontWeight', 'bold');
    hAreaOfInterest = plot( ...
        polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
        'FaceColor', areaOfInterestColor);
    hGridPts = plot(simState.mapGridXYPts(:,1), ...
        simState.mapGridXYPts(:,2), '.', 'MarkerSize', 2.5, ...
        'Color', darkBlue);
    adjustFigSizeByContent(hFigAreaOfInterest, [], 'height', 0.9);
    axis equal; view(2); %grid on; grid minor;
    hLeg = legend([hAreaOfInterest, hGridPts], ...
        'Area of interest', 'User location');
    xlabel('UTM x (m)'); ylabel('UTM y (m)'); box on;
    % Adjust legend and the exponent label for y axis.
    switch simConfigs.CURRENT_SIMULATION_TAG
        case 'Tipp'
            set(hLeg, 'Position', [0.4789, 0.8062, 0.4258, 0.0966]);
        case 'ShrinkedWHIN'
            set(hLeg, 'Location', 'northwest');
            transparentizeCurLegends;

            % annotation(hFigAreaOfInterest, 'textbox',...
            %     [0.1609 0.8369 0.1785 0.0945],...
            %      'String', ['\times10^', ...
            %     num2str(hCurAxis.YAxis.Exponent)],...
            %      'FontWeight', hCurAxis.YAxis.FontWeight,...
            %     'FontSize', hCurAxis.YAxis.FontSize,...
            %      'EdgeColor', 'none');
            % yticks(yticks);
            %  yticklabels(yticklabels);

            set(hLeg,'visible','off');
        case 'ShrinkedIN'
            set(hLeg, 'Location', 'SouthEast');

            set(hLeg,'visible','off');
    end
    curDirToSave = fullfile(pathToSaveResults, ...
        ['Overview_UserLocGrid-', preset]);
    saveEpsFigForPaper(hFigAreaOfInterest, curDirToSave);
end

%% Cleanup

diary off;

% EOF