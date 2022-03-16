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
% CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700, 7000, 13000, 28000};
CARRIER_FREQUENCIES_IN_MHZ = {1900};

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

% Orange - red - black.
customHot = hot(256); customHot = customHot(168:-1:1, :);
% Yellow - red - black.
customHotLong = hot(256); customHotLong = customHotLong(192:-1:1, :);

% For scale legend adjustment.
largeZValue = 10^7;
tippScalePos = [-86.7638, 40.235, largeZValue];
tippScaleFontS = 9;
whinScalePos = [-86.3616, 39.945, largeZValue];
whinScaleFontS = 9;
inScalePos = [-85.5442, 38.3774, largeZValue];
inScaleFontS = 9;

%% Cell Towers to Consider on Roadmaps + User Location Grid

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating simulation overview plots ...'])

% The figure size was gotten by:
%   defaultFigPos = get(0,'defaultfigureposition');
%    customFigSize = defaultFigPos(3:4);
curCustomFigSize = [560, 420]*0.9;

% We will use the 1900 MHz results.
numOfPresets = length(PRESETS);
freqInMhz = 1900;
for idxPreset = 1:numOfPresets
    preset = PRESETS{idxPreset};
    disp(['        [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(freqInMhz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    clearvars simConfigs simState;
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
            hLeg = legend(hPolyWhin, 'WHIN boundary');
            set(hLeg, 'Position', [0.1686, 0.8792, 0.4669, 0.0476]);
        case 'ShrinkedIN'
            h = makescale(3.15, 'se', 'units', 'si');
            hPolyIn = plot(inBoundaryLons, inBoundaryLats, ...
                '-.', 'LineWidth', 2.3, 'Color', 'k');
            hLeg = legend(hPolyIn, 'Indiana boundary');
            set(hLeg, 'Position', [0.1969, 0.8792, 0.5934, 0.0476]);
    end

    curDirToSave = fullfile(pathToSaveResults, ...
        ['Overview_CellularTowersToConsider_RoadMap-', preset]);
    saveEpsFigForPaper(hFigCellOverview, curDirToSave, false);

    %% User Location Grid

    hFigAreaOfInterest = figure('Position', [0,0,curCustomFigSize]);
    hCurAxis = gca;
    hold on; set(hCurAxis, 'fontWeight', 'bold');
    hAreaOfInterest = plot( ...
        polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
        'FaceColor', areaOfInterestColor);
    hGridPts = plot(simState.mapGridXYPts(:,1), ...
        simState.mapGridXYPts(:,2), '.', 'MarkerSize', 3, ...
        'Color', darkBlue);
    % axis equal;
    adjustFigSizeByContent(hFigAreaOfInterest, [], 'height', 0.9);
    axis equal; view(2); %grid on; grid minor;
    hLeg = legend([hAreaOfInterest, hGridPts], ...
        'Area of interest', 'User location');
    xlabel('UTM x (m)'); ylabel('UTM y (m)'); box on;
    % Tighten the figure twice.
    %   tightfig(hFigAreaOfInterest); tightfig(hFigAreaOfInterest);
    % Adjust legend and the exponent label for y axis.
    switch simConfigs.CURRENT_SIMULATION_TAG
        case 'Tipp'
            set(hLeg, 'Position', ... [4.9648, 8.0972, 3.5719, 0.9128]);
                [0.5125, 0.8057, 0.3922, 0.0966]);
        case 'ShrinkedWHIN'
            set(hLeg, 'Location', 'northwest');
            transparentizeCurLegends;

            % Relocate the x10^N label for y axis.
            %  annotation(hFigAreaOfInterest, 'textbox',...
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
    % Tighten the figure.
    tightfig(hFigAreaOfInterest);
    curDirToSave = fullfile(pathToSaveResults, ...
        ['Overview_UserLocGrid-', preset]);
    saveEpsFigForPaper(hFigAreaOfInterest, curDirToSave);
end

close all;
disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Simulation Parameters as an Overview

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Extracting key simulation parameters as an overview ...'])

% We will use the 1900 MHz results.
numOfPresets = length(PRESETS);
freqInMhz = 1900;

% The key parameters we will extract include (please refer to the TWC paper
% for more details):
%   Area of interest (km^2), simulation area (km^2), N_{Tower}, length of
%   the long side (km), N_{Samp}, delta_d (km), N_{User}, simulation time
%   (day).
% Simulation time needs to be estimated manually from the diary files.
[AreaOfIntInKm2, SimAreaInKm2, ...
    NTower, NTowerInArea, towerDenInNumPer1000SquKm, ...
    LengthOfLongerSideInKm, NSamp, DeltaDInKm, NUser] ...
    = deal(nan(numOfPresets, 1));
for idxPreset = 1:numOfPresets
    preset = PRESETS{idxPreset};
    disp(['        [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(freqInMhz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    clearvars simConfigs simState;
    load(fullfile(pathToReadResults, 'simConfigs.mat'));
    load(fullfile(pathToReadResults, 'simState.mat'));

    % For the extended area.
    utmXYBoundaryOfExtendedArea = extendUtmXYBoundOfInt( ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ...
        simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);

    % For the grid resolution.
    mapMinX = min(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
    mapMaxX = max(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
    mapMinY = min(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
    mapMaxY = max(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
    mapWidthInM = mapMaxX-mapMinX;
    mapHeightInM = mapMaxY-mapMinY;
    curLengthOfLongerSideInKm = max([mapWidthInM, mapHeightInM])/1000;
    gridResolutionInKm = curLengthOfLongerSideInKm ...
        ./simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE;

    % Store the parameters.
    AreaOfIntInKm2(idxPreset) ...
        = polyarea(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2)) / (1000^2);
    SimAreaInKm2(idxPreset) ...
        = polyarea(utmXYBoundaryOfExtendedArea(:,1), ...
        utmXYBoundaryOfExtendedArea(:,2)) / (1000^2);
    NTower(idxPreset) ...
        = size(simState.CellAntsEffectiveIds, 1);
    NTowerInArea(idxPreset) = sum(inpolygon( ...
        simState.CellAntsXyhEffective(:,1), ...
        simState.CellAntsXyhEffective(:,2), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2)));
    towerDenInNumPer1000SquKm(idxPreset) ...
        = NTowerInArea(idxPreset)/AreaOfIntInKm2(idxPreset)*1000;
    LengthOfLongerSideInKm(idxPreset) ...
        = curLengthOfLongerSideInKm;
    NSamp(idxPreset) ...
        = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE;
    DeltaDInKm(idxPreset) ...
        = gridResolutionInKm;
    NUser(idxPreset) ...
        = size(simState.mapGridXYPts, 1);
end

% Export the raw results to a .csv file.
curPathToSaveCsv = fullfile(pathToSaveResults, ...
    'keySimParameters_Raw.csv');
writetable(table(AreaOfIntInKm2, NTowerInArea, ...
    towerDenInNumPer1000SquKm, SimAreaInKm2, NTower, ...
    LengthOfLongerSideInKm, NSamp, DeltaDInKm, NUser), ...
    curPathToSaveCsv);

% Export the rounded parameters to another .csv file.
nthDigitToRoundTo = 1;
AreaOfIntInKm2 = round(AreaOfIntInKm2, nthDigitToRoundTo);
SimAreaInKm2 = round(SimAreaInKm2, nthDigitToRoundTo);
LengthOfLongerSideInKm = round(LengthOfLongerSideInKm, nthDigitToRoundTo);
towerDenInNumPer1000SquKm ...
    = round(towerDenInNumPer1000SquKm, nthDigitToRoundTo);

curPathToSaveCsv = fullfile(pathToSaveResults, ...
    'keySimParameters.csv');
writetable(table(AreaOfIntInKm2, NTowerInArea, ...
    towerDenInNumPer1000SquKm, SimAreaInKm2, NTower, ...
    LengthOfLongerSideInKm, NSamp, DeltaDInKm, NUser), ...
    curPathToSaveCsv);

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Recorded Simulation Time

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Estimating simulation time ...'])

% We will go through all carrier frequencies and presets.
numOfFreqs = length(CARRIER_FREQUENCIES_IN_MHZ);
numOfPresets = length(PRESETS);

% A matrix for the results. Each row corresponds to the time values for all
% presets with one carrier frequency.
%   - SimTimeInDay
%     This is ased on the tic/toc results recorded during the simulation,
%     which appear to be a way-under-estimate, possibly because some of the
%     workers may die during the simulation.
%   - SimTimeInDayBasedOnDiary
%     This is based on the diary log, which is still an under-estimate but
%     should be more accurate.
[SimTimeInDay, SimTimeInDayBasedOnDiary, SimTimeInDayBasedOnDiaryMore, ...
    numOfWorkers, minPathlossInDb, maxPathlossInDb, ...
    minBlockDistInM, maxBlockDistInM] ...
    = deal(nan(numOfFreqs, numOfPresets));
simMachineNames = cell(numOfFreqs, numOfPresets);
for idxFreq = 1:numOfFreqs
    freqInMhz = CARRIER_FREQUENCIES_IN_MHZ{idxFreq};
    for idxPreset = 1:numOfPresets
        preset = PRESETS{idxPreset};
        disp(['        [', datestr(now, datetimeFormat), ...
            '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
            '/', num2str(numOfPresets), ') at ', ...
            num2str(freqInMhz), ' MHz (carrier #', num2str(idxFreq), ...
            '/', num2str(numOfFreqs), ') ...'])

        pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'PostProcessingResults', ['Simulation_', preset, ...
            '_Carrier_', num2str(freqInMhz), ...
            'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);
        pathToSimDiary = fullfile(pathToReadResults, 'diary.txt' );

        % Note that the cache file contains simState.
        clearvars simConfigs simState;
        pathToReadCache = fullfile(pathToReadResults, ...
            'CovAnalysisCache_*.mat' );
        cacheFile = dir(pathToReadCache);

        % The full path to the cache file seems to cause trouble in the
        % loading process. This is a workaround.
        pathToCacheFileCopy = fullfile(cacheFile.folder, 'deleteme.mat');
        copyfile(fullfile(cacheFile.folder, cacheFile.name), ...
            pathToCacheFileCopy);
        load(pathToCacheFileCopy);
        delete(pathToCacheFileCopy);

        % Number of workers based on the cache file.
        numOfWorkers(idxFreq, idxPreset) ...
            = length(locIndicesForAllWorkersForAllCellsEff{1});

        % Host name of the machine on which the simulation was carried out.
        simMachineNames{idxFreq, idxPreset} ...
            = parseCacheFilename(cacheFile.name);

        % The final version of simState.
        clearvars simConfigs simState;
        load(fullfile(pathToReadResults, 'simState.mat'));

        % Time used for all heights and pixels, according to the simState
        % and cache files.
        timeUsedForAllMapSets = [simState.TimeUsedInSForEachPixel{:}];
        SimTimeInS = [timeUsedForAllMapSets{:}]';
        SimTimeInDay(idxFreq, idxPreset) ...
            = sum(SimTimeInS)/60/60/24/numOfWorkers(idxFreq, idxPreset);

        % Time used for all presets at each carrier, according to the diary
        % log.
        numOfEffeTowers = size(simState.CellAntsEffectiveIds, 1);
        [SimTimeInDayBasedOnDiary(idxFreq, idxPreset), ~, ...
            SimTimeInDayBasedOnDiaryMore(idxFreq, idxPreset)] ...
            = extractSimTimeInDayFromDiary( ...
            pathToSimDiary, numOfEffeTowers);

        % Other parameters of interest.
        minPathlossInDb(idxFreq, idxPreset) ...
            = min(vertcat(simState.coverageMaps{:}));
        maxPathlossInDb(idxFreq, idxPreset) ...
            = max(vertcat(simState.coverageMaps{:}));
        minBlockDistInM(idxFreq, idxPreset) ...
            = min(vertcat(simState.blockageDistMaps{:}));
        maxBlockDistInM(idxFreq, idxPreset) ...
            = max(vertcat(simState.blockageDistMaps{:}));
    end
end

% Sim time in days and in hours based on simState.
exportSimTimeToCsv(SimTimeInDay, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDays.csv'), 1);
exportSimTimeToCsv(SimTimeInDay, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDays_Raw.csv'));

exportSimTimeToCsv(SimTimeInDay*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInH.csv'), 1);
exportSimTimeToCsv(SimTimeInDay*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInH_Raw.csv'));

% Sim time in days and in hours based on diary logs.
exportSimTimeToCsv(SimTimeInDayBasedOnDiary, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDaysBasedOnDiary.csv'), 1);
exportSimTimeToCsv(SimTimeInDayBasedOnDiary, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDaysBasedOnDiary_Raw.csv'));

exportSimTimeToCsv(SimTimeInDayBasedOnDiary*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInHBasedOnDiary.csv'), 1);
exportSimTimeToCsv(SimTimeInDayBasedOnDiary*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInHBasedOnDiary_Raw.csv'));

% Longer sim time in days and in hours based on diary logs.
exportSimTimeToCsv(SimTimeInDayBasedOnDiaryMore, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDaysBasedOnDiaryMore.csv'), 1);
exportSimTimeToCsv(SimTimeInDayBasedOnDiaryMore, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInDaysBasedOnDiaryMore_Raw.csv'));

exportSimTimeToCsv(SimTimeInDayBasedOnDiaryMore*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInHBasedOnDiaryMore.csv'), 1);
exportSimTimeToCsv(SimTimeInDayBasedOnDiaryMore*24, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'SimTimeInHBasedOnDiaryMore_Raw.csv'));

% Parameter ranges.
exportSimTimeToCsv(minPathlossInDb, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'minPathlossInDb.csv'), 1);
exportSimTimeToCsv(maxPathlossInDb, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'maxPathlossInDb.csv'), 1);
exportSimTimeToCsv(minBlockDistInM, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'minBlockDistInM.csv'), 1);
exportSimTimeToCsv(maxBlockDistInM, ...
    PRESETS, CARRIER_FREQUENCIES_IN_MHZ, ...
    fullfile(pathToSaveResults, 'maxBlockDistInM.csv'), 1);

% Hostnames and worker numbers.
writecell(simMachineNames, ...
    fullfile(pathToSaveResults, 'SimTimeHostname.csv'));
writematrix(numOfWorkers, ...
    fullfile(pathToSaveResults, 'SimTimeWorkerNums.csv'));

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Example Blockage Figs (Status and Distance for Tipp + WHIN + IN)

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating example blockage status and distance maps ...'])

curFlagGenFigsSilently = false;
curFlagZoomIn = true;

% We will use the 1900 MHz results.
curPresets = PRESETS;
numOfPresets = length(curPresets);
freqInMhz = 1900;
hightsToInspect = [1.5, 3, 5, 7.5, 10, 50, 100];
indicesH = [1:5, 9, 14];
% For legends on maps for Tipp and IN.
staLegendsShown = [false, false];
% For distance scales on maps for Tipp, WHIN, and IN.
distLegendsShown = [false, false, false];

% This should be big enough for all the simulations.
maxBlockDistInM = 1000;
for idxPreset = 1:numOfPresets
    preset = curPresets{idxPreset};
    disp(['        [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(freqInMhz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    clearvars simConfigs simState;
    load(fullfile(pathToReadResults, 'simConfigs.mat'));
    load(fullfile(pathToReadResults, 'simState.mat'));

    [effeCellAntLats, effeCellAntLons] ...
        = simConfigs.utm2deg_speZone( ...
        simState.CellAntsXyhEffective(:, 1), ...
        simState.CellAntsXyhEffective(:, 2));
    effeCellAntLonLats = [effeCellAntLons, effeCellAntLats];

    mapGridLonLats = simState.mapGridLatLonPts(:, [2,1]);

    for idxH = indicesH
        rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
        assert(ismember(rxAntH, hightsToInspect), ...
            'idxH is wrong according to hightsToInspect!');

        %---------------
        % For blockage status maps.
        %---------------
        flagRiseTxToTop = true;
        surfFaceAlpha = 0.5;
        % Smaller maps for publication.
        %   - [500, 500].*0.6 was used for the ICC 2020 paper.
        curBlockStaFigSize = [500, 500].*0.6;
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            flagRiseTxToTop = false;
            surfFaceAlpha = 0.6;
            curBlockStaFigSize = [500, 500].*0.75;
        end
        [ hCurBlStaMap ] ...
            = plotBlockageMap( ...
            [mapGridLonLats, simState.blockageMaps{idxH}], ...
            effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
            curFlagZoomIn, curBlockStaFigSize, [1 1 0; 1 0 0], ...
            flagRiseTxToTop, surfFaceAlpha);
        hLeg = findobj(hCurBlStaMap, 'Type', 'Legend');
        % Tighten the figure.
        xlabel(''); ylabel('');
        tightfig(hCurBlStaMap);

        % Hide the legend except for the first Tipp figure and the first IN
        % figure.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
            if (~staLegendsShown(1))
                set(hLeg, 'Position', [3.1070, 5.3797, 2.8840, 1.2435]);
                staLegendsShown(1) = true;
            else
                legend off;
            end
        elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            if (~staLegendsShown(2))
                set(hLeg, 'Position', [1.6360, 0.1522, 2.8046, 1.2435]);
                staLegendsShown(2) = true;
            else
                legend off;
            end
        else
            legend off;
        end

        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageStatusMap_', simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        saveas(hCurBlStaMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurBlStaMap, [pathToSaveFig, '.png']);
        saveas(hCurBlStaMap, [pathToSaveFig, '.fig']);

        %---------------
        % For blockage distance maps.
        %---------------
        flagRiseTxToTop = true;
        txMarkerSize = 6;
        curBlockDistFigSize = [500, 500].*0.6;
        flagScatterForTx = false;
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            % flagRiseTxToTop = true; txMarkerSize = 3;
            txMarkerAlpha = 0.75;
            curBlockDistFigSize = [500, 500].*0.75;
            flagScatterForTx = true;
        end
        [ hCurDistMap, hTxTs, hCb ] = plotPathLossMap( ...
            [mapGridLonLats, simState.blockageDistMaps{idxH}], ...
            effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
            curFlagZoomIn, 'griddatasurf', curBlockDistFigSize, ...
            customHot, flagRiseTxToTop, 'k', ...
            txMarkerSize, flagScatterForTx);

        % Resize the figure for WHIN and IN plots.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedwhin')
            set(hCurDistMap, 'Position', [0, 0, 245, 300]);
            axis([-87.17148148, -86.13560255, 39.83917848, 41.19774935]);
        elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            set(hCurDistMap, 'Position', [0, 0, 275, 375]);
            axis([-87.49255328, -85.09043091, 38.11187646, 41.55153020]);

            % Also transparentize the Tx tower markers.
            alpha(hTxTs, txMarkerAlpha);
        end

        % Tighten the figure.
        xlabel(''); ylabel('');
        pause(1);
        tightfig(hCurDistMap);

        % Hide the legend except for the first figure of Tipp and the first
        % figure of IN.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
            if ~distLegendsShown(1)
                hLeg = findobj(hCurDistMap, 'Type', 'Legend');
                set(hLeg, 'Position', [2.8301, 6.1344, 2.8840, 0.4762], ...
                    'AutoUpdate', 'off');

                distLegendsShown(1) = true;
            else
                legend off;
            end
        elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            if ~distLegendsShown(3)
                hLeg = findobj(hCurDistMap, 'Type', 'Legend');
                set(hLeg, 'Position', [0.1468, 7.7788, 2.8046, 0.4762], ...
                    'AutoUpdate', 'off');

                distLegendsShown(3) = true;
            else
                legend off;
            end
        else
            legend off;
        end

        % Always show the distance scale.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
            h = makescale(3.5, 'se', 'units', 'si');
            % The scale will be blocked by the plot if not adjusted along
            % the z axis.
            h(1).ZData = ones(4,1).*largeZValue;
            h(2).ZData = ones(2,1).*largeZValue;
            set(h(3), 'Position', tippScalePos);
            h(3).FontSize = tippScaleFontS;
        end

        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedwhin')
            % Show distance scale for the first figure of WHIN. if
            % (~distLegendsShown(2))&&(strcmpi( ...
            %         simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedwhin'))
            h = makescale(2.9, 'se', 'units', 'si');
            % The scale will be blocked by the plot if not adjusted along
            % the z axis.
            h(1).ZData = ones(4,1).*largeZValue;
            h(2).ZData = ones(2,1).*largeZValue;
            set(h(3), 'Position', whinScalePos);
            h(3).FontSize = whinScaleFontS;
        end

        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            % Show distance scale for the first figure of IN. if
            % (~distLegendsShown(2))&&(strcmpi( ...
            %         simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedwhin'))
            h = makescale(3.45, 'se', 'units', 'si');
            % The scale will be blocked by the plot if not adjusted along
            % the z axis.
            h(1).ZData = ones(4,1).*largeZValue;
            h(2).ZData = ones(2,1).*largeZValue;
            set(h(3), 'Position', inScalePos);
            h(3).FontSize = inScaleFontS;
        end

        % Mannual figure adjustments.
        %  if strcmpi( ...
        %         simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedwhin')
        %     pause(1);
        %      set(hCurDistMap, 'Position', [0, 0, 240, 275]);
        %     set(gca, 'Position', [0.0601, 0.0208, 0.6478, 0.8964]);
        %      set(hCb, 'Position', [4.8502, 0.1500, 0.3360, 6.4823]);
        %  end

        % Move the colorbar title a little to the left so that the unit can
        % be seen.
        title(hCb, 'Blockage Distance (m)            ');

        % Hide the color bar except for the first figure of WHIN.
        %  if (~distLegendsShown(2)) ...
        %         && (strcmpi( ... simConfigs.CURRENT_SIMULATION_TAG,
        %         'shrinkedwhin'))
        %     set(hCb, 'Visible', 'off');
        %      distLegendsShown(2) = true;
        %  end

        % With different color bar for each plot.
        pause(1);
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageDistMap_', simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
        set(hCurDistMap, 'PaperPositionMode', 'auto');
        saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurDistMap, [pathToSaveFig, '.png']);
        saveas(hCurDistMap, [pathToSaveFig, '.fig']);

        % The same color bar is shared in different figures.
        caxis([0, maxBlockDistInM]);

        pause(1);
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageDistMap_SameCB_', simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
        saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurDistMap, [pathToSaveFig, '.png']);
        saveas(hCurDistMap, [pathToSaveFig, '.fig']);

        % The same color bar is shared in different figures but hidden.
        colorbar off;

        pause(1);
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageDistMap_SameCBHidden_', ...
            simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
        saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurDistMap, [pathToSaveFig, '.png']);
        saveas(hCurDistMap, [pathToSaveFig, '.fig']);

        % Make a copy of the figure;
        a1 = gca;
        hCurDistMapCopy = figure('Position', get(hCurDistMap, 'Position'));
        a2 = copyobj(a1,hCurDistMapCopy);
        colormap(customHot);

        % The same color bar is shared in different figures but hidden;
        % figure is tighten, too.
        figure(hCurDistMap);
        title('');
        % There is a bug with tightfig; an unintended empty area shows up
        % at the top. This is a workaround for that. Note that the figure
        % resizing will take some time. The pause command will make sure
        % things take effect as expected.
        pause(1);
        tightfig(hCurDistMap);
        pause(1);
        set(gca, 'Unit', 'pixel'); set(gcf, 'Unit', 'pixel');
        curAxesPos = get(gca, 'Position');
        set(gcf, 'Position', [0,0, curAxesPos(3:4)+10]);

        pause(1);
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageDistMap_SameCBHiddenTighten_', ...
            simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        % We may need to re-adjust the scale.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
            set(h(3), 'Position', tippScalePos);
            h(3).FontSize = tippScaleFontS;
        end
        set(hCurDistMap, 'PaperPositionMode', 'auto');
        pause(1);
        export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
        saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurDistMap, [pathToSaveFig, '.png']);
        saveas(hCurDistMap, [pathToSaveFig, '.fig']);

        % The same color bar is shared in different figures but hidden;
        % figure is tighten only vertically.
        curFigPos = get(hCurDistMap, 'Position');
        oldFigPos = get(hCurDistMapCopy, 'Position');

        figure(hCurDistMapCopy);
        set(gca, 'Unit', 'pixel');
        oldAxesPos = get(gca, 'Position');
        set(hCurDistMapCopy, 'Position', [oldFigPos(1:3), curFigPos(4)])
        set(gca, 'Position', oldAxesPos);

        pause(1);
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageDistMap_SameCBHiddenTightenVert_', ...
            simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        % We may need to re-adjust the scale.
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
            set(h(3), 'Position', tippScalePos);
            h(3).FontSize = tippScaleFontS;
        end
        % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
        saveas(hCurDistMapCopy, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurDistMapCopy, [pathToSaveFig, '.png']);
        saveas(hCurDistMapCopy, [pathToSaveFig, '.fig']);
    end
end

close all;
disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Example Path Loss Figs (Different F_C for Tipp + WHIN + IN)

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating example path loss maps ...'])

curFlagGenFigsSilently = false;
curFlagZoomIn = true;
defaultCmdToPlotPLMaps = 'surf';

% We will use the h_R = 1.5m results.
curPresets = PRESETS;
numOfPresets = length(curPresets);
freqsToInspectInMhz = [1900, 4700, 7000, 28000];
hightsToInspect = [1.5, 3, 5, 7.5, 10, 50, 100];
indicesH = [1:5, 9, 14];

% For better coloring effects.
customCAxis = [120, 150];
txMarkerAlpha = 0.75;

% For legends on maps for Tipp and IN.
txLegendsShown = [false, false];

% Only inspect frequencies specified by CARRIER_FREQUENCIES_IN_MHZ.
freqsToInspectInMhz = freqsToInspectInMhz( ...
    ismember(freqsToInspectInMhz, [CARRIER_FREQUENCIES_IN_MHZ{:}]));
numOfFs = length(freqsToInspectInMhz);
for idxF = 1:numOfFs
    fcInMHz = freqsToInspectInMhz(idxF);
    if ismember(fcInMHz, [CARRIER_FREQUENCIES_IN_MHZ{:}])
        assert(ismember(fcInMHz, freqsToInspectInMhz), ...
            'idxF is wrong according to freqsToInspectInMhz!');

        disp(['        [', datestr(now, datetimeFormat), ...
            '] Frequency #', num2str(idxF), ...
            '/', num2str(numOfFs), ': ', num2str(fcInMHz), ' MHz ...'])

        for idxPreset = 1:numOfPresets
            preset = curPresets{idxPreset};
            disp(['            [', datestr(now, datetimeFormat), ...
                '] Processing ', preset, ' (preset #', ...
                num2str(idxPreset), '/', num2str(numOfPresets), ') ...'])

            pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
                'PostProcessingResults', ['Simulation_', preset, ...
                '_Carrier_', num2str(fcInMHz), ...
                'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

            clearvars simConfigs simState;
            try
                load(fullfile(pathToReadResults, 'simConfigs.mat'));
            catch
                load(fullfile(pathToReadResults, 'simConfigs_Raw.mat'));
            end
            load(fullfile(pathToReadResults, 'simState.mat'));

            for countH = 1:length(indicesH)
                idxH = indicesH(countH);
                hightToInspect = hightsToInspect(countH);

                [effeCellAntLats, effeCellAntLons] ...
                    = simConfigs.utm2deg_speZone( ...
                    simState.CellAntsXyhEffective(:, 1), ...
                    simState.CellAntsXyhEffective(:, 2));
                effeCellAntLonLats = [effeCellAntLons, effeCellAntLats];

                mapGridLonLats = simState.mapGridLatLonPts(:, [2,1]);

                rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
                assert(rxAntH==hightToInspect, ...
                    'idxH is wrong according to hightToInspect!');

                %---------------
                % For path loss maps.
                %---------------
                % Smaller maps for publication.
                %   - [500, 500].*0.6 was used for the ICC 2020 paper.
                curPLFigSize = [500, 500].*0.6;
                flagScatterForTx = false;
                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, ...
                        'shrinkedin')
                    curPLFigSize = [500, 500].*0.75;
                    flagScatterForTx = true;
                end
                [ hCurPLMap, hTxTs, hCb] = plotPathLossMap( ...
                    [mapGridLonLats, simState.coverageMaps{idxH}], ...
                    effeCellAntLonLats, simConfigs, ...
                    ~curFlagGenFigsSilently, curFlagZoomIn, ...
                    defaultCmdToPlotPLMaps, curPLFigSize, ...
                    customHotLong, true, 'k', 6, flagScatterForTx);

                caxis(customCAxis);
                xlabel(''); ylabel('');
                tightfig(hCurPLMap);
                pause(3);

                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, ...
                        'shrinkedwhin')
                    set(hCurPLMap, 'Position', [0, 0, 215, 275]);
                    pause(1);
                    axis([-87.18098142, -86.12610260, ...
                        39.83917848, 41.19774935]);
                    set(hCb, 'Position', [4.2946, 0.1500, 0.3360, 6.51]);
                elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, ...
                        'shrinkedin')
                    alpha(hTxTs, txMarkerAlpha);
                    % Resize the figure for IN plots.
                    set(hCurPLMap, 'Position', [0, 0, 225, 335]);
                    pause(1);
                    axis([-87.47910844, -85.12814685, ...
                        38.11617027, 41.53260018]);
                    set(hCb, 'Position', [4.6011, 0.15, 0.3915, 8.0962]);
                end

                % Hide the legend except for the first Tipp figure and the
                % first IN figure.
                if idxF==1
                    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
                        if ~txLegendsShown(1)
                            hLeg = findobj(hCurPLMap, 'Type', 'Legend');
                            set(hLeg, 'Position', ...
                                [2.8254, 6.11, 2.8840, 0.4762], ...
                                'AutoUpdate', 'off');

                            txLegendsShown(1) = true;
                        else
                            legend off;
                        end
                    elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, ...
                            'shrinkedin')
                        if ~txLegendsShown(2)
                            hLeg = findobj(hCurPLMap, 'Type', 'Legend');
                            set(hLeg, 'Position', ...
                                [0.1486, 7.7788, 2.8045, 0.4762], ...
                                'AutoUpdate', 'off');

                            txLegendsShown(2) = true;
                        else
                            legend off;
                        end
                    else
                        legend off;
                    end
                else
                    legend off;
                end

                % Always show the distance scale.
                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
                    h = makescale(3.5, 'se', 'units', 'si');
                    % The scale will be blocked by the plot if not adjusted
                    % along the z axis.
                    h(1).ZData = ones(4,1).*largeZValue;
                    h(2).ZData = ones(2,1).*largeZValue;
                    set(h(3), 'Position', tippScalePos);
                    h(3).FontSize = tippScaleFontS;
                end

                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, ...
                        'shrinkedwhin')
                    % Show distance scale for the first figure of WHIN. if
                    % (~distLegendsShown(2))&&(strcmpi( ...
                    %         simConfigs.CURRENT_SIMULATION_TAG,
                    %         'shrinkedwhin'))
                    h = makescale(3, 'se', 'units', 'si');
                    % The scale will be blocked by the plot if not adjusted
                    % along the z axis.
                    h(1).ZData = ones(4,1).*largeZValue;
                    h(2).ZData = ones(2,1).*(largeZValue+1);
                    set(h(3), 'Position', whinScalePos);
                    h(3).FontSize = whinScaleFontS;
                end

                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
                    % Show distance scale for the first figure of WHIN. if
                    % (~distLegendsShown(2))&&(strcmpi( ...
                    %         simConfigs.CURRENT_SIMULATION_TAG,
                    %         'shrinkedwhin'))
                    h = makescale(3.35, 'se', 'units', 'si');
                    % The scale will be blocked by the plot if not adjusted
                    % along the z axis.
                    h(1).ZData = ones(4,1).*largeZValue;
                    h(2).ZData = ones(2,1).*(largeZValue+1);
                    set(h(3), 'Position', inScalePos);
                    h(3).FontSize = whinScaleFontS;
                end

                pathToSaveFig = fullfile(pathToSaveResults, ...
                    ['PathLossMap_', simConfigs.CURRENT_SIMULATION_TAG, ...
                    '_Fc_', num2str(fcInMHz), 'MHz_RxHeight_', ...
                    strrep(num2str(rxAntH), '.', '_')]);
                set(hCurPLMap, 'PaperPositionMode', 'auto');
                pause(1);
                % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
                saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
                saveas(hCurPLMap, [pathToSaveFig, '.png']);
                saveas(hCurPLMap, [pathToSaveFig, '.fig']);

                % The version with colorbar hidden.
                colorbar off;

                % Make a copy of the figure;
                a1 = gca;
                hCurPLMapCopy = figure('Position', ...
                    get(hCurPLMap, 'Position'));
                a2 = copyobj(a1,hCurPLMapCopy);
                colormap(customHot);

                % The same color bar is shared in different figures but
                % hidden; figure is tighten, too.
                figure(hCurPLMap);
                title('');
                % There is a bug with tightfig; an unintended empty area
                % shows up at the top. This is a workaround for that. Note
                % that the figure resizing will take some time. The pause
                % command will make sure things take effect as expected.
                pause(1);
                tightfig(hCurPLMap);
                pause(1);
                set(gca, 'Unit', 'pixel'); set(gcf, 'Unit', 'pixel');
                curAxesPos = get(gca, 'Position');
                set(gcf, 'Position', [0,0, curAxesPos(3:4)+10]);
                % We may need to re-adjust the scale.
                if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'tipp')
                    set(h(3), 'Position', tippScalePos);
                    h(3).FontSize = tippScaleFontS;
                end
                pause(1);

                pathToSaveFig = fullfile(pathToSaveResults, ...
                    ['PathLossMap_CBHidden_', ...
                    simConfigs.CURRENT_SIMULATION_TAG, ...
                    '_Fc_', num2str(fcInMHz), 'MHz_RxHeight_', ...
                    strrep(num2str(rxAntH), '.', '_')]);

                % To avoid resizing figure before saving.
                set(hCurPLMap, 'PaperPositionMode', 'auto');
                % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
                saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
                saveas(hCurPLMap, [pathToSaveFig, '.png']);
                saveas(hCurPLMap, [pathToSaveFig, '.fig']);

                % The same color bar is shared in different figures but
                % hidden; figure is tighten only vertically.
                curFigPos = get(hCurPLMap, 'Position');
                oldFigPos = get(hCurPLMapCopy, 'Position');

                figure(hCurPLMapCopy);
                set(gca, 'Unit', 'pixel');
                oldAxesPos = get(gca, 'Position');
                set(hCurPLMapCopy, 'Position', ...
                    [oldFigPos(1:3), curFigPos(4)])
                set(gca, 'Position', oldAxesPos);
                pause(1);

                pathToSaveFig = fullfile(pathToSaveResults, ...
                    ['PathLossMap_CBHidden_TightenVert_', ...
                    simConfigs.CURRENT_SIMULATION_TAG, ...
                    '_Fc_', num2str(fcInMHz), 'MHz_RxHeight_', ...
                    strrep(num2str(rxAntH), '.', '_')]);
                % export_fig([pathToSaveFig, '_export_fig.eps'], '-eps');
                saveas(hCurPLMapCopy, [pathToSaveFig, '.eps'], 'epsc');
                saveas(hCurPLMapCopy, [pathToSaveFig, '.png']);
                saveas(hCurPLMapCopy, [pathToSaveFig, '.fig']);

                close(hCurPLMap);
                close(hCurPLMapCopy);
            end
        end
    end
end

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Coverage Ratio Plots (for Tipp, WHIN, and IN at 1900MHz)

% Coverage Ratio Plots - Blockage
%   1. LoS coverage ratio over inspected heights.
%    2. Coverage ratio over max allowed block dist.
%   3. Coverage ratio gain over max allowed block dist.

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating coverage ratio over relay height plot: blockage ...'])

curPresets = PRESETS;
numOfPresets = length(curPresets);
freqsToInspectInMhz = 1900;
numOfFs = length(freqsToInspectInMhz);
curFlagGenFigsSilently = false;
lineStyles = {'-o', '--^', ':s'};
flagResizeFigForPublication = true;

% For coverage ratio plots.
%
% Before adding new heights below 10 m:
%      heightsInM = [1.5, 10, 20, 30, 40, 50, 80, 125];
%     heightsIndices = [1, 2, 3, 4, 5, 6, 9, 14];
%
%      heightsInM = [1.5, 10, 30, 50, 70, 100, 125];
%     heightsIndices = [1, 2, 4, 6, 8, 11, 14];
%
%      heightsInM = [1.5, 10:10:120, 125];
%     heightsIndices = 1:14;
%
%      heightsInM = [1.5, 10, 20:20:120, 125];
%     heightsIndices = [1:2 3:2:13, 14];
%
% After:
%       simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
%           = [1.5; 3; 5; 7.5; (10:10:120)'; 125];
%
heightsInM = [1.5, 3, 5, 7.5, 10, 20:20:120, 125];
heightsIndices = [1:5, 6:2:16, 17];

visibleBlockDistRange = [1, 1000];
curManualXLim = true;
curShowFSPL = true;

hLosCovRatOverHFig= figure('Position', [0, 0, 500, 180]);
hold on; set(gca, 'fontWeight', 'bold');
xlabel('Relay Height (m)'); ylabel('LoS Coverage Ratio');
grid on; grid minor;

idxF = 1;
fcInMHz = freqsToInspectInMhz(idxF);
assert(ismember(fcInMHz, freqsToInspectInMhz), ...
    'idxF is wrong according to freqsToInspectInMhz!');

disp(['        [', datestr(now, datetimeFormat), ...
    '] Frequency #', num2str(idxF), ...
    '/', num2str(numOfFs), ': ', num2str(fcInMHz), ' MHz ...'])

% For debugging, cache the intermediate results.
[hLegs, covRatios, ...
    coverageMapsCovRatioMetas, blockageDistMapsCovRatioMetas] ...
    = deal(cell(numOfPresets, 1));
for idxPreset = 1:numOfPresets
    preset = curPresets{idxPreset};
    disp(['            [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(fcInMHz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    clearvars simConfigs simState;
    try
        load(fullfile(pathToReadResults, 'simConfigs.mat'));
    catch
        load(fullfile(pathToReadResults, 'simConfigs_Raw.mat'));
    end
    load(fullfile(pathToReadResults, 'simState.mat'));

    assert( all( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(heightsIndices) ...
        == heightsInM'), 'heightsIndices does not match with heightsInM!');

    %---------------
    % LoS coverage ratio over inspected heights.
    %---------------
    [curCovRatioVsHFig, covRatios{idxPreset}] ...
        = plotCovRatioVsInspectedHeight(simState.blockageMaps, ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M, ...
        ~curFlagGenFigsSilently);
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['CoverageRatioVsDroneHeight_Blockage_', ...
        simConfigs.CURRENT_SIMULATION_TAG, ...
        '_Fc_', num2str(fcInMHz), 'MHz']);
    figure(curCovRatioVsHFig);
    saveEpsFigForPaper(curCovRatioVsHFig,  pathToSaveFig);
    close(curCovRatioVsHFig);

    figure(hLosCovRatOverHFig);
    hLegs{idxPreset} = plot( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M, ...
        covRatios{idxPreset}, lineStyles{idxPreset}, ...
        'MarkerSize', 6, 'LineWidth', 1);

    %---------------
    % Coverage ratio over max allowed block dist.
    %---------------
    mapType = 'BlockageDist';

    clearvars tempSimState tempSimConfigs;
    tempSimState.blockageDistMaps ...
        = simState.blockageDistMaps(heightsIndices);
    tempSimState.mapGridXYPts = simState.mapGridXYPts;
    tempSimState.CellAntsXyhEffective = simState.CellAntsXyhEffective;
    tempSimConfigs.CARRIER_WAVELENGTH_IN_M ...
        = simConfigs.CARRIER_WAVELENGTH_IN_M;
    tempSimConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
        = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(heightsIndices);

    [curEmpCdfFig, blockageDistMapsCovRatioMetas] ...
        = plotEmpiricalCdfForCoverage(tempSimState, tempSimConfigs, ...
        mapType, ~curFlagGenFigsSilently);
    xlim(visibleBlockDistRange);
    if idxPreset==1
        hLeg = findobj(curEmpCdfFig, 'Type', 'Legend');
        set(hLeg, 'Position', [0.6367, 0.1522, 0.2700, 0.6717]);
        transparentizeCurLegends;
    else
        legend off;
    end
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType, ...
        simConfigs.CURRENT_SIMULATION_TAG, ...
        '_Fc_', num2str(fcInMHz), 'MHz']);
    saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
    close(curEmpCdfFig);

    %---------------
    % Coverage ratio gain over max allowed block dist.
    %---------------
    [ curDistCovRatioGainFig ] ...
        = plotCoverageRatioGain(tempSimState, tempSimConfigs, ...
        mapType, ~curFlagGenFigsSilently);
    xlim(visibleBlockDistRange);
    if idxPreset==1
        hLeg = findobj(curDistCovRatioGainFig, 'Type', 'Legend');
        set(hLeg, 'Position', [0.6367, 0.3111, 0.2700, 0.6167]);
        transparentizeCurLegends;
    else
        legend off;
    end
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType, 'CovRatGain_', ...
        simConfigs.CURRENT_SIMULATION_TAG, ...
        '_Fc_', num2str(fcInMHz), 'MHz']);
    saveEpsFigForPaper(curDistCovRatioGainFig, pathToSaveFig);
    close(curDistCovRatioGainFig);
end

figure(hLosCovRatOverHFig);
xlim([0, 50]); ylim([0, 1]);
legend([hLegs{:}], 'Tippecanoe', 'WHIN', 'IN', 'Location', 'se');
pathToSaveFig = fullfile(pathToSaveResults, ...
    'CoverageRatioVsDroneHeight_Blockage');
saveEpsFigForPaper(hLosCovRatOverHFig,  pathToSaveFig);
close(hLosCovRatOverHFig);

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

% Coverage Ratio Plots - Path Loss
%   1.Coverage ratio over max allowed path loss.
%    2. Coverage ratio gain over max allowed path loss.

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating coverage ratio over relay height plot: path loss ...'])

% eHata results for heights over 10 m are less trustworthy.
heightsInM = [1.5, 3, 5, 7.5, 10, 50, 100, 125];
heightsIndices = [1:5, 9, 14, 17];

visiblePathLossRange = [100, 140];

% FSPL calculated for clear LoS links at 50 m as the best performance
% limit.
refBlockageMapIdx = 9;
refBlockageMapRxHInM = 50;

idxF = 1;
fcInMHz = freqsToInspectInMhz(idxF);
assert(ismember(fcInMHz, freqsToInspectInMhz), ...
    'idxF is wrong according to freqsToInspectInMhz!');

disp(['        [', datestr(now, datetimeFormat), ...
    '] Frequency #', num2str(idxF), ...
    '/', num2str(numOfFs), ': ', num2str(fcInMHz), ' MHz ...'])

for idxPreset = 1:numOfPresets
    preset = curPresets{idxPreset};
    disp(['            [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(fcInMHz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    clearvars simConfigs simState;
    try
        load(fullfile(pathToReadResults, 'simConfigs.mat'));
    catch
        load(fullfile(pathToReadResults, 'simConfigs_Raw.mat'));
    end
    load(fullfile(pathToReadResults, 'simState.mat'));

    assert( all( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(heightsIndices) ...
        == heightsInM'), 'heightsIndices does not match with heightsInM!');

    %---------------
    % Coverage ratio over max allowed path loss.
    %---------------
    % For coverage maps.
    mapType = 'Coverage';

    clearvars tempSimState tempSimConfigs;
    tempSimState.coverageMaps = simState.coverageMaps(heightsIndices);
    tempSimState.mapGridXYPts = simState.mapGridXYPts;
    tempSimState.CellAntsXyhEffective = simState.CellAntsXyhEffective;
    tempSimConfigs.CARRIER_WAVELENGTH_IN_M ...
        = simConfigs.CARRIER_WAVELENGTH_IN_M;
    tempSimConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M ...
        = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(heightsIndices);
    tempSimConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB = visiblePathLossRange;

    % FSPL calculated for clear LoS links at specified height as the best
    % performance limit.
    assert(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M( ...
        refBlockageMapIdx) ...
        == refBlockageMapRxHInM, ...
        ['refBlockageMapIdx does not match with ', ...
        'refBlockageMapRxHInM!'])
    curRefFsplVs = simState.blockageMaps{refBlockageMapIdx};

    [curEmpCdfFig, coverageMapsCovRatioMetas{idxPreset}] ...
        = plotEmpiricalCdfForCoverage(tempSimState, tempSimConfigs, ...
        mapType, ~curFlagGenFigsSilently, curManualXLim, ...
        curShowFSPL, curRefFsplVs);
    if idxPreset==1
        hLeg = findobj(curEmpCdfFig, 'Type', 'Legend');
        set(hLeg, 'Position', [0.1293, 0.3636, 0.2880, 0.5617]);
        transparentizeCurLegends;
    else
        legend off;
    end
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType, ...
        simConfigs.CURRENT_SIMULATION_TAG, ...
        '_Fc_', num2str(fcInMHz), 'MHz']);
    saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
    close(curEmpCdfFig);

    %---------------
    % Coverage ratio gain over max allowed path loss.
    %---------------
    [ curCovRatioGainFig ] ...
        = plotCoverageRatioGain(tempSimState, tempSimConfigs, ...
        mapType, ~curFlagGenFigsSilently, curManualXLim, ...
        curShowFSPL, curRefFsplVs);
    if idxPreset==1
        hLeg = findobj(curCovRatioGainFig, 'Type', 'Legend');
        set(hLeg, 'Position', [0.1293, 0.4211, 0.2880, 0.5067]);
        transparentizeCurLegends;
    else
        legend off;
    end
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType, 'CovRatGain_', ...
        simConfigs.CURRENT_SIMULATION_TAG, ...
        '_Fc_', num2str(fcInMHz), 'MHz']);
    saveEpsFigForPaper(curCovRatioGainFig, pathToSaveFig);
    close(curCovRatioGainFig);
end

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Blockage and Path Loss Maps for Indiana
% We will focus on the 1900 MHz case with selected relay heights.

%% Tower Density for Each County

%% Also Plots for Future Paper
% We need two selected areas. One for urban environment and the other for
% rual environment.

%% Cleanup

close all;
diary off;

% EOF