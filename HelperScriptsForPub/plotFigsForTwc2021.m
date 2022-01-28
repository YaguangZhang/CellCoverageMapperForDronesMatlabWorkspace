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
    disp(['        [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

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
    % Tighten the figure.
    %   tightfig(hFigCellOverview);

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
        simState.mapGridXYPts(:,2), '.', 'MarkerSize', 3, ...
        'Color', darkBlue);
    adjustFigSizeByContent(hFigAreaOfInterest, [], 'height', 0.9);
    axis equal; view(2); %grid on; grid minor;
    hLeg = legend([hAreaOfInterest, hGridPts], ...
        'Area of interest', 'User location');
    xlabel('UTM x (m)'); ylabel('UTM y (m)'); box on;
    % Adjust legend and the exponent label for y axis.
    switch simConfigs.CURRENT_SIMULATION_TAG
        case 'Tipp'
            set(hLeg, 'Position', [0.4789, 0.8062, 0.4258, 0.0965]);
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
[AreaOfIntInKm2, SimAreaInKm2, NTower, ...
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
writetable(table(AreaOfIntInKm2, SimAreaInKm2, NTower, ...
    LengthOfLongerSideInKm, NSamp, DeltaDInKm, NUser), ...
    curPathToSaveCsv)

% Export the rounded parameters to another .csv file.
nthDigitToRoundTo = 1;
AreaOfIntInKm2 = round(AreaOfIntInKm2, nthDigitToRoundTo);
SimAreaInKm2 = round(SimAreaInKm2, nthDigitToRoundTo);
LengthOfLongerSideInKm = round(LengthOfLongerSideInKm, nthDigitToRoundTo);

curPathToSaveCsv = fullfile(pathToSaveResults, ...
    'keySimParameters.csv');
writetable(table(AreaOfIntInKm2, SimAreaInKm2, NTower, ...
    LengthOfLongerSideInKm, NSamp, DeltaDInKm, NUser), ...
    curPathToSaveCsv)

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
[SimTimeInDay, SimTimeInDayBasedOnDiary, numOfWorkers] ...
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

        load(fullfile(pathToReadResults, 'simState.mat'));

        pathToReadCache = fullfile(pathToReadResults, ...
            'CovAnalysisCache_*.mat' );
        cacheFile = dir(pathToReadCache);

        pathToSimDiary = fullfile(pathToReadResults, 'diary.txt' );

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
        % Time used for all heights and pixels, according to the simState
        % and cache files.
        timeUsedForAllMapSets = [simState.TimeUsedInSForEachPixel{:}];
        SimTimeInS = [timeUsedForAllMapSets{:}]';
        SimTimeInDay(idxFreq, idxPreset) ...
            = sum(SimTimeInS)/60/60/24/numOfWorkers(idxFreq, idxPreset);

        % Time used for all presets at each carrier, according to the diary
        % log.
        numOfEffeTowers = size(simState.CellAntsEffectiveIds, 1);
        SimTimeInDayBasedOnDiary(idxFreq, idxPreset) ...
            = extractSimTimeInDayFromDiary( ...
            pathToSimDiary, numOfEffeTowers);

        % Host name of the machine on which the simulation was carried out.
        simMachineNames{idxFreq, idxPreset} ...
            = parseCacheFilename(cacheFile.name);
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

% Hostnames and worker numbers.
writecell(simMachineNames, ...
    fullfile(pathToSaveResults, 'SimTimeHostname.csv'));
writematrix(numOfWorkers, ...
    fullfile(pathToSaveResults, 'SimTimeWorkerNums.csv'));

disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Blockage Status Example Figs (Tipp + WHIN)

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating example blockage status maps ...'])

curFlagGenFigsSilently = false;
curFlagZoomIn = true;

% We will use the 1900 MHz results.
curPresets = PRESETS(1:2);
numOfPresets = length(curPresets);
freqInMhz = 1900;
hightsToInspect = [1.5, 10, 100];
legendsShown = false;
for idxPreset = 1:numOfPresets
    preset = PRESETS{idxPreset};
    disp(['        [', datestr(now, datetimeFormat), ...
        '] Processing ', preset, ' (preset #', num2str(idxPreset), ...
        '/', num2str(numOfPresets), ') ...'])

    pathToReadResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', preset, ...
        '_Carrier_', num2str(freqInMhz), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_TO_USE]);

    load(fullfile(pathToReadResults, 'simConfigs.mat'));
    load(fullfile(pathToReadResults, 'simState.mat'));

    [effeCellAntLats, effeCellAntLons] ...
        = simConfigs.utm2deg_speZone( ...
        simState.CellAntsXyhEffective(:, 1), ...
        simState.CellAntsXyhEffective(:, 2));
    effeCellAntLonLats = [effeCellAntLons, effeCellAntLats];

    mapGridLonLats = simState.mapGridLatLonPts(:, [2,1]);

    % Smaller maps for publication.
    if ~isfield(simState,'flagResizeFigForPublication')
        % By default, we do not need to resize figures.
        evalin('base', 'flagResizeFigForPublication = false;');
    else
        evalin('base', ...
            ['flagResizeFigForPublication=', ...
            'simState.RESIZE_FIG_FOR_PUBLICATION']);
    end
    % This works for Tipp and TippExtended.
    flagResizeFigForPublication ...
        = evalin('base', 'flagResizeFigForPublication');
    if flagResizeFigForPublication
        % [500, 500].*0.6 was used for the ICC 2020 paper.
        customFigSize = [500, 500].*0.75;
    else
        defaultFigPos = get(0, 'defaultfigureposition');
        customFigSize = defaultFigPos(3:4);
    end

    for idxH = [1, 2, 11]
        assert(ismember( ...
            simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH), ...
            hightsToInspect), ...
            'idxH is wrong according to hightsToInspect!');
        rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);

        [ hCurBlMap ] ...
            = plotBlockageMap( ...
            [mapGridLonLats, simState.blockageMaps{idxH}], ...
            effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
            curFlagZoomIn, customFigSize);
        hLeg = findobj(hCurBlMap, 'Type', 'Legend');
        switch lower(simConfigs.CURRENT_SIMULATION_TAG)
            case 'tipp'
                set(hLeg, 'Position', [0.6119, 0.1094, 0.2922, 0.1262]);
            case 'shrinkedwhin'
                set(hLeg, 'Position', [0.4809, 0.1104, 0.4258, 0.1119]);
        end
        % Tighten the figure.
        xlabel(''); ylabel('');
        tightfig(hCurBlMap);

        % transparentizeCurLegends;

        % Hide the legend except for the first figure.
        if legendsShown
            legend off;
        else
            legendsShown = true;
        end

        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['BlockageStatusMap_', simConfigs.CURRENT_SIMULATION_TAG, ...
            '_RxHeight_', strrep(num2str(rxAntH), '.', '_')]);
        saveas(hCurBlMap, [pathToSaveFig, '.eps'], 'epsc');
        saveas(hCurBlMap, [pathToSaveFig, '.png']);
        saveas(hCurBlMap, [pathToSaveFig, '.fig']);
    end
end

close all;
disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Cleanup

diary off;

% EOF