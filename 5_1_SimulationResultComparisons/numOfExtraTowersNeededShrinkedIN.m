% NUMOFEXTRATOWERSNEEDEDSHRINKEDIN Estimate the number of towers needed to
% cover the remaining part of the area of interest (for preset ShrinkedIN),
% with the users located at a 1.5 m height.
%
% Yaguang Zhang, Purdue, 07/15/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Parameters

pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
% Ref: https://en.wikipedia.org/wiki/Indiana
areaInSqKmIN = 94321;

ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs', ...
    'NtiaLayoutMergedWithHifldCellTs_Threshold_1000m_LatLonH.csv');

% Presets.
PRESET = 'ShrinkedIN';
% Carrier frequencies.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%     For broadband wireless backhaul.
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CARRIER_FREQUENCIES_IN_MHZ = {1900, 3700, 4700, 7000, 13000, 28000};
% Lidar data set used in the simulations.
LIDAR_DATA_SET_USED = 'IN_DSM_2019';

% Max allowed path loss.
maxAllowedPathLossInDb = 140;

% Expected user height.
expectedRxHInM = 1.5;

% 'FSPL' or 'eHata'
cellSizeEstiMethod = 'eHata';

% The absolute paths to simulation results.
numOfCarriers = length(CARRIER_FREQUENCIES_IN_MHZ);
absPathsToSimResults = cell(numOfCarriers, 1);

for idxCarrier = 1:numOfCarriers
    absPathsToSimResults{idxCarrier} ...
        = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', ['Simulation_', PRESET, ...
        '_Carrier_', num2str(CARRIER_FREQUENCIES_IN_MHZ{idxCarrier}), ...
        'MHz_LiDAR_', LIDAR_DATA_SET_USED]);
end

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'Simulation_Comparisons');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

% The absolute path to save the summary table of the results.
dirToSaveSummaryTable = fullfile(pathToSaveResults, ...
    ['numOfExTs_', PRESET, ...
    '_Lidar_', LIDAR_DATA_SET_USED, ...
    '_maxPL_', num2str(maxAllowedPathLossInDb), ...
    'dB_CellRangeBy_', cellSizeEstiMethod, '.csv']);

% The absolute path to save an illustration plot of the area of interest.
dirToSaveIllu = fullfile(pathToSaveResults, ...
    ['areaOfInt_', PRESET]);

%% Get Functions for Pathloss Value Evaluation
% This is for estimating cell size.
switch cellSizeEstiMethod
    case 'FSPL'
        % Free Space Path Loss (FSPL) inverse.
        fsplInDb2distInM = @(fInMz, lossInDb) 10.^(lossInDb./20) ...
            .*(physconst('LightSpeed')./(fInMz.*(10^6))) ./ 4 ./ pi;
    case 'eHata'
        % Load the NTIA eHata library first, if necessary, to avoid the
        % "unable to find ehata" error. Note that a compiler is needed:
        %    https://www.mathworks.com/support/requirements/supported-compilers.html
        if ~libisloaded('ehata')
            loadlibrary('ehata');
        end

        % The distance values to inspect for creating the cell radius vs
        % path loss data set.
        distsInMToInspect = 1:20000;

        expectedTowerHInM = 50;
    otherwise
        error(['Unsupported cell size estimation method ', ...
            cellSizeEstiMethod, '!']);
end

%% Load Data

% GPS boundary for IN.
inBoundary = load(fullfile(pathToStateInfoFolder, ...
    'IN', 'boundary.mat'));
estAreaInSqKmIN = polyarea( ...
    inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2))./(10^6);

% Cell tower locs.
cellAntsLatLonH = csvread(ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1);
% Preassign memory for the UTM version of cell tower records.
numAnts = size(cellAntsLatLonH, 1);
cellAntsXYH = nan(numAnts, 3);

%% Compute Uncovered Area for Each Simulation

presets = cell(numOfCarriers, 1);
[areasInSqKm, coveredAreasInSqKm, ...
    coveredRadiiPerCellInKm, coveredAreasPerCellInSqKm, ...
    numsOfExistingTs, numsOfExtraTs] ...
    = deal(nan(numOfCarriers, 1));
totalAreaOfIndianaInSqKm = ones(numOfCarriers, 1).*areaInSqKmIN;
maxAllowedPLsInDb = ones(numOfCarriers, 1).*maxAllowedPathLossInDb;
for idxCarrier = 1:numOfCarriers
    presets{idxCarrier} = PRESET;
    curAbsPathToSimResMat = fullfile(absPathsToSimResults{idxCarrier}, ...
        'simState.mat');
    curAbsPathToSimConfMat = fullfile(absPathsToSimResults{idxCarrier}, ...
        'simConfigs_Raw.mat');
    clearvars simState simConfig;
    load(curAbsPathToSimResMat);
    load(curAbsPathToSimConfMat);

    % For GPS and UTM conversions.
    [deg2utm_speZone, utm2deg_speZone] ...
        = genUtmConvertersForFixedZone(simConfigs.UTM_ZONE);

    % Compute the area for the region of interest considering the error
    % caused by reference system conversion.
    areasInSqKm(idxCarrier) = polyarea( ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2))./(10^6) ...
        ./estAreaInSqKmIN.*areaInSqKmIN;

    assert( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(1)==expectedRxHInM, ...
        ['The first RX ant height expected (', ...
        num2str(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(1)), ...
        ' m) is not ', num2str(expectedRxHInM), ' m!']);
    coveredAreasInSqKm(idxCarrier) = ...
        sum(simState.coverageMaps{1}<=maxAllowedPathLossInDb) ...
        ./length(simState.coverageMaps{1}).*areasInSqKm(idxCarrier);

    switch cellSizeEstiMethod
        case 'FSPL'
            coveredRadiiPerCellInKm(idxCarrier) ...
                = fsplInDb2distInM(CARRIER_FREQUENCIES_IN_MHZ{idxCarrier}, ...
                maxAllowedPathLossInDb)/(10^3);
        case 'eHata'
            numOfDistToInspect = length(distsInMToInspect);
            preCompPathLossRecords = nan(numOfDistToInspect, 1);
            for idxD = 1:numOfDistToInspect
                curDInM = distsInMToInspect(idxD);
                preCompPathLossRecords(idxD) = ...
                    computeCoveragePL([0,0,expectedTowerHInM], ...
                    [curDInM,0,expectedRxHInM], ...
                    [9, curDInM/9, zeros(1,10)], simConfigs, zeros(1,3));
            end
            coveredRadiiPerCellInKm(idxCarrier) = ...
                interp1(preCompPathLossRecords, ...
                distsInMToInspect, maxAllowedPathLossInDb)/1000;
    end

    coveredAreasPerCellInSqKm(idxCarrier) ...
        = pi*((coveredRadiiPerCellInKm(idxCarrier))^2);

    [cellAntsXYH(:,1), cellAntsXYH(:,2)] ...
        = deg2utm_speZone(cellAntsLatLonH(:,1), cellAntsLatLonH(:,2));

    boolsTsInAreaOfInt = inpoly2(cellAntsXYH(:,1:2), ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST);
    numsOfExistingTs(idxCarrier) = sum(boolsTsInAreaOfInt);

    numsOfExtraTs(idxCarrier) = ...
        (areasInSqKm(idxCarrier) ...
        - coveredAreasInSqKm(idxCarrier)) ...
        /coveredAreasPerCellInSqKm(idxCarrier);

    % Only need to generate the overview plot once.
    if ~exist([dirToSaveIllu, '.jpg'], 'file')
        % Keep only the cell towers which can cover some part of the area
        % of interest.
        utmXYBoundaryPolyToKeepCellTowers ...
            = polybuffer( ...
            polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
            simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);

        assert(utmXYBoundaryPolyToKeepCellTowers.NumRegions==1, ...
            'The extended area of interest have more than one region!');

        utmXYBoundaryToKeepCellTowers = ...
            utmXYBoundaryPolyToKeepCellTowers.Vertices;
        boolsCellAntsToKeep = inpoly2(cellAntsXYH(:,1:2), ...
            utmXYBoundaryToKeepCellTowers);

        % Effective cellular towers.
        effeCellAntsXYH = cellAntsXYH(boolsCellAntsToKeep, :);
        inEffeCellAntsXYH = cellAntsXYH(~boolsCellAntsToKeep, :);

        [effeCellAntsLats, effeCellAntsLons] ...
            = simConfigs.utm2deg_speZone(effeCellAntsXYH(:,1), ...
            effeCellAntsXYH(:,2));
        [gpsLatsBoundaryToKeepCellTowers, gpsLonsBoundaryToKeepCellTowers] ...
            = simConfigs.utm2deg_speZone( ...
            utmXYBoundaryToKeepCellTowers(:,1), ...
            utmXYBoundaryToKeepCellTowers(:,2));
        [gpsLatsBoundaryOfInterest, gpsLonsBoundaryOfInterest] ...
            = simConfigs.utm2deg_speZone( ...
            simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
            simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
        [inEffeCellAntsLats, inEffeCellAntsLons] ...
            = simConfigs.utm2deg_speZone(inEffeCellAntsXYH(:,1), ...
            inEffeCellAntsXYH(:,2));

        % Convert IN boundary in the UTM system to GPS.
        assert( ...
            strcmp(inBoundary.boundary.UTM_ZONE, simConfigs.UTM_ZONE), ...
            ['UTM zone used in the simulation is different ', ...
            'from that for the IN boundary!']);
        [inBoundaryLats, inBoundaryLons] = ...
            utm2deg_speZone( ...
            inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
            inBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

        curCustomFigSize = [500, 500];
        % For plotting.
        areaOfInterestColor = [0.9290 0.6940 0.1250];
        lightBlue = [0.3010 0.7450 0.9330];
        darkBlue = [0 0.4470 0.7410];
        colorEffectiveTowers = 'b';
        markerEffectiveTowers = '.';
        markerSizeEffectiveTowers = 12;
        colorIneffectiveTowers = lightBlue;
        markerIneffectiveTowers = '.';
        lineWidthIneffectiveTowers = 1.5;

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
        hAreaOfInterest = plot(polyshape( ...
            [gpsLonsBoundaryOfInterest, gpsLatsBoundaryOfInterest]), ...
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
            [min(inBoundaryLons), ...
            max(inBoundaryLons), ...
            min(inBoundaryLats), ...
            max(inBoundaryLats)], extensionFactor, simConfigs);
        adjustFigSizeByContent(hFigCellOverview, axisLonLatToSet, ...
            'height', weightForWidth.*0.9);
        view(2); plot_google_map;
        grid on; grid minor;
        xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
        box on; makescale('sw', 'units', 'si');
        hLeg = legend( ...
            [hAreaOfInterest, hEffeCells], ...
            'Simulation area', 'Effective towers', 'AutoUpdate','off');
        if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'ExtendedTipp')
            % Manually adjust the figure for publication.
            set(hLeg, 'Position', [0.4043, 0.7667, 0.5091, 0.1680]);
            axis([-88.1448722958386, -85.5528198028410, ...
                39.3771676547221, 41.7336964261903]);
        elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
            set(hLeg, 'Location', 'SouthEast');
            hMS = makescale('nw', 'units', 'si');
            hMS(3).FontSize = 9;
            hMS(3).FontWeight = 'bold';
        end
        hPolyIn = plot3(inBoundaryLons, inBoundaryLats, ...
            ones(size(inBoundaryLons)), ...
            'r-.', 'LineWidth', 3);

        saveas(hFigCellOverview, [dirToSaveIllu, '.fig']);
        saveas(hFigCellOverview, [dirToSaveIllu, '.jpg']);
        saveas(hFigCellOverview, [dirToSaveIllu, '.eps'], 'epsc');

        set(hFigCellOverview, 'Color', 'w');
        export_fig(hFigCellOverview, [dirToSaveIllu, '_alt.eps']);
    end
end

%% Put Results in a Table

carrierFreqInMHz = CARRIER_FREQUENCIES_IN_MHZ';
summaryTable = table(presets, areasInSqKm, totalAreaOfIndianaInSqKm, ...
    carrierFreqInMHz, ...
    maxAllowedPLsInDb, coveredAreasInSqKm, ...
    coveredRadiiPerCellInKm, coveredAreasPerCellInSqKm, ...
    numsOfExtraTs, numsOfExistingTs);
disp(summaryTable);

% Export the results to a file.
writetable(summaryTable, dirToSaveSummaryTable);

% Export the results to a file as, roughly speaking, a Latex snippet. Note
% that the column order here is adjusted to match the Asilomar paper.
summaryTableRounded = table(presets, ...
    round(totalAreaOfIndianaInSqKm, 0), ...
    round(areasInSqKm, 0), ...
    carrierFreqInMHz, maxAllowedPLsInDb, ...
    round(coveredAreasInSqKm, 0), ...
    round(coveredRadiiPerCellInKm, 1), ...
    round(coveredAreasPerCellInSqKm, 1), ...
    ceil(numsOfExtraTs), ...
    numsOfExistingTs);
disp(summaryTableRounded);

% Export new values to a file. Note that Delimiter must be one of these
% characters: ' ', '\t', ',', ';', '|', or its corresponding character
% name: 'space', 'tab', 'comma', 'semi', or 'bar'.
[logFilePath, logFileName, ~] = fileparts(dirToSaveSummaryTable);
writetable(summaryTableRounded, ...
    [fullfile(logFilePath, logFileName), '_Latex.csv']);

% EOF