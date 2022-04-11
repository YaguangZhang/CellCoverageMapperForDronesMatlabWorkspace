function [ simState ] ...
    = simulateCoverage(pathToSaveResults, ...
    lidarFileAbsDirs, lidarFileXYCoveragePolyshapes, ...
    cellAntsXYH, simConfigs)
%SIMULATECOVERAGE Generate coverage and blockage maps via channel
%simulations.
%
% This function is based on:
%   lib/analyzeCoverageViaChannelSimulation.m
%
% Coverage map:
%   - For locations with RX-to-TX distance>=1km
%     We use the eHata model (NTIA C++ implementation).
%	- For locations with RX-to-TX distance<1km.
%     We use eHata weighted with line-of-sight (determined by partial
%     clearance of 1st Fresnel Zone) FSPL.
%
% Blockage map:
%   Line-of-sight (determined by partial clearance of 1st Fresnel Zone)
%   FSPL is used.
%
% Inputs:
%   - pathToSaveResults
%     The absolute directory to save resulting .mat files and plots.
%   - lidarFileAbsDirs, lidarFileXYCoveragePolyshapes
%     The outputs of the function preprocessIndianaLidarDataSet, which are
%     two cell arrays holding the absolute directories to all the LiDAR
%     data available as raw .img/.tif files, and their polygon boundries in
%     the UTM (x, y) system.
%   - cellAntsXYH
%     The (x, y, height) of the cellular towers. Note: we use "height" to
%     indicate the vertical distance from the ground to the antenna;
%     "elevation" to indicate the ground elevation; and "altitude" to
%     indicate elevation+height.
%   - simConfigs
%     A structure for the configuration parameters for the simulation,
%     including fields
%       - CURRENT_SIMULATION_TAG
%         A string label to identify this simulation.
%       - UTM_ZONE
%         The zone label to use in the UTM (x, y) system, e.g. '16 T'.
%       - UTM_X_Y_BOUNDARY_OF_INTEREST
%         The UTM (x, y) polygon representing the area of interest for
%         generating the coverage maps.
%       - NUM_OF_PIXELS_FOR_LONGER_SIDE
%         We will use this number of pixels for the longer side
%         (width/height) of the map; the number of pixels for the other
%         side will be proportional to its length.
%       - MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M
%         The guaranteed spacial resolution for terrain profiles; a larger
%         value will decrease the simulation time but small obstacles may
%         get ingored.
%       - MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M
%         Similarly, the guaranteed spacial resolution for LiDAR profiles.
%       - MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE
%         The guaranteed minimum number of LiDAR z (or possibly elevation)
%         elements in one terrain profile; this will ensure non-empty
%         terrain profiles.
%       - deg2utm_speZone, utm2deg_speZone
%         The functions to convert GPS (lat, lon) to (and back from) UTM
%         (x, y) for the LiDAR data area.
%
% Output:
%   - simState
%     A structure for the state of the simulation, including fields
%       - mapGridXYPts, mapGridLatLonPts
%         The UTM (x,y) and GPS (lat, lon) locations for the map grid
%         points, respectively. Note that the shape of the map is
%         determined by simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST. This means
%         that depending on the specificed area of interest, the shape of
%         all the grid points may be arbitrary, so all the grid point
%         information, like locations mapGridXYPts and mapGridLatLonPts,
%         are stored sequentially.
%       - blockageMaps, coverageMaps
%         Two cells holding the final aggregated path loss maps, for the
%         blockage (status) maps and the coverage (path loss) maps,
%         respectively, for all the RX height that have been inspected.
%       - blockageDistMaps, blockageByTerrainDistMaps
%         Two cells holding the final aggregated maps for the accumulated
%         blockage distances evaluated with LiDAR data and ground elevation
%         data, respectively.
%       - blockageByVegDistMaps
%         A cell stores the accumulative blocage distances caused by
%         vegetation that are estimated by comparing
%         simState.blockageDistMapsForEachCell and
%         simState.blockageByTerrainDistMapsForEachCell. Note that the
%         element values could be negative.
%       - pathLossWithVegMaps
%         A cell holding the final aggregated path loss map considering
%         vegetation. It is gotten from the NTIA eHata model + ITU
%         obtruction by woodland model. Ref: Recommendation ITU-R P.833-7
%         (02/2012) Attenuation in vegetation.
%       - blockageMapsForEachCell, coverageMapsForEachCell
%         Two cells for the blockage and coverage path loss maps (for all
%         cell antennas affecting the simulation for all drone heights),
%         respectively. For example,
%               simState.blockageMapsForEachCell{idxCell}{idxDroneHeight}
%         is the path loss vector (for locations in simState.mapGridXYPts)
%         for the idxCell-th cellular tower that has effect on the
%         simulation, at the idxDroneHeight-th drone height that needs to
%         be inspected.
%       - blockageDistMapsForEachCell, blockageByTerrainDistMapsForEachCell
%         Cells holding the (accumulative) blockage distance in meters with
%         a structure similar to that for blockageMapsForEachCell. The
%         blockage is evaluated with the LiDAR DSM data and the ground
%         elevation data, respectively.
%       - pathLossWithVegMapsForEachCell
%         A cell holding the path loss values considering the propogation
%         via vegetation with a structure similar to that for
%         blockageMapsForEachCell.
%       - TimeUsedInSForEachPixel
%         For debugging and evaluating computation performance of the
%         simulator, we also record the time used in second for each map
%         pixel. For example,
%               simState.TimeUsedInSForEachPixel{idxCell}{idxDroneHeight}
%         is the processing time vector (for locations in
%         simState.mapGridXYPts) for the idxCell-th cellular tower that has
%         effect on the simulation, at the idxDroneHeight-th drone height
%         that needs to be inspected.
%
%         We expect the time used for a pixel at its first drone height
%         will be way longer than for other heights, because for the other
%         heights, we will reuse the terrain profiles generated for the
%         first height.
%       - CellAntsXyhAll, CellAntsXyhEffective
%         For debugging and plotting, we record the location information
%         for (1) all the cellular towers and (2) those that are close
%         enough to the area of interest to affect the final results.
%       - CellAntsEffectiveIds
%         The indices for the effective cellular towers, used as their
%         unique IDs in the simulation.
%       - MaxAllowedPathlossInDb
%         For plotting, we need an expected value for the maximum allowed
%         path loss. With 20 MHz bandwidth, 9 dB RX noise figure, and 100 W
%         TX power, we have a maximum path loss of ~142 dB.
%
% Yaguang Zhang, Purdue, 05/07/2021

%% Before Processing the Data

% Close the pool if it exists to make sure all resource/workers are
% available for the simulation.
previousPool = gcp('nocreate');
if ~isempty(previousPool)
    delete(previousPool);
end

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% Extra information on the LiDAR data set.
%   - Overall boundries for the area covered by the LiDAR data set in UTM.
% lidarFilesXYCoveragePolyshape ...
%     = mergePolygonsForAreaOfInterest(lidarFileXYCoveragePolyshapes, 1);
%   - Centroids for the LiDAR files in UTM.
lidarFileXYCentroids ...
    = extractCentroidsFrom2DPolyCell(lidarFileXYCoveragePolyshapes);
%   - The .mat copies for the LiDAR data. For the 2019 dataset, they are
%   stored in a cache folder.
lidarMatFileAbsDirs = lidarFileAbsDirs;
for idxMatF = 1:length(lidarMatFileAbsDirs)
    [lidarMatFPath, lidarMatFName, ~] ...
        = fileparts(lidarMatFileAbsDirs{idxMatF});
    lidarMatFileAbsDirs{idxMatF} = fullfile(lidarMatFPath, '..', ...
        'MatlabCache', [lidarMatFName, '.mat']);
end

% We have not yet worked with any area of interest with multiple regions.
assert(length(regions( ...
    polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST)))==1, ...
    'The area of interest should have only one region!');

% For plotting.
areaOfInterestColor = [0.9290 0.6940 0.1250];
lightBlue = [0.3010 0.7450 0.9330]; %#ok<NASGU>
darkBlue = [0 0.4470 0.7410];
colorEffectiveTowers = 'b';
markerEffectiveTowers = '.';
markerSizeEffectiveTowers = 12;
colorIneffectiveTowers = 'r';
markerIneffectiveTowers = 'x';
lineWidthIneffectiveTowers = 1.5;

if simConfigs.RESIZE_FIG_FOR_PUBLICATION
    customFigSize = [500, 500];
else
    defaultFigPos = get(0,'defaultfigureposition');
    customFigSize = defaultFigPos(3:4);
end

% For each cell, if this ratio (or more) of workers will not be assigned
% with any pixels to process, normal for loop (instead of parfor) will be
% used. For example, set it to be 0 to always use for; set it to 1 to
% always use parfor.
MIN_RATIO_OF_EMPTY_WORKERS_TO_AVOID_PAR = 1;

% Define the ITU obstruction by woodland model.
defineModelItuObsByWoodland;
ituGamma = modelItuObsByWoodland.estimateGamma( ...
    simConfigs.CARRIER_FREQUENCY_IN_MHZ);
ituAm = modelItuObsByWoodland.estimateAm( ...
    simConfigs.CARRIER_FREQUENCY_IN_MHZ);
estimateExcessPL = @(dInM) modelItuObsByWoodland.excessLossFormula( ...
    max(dInM, 0), ituGamma, ituAm);

%% Create the Map Grid

switch simConfigs.CURRENT_SIMULATION_TAG
    case 'WHIN_WEATHER_STATIONS'
        [mapGridLats, mapGridLons] ...
            = loadWhinWeatherStationInfo();
        simState.mapGridLatLonPts = [mapGridLats, mapGridLons];

        % Convert (lat, lon) to UTM (x, y).
        [mapGridXs, mapGridYs] = simConfigs.deg2utm_speZone( ...
            simState.mapGridLatLonPts(:,1), ...
            simState.mapGridLatLonPts(:,2));
        simState.mapGridXYPts = [mapGridXs, mapGridYs];
    otherwise
        [simState.mapGridXYPts, ] ...
            = buildSimGrid(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ...
            simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE);

        % Convert UTM (x, y) to (lat, lon).
        [mapGridLats, mapGridLons] = simConfigs.utm2deg_speZone( ...
            simState.mapGridXYPts(:,1), simState.mapGridXYPts(:,2));
        simState.mapGridLatLonPts = [mapGridLats, mapGridLons];
end

% Plot. Only generate figures if they are not there.
curDirToSave = fullfile(pathToSaveResults, 'Overview_RxLocGrid');
if ~exist([curDirToSave, '.eps'], 'file')
    if simConfigs.RESIZE_FIG_FOR_PUBLICATION
        curCustomFigSize = customFigSize.*0.6;
    else
        curCustomFigSize = customFigSize;
    end
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
        'Area of interest', 'UAV location grid');
    xlabel('UTM x (m)'); ylabel('UTM y (m)');
    % Adjust legend and the exponent label for y axis.
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'tipp'
            annotation(hFigAreaOfInterest, 'textbox',...
                [0.1409 0.8502 0.1785 0.0945],...
                'String', ['\times10^', ...
                num2str(hCurAxis.YAxis.Exponent)],...
                'FontWeight', hCurAxis.YAxis.FontWeight,...
                'FontSize', hCurAxis.YAxis.FontSize,...
                'EdgeColor', 'none');
            yticks(yticks);
            yticklabels(yticklabels);
        case 'extendedtipp'
            set(hLeg, 'Location', 'northwest');
            transparentizeCurLegends;

            annotation(hFigAreaOfInterest, 'textbox',...
                [0.1609 0.8369 0.1785 0.0945],...
                'String', ['\times10^', ...
                num2str(hCurAxis.YAxis.Exponent)],...
                'FontWeight', hCurAxis.YAxis.FontWeight,...
                'FontSize', hCurAxis.YAxis.FontSize,...
                'EdgeColor', 'none');
            yticks(yticks);
            yticklabels(yticklabels);
        case 'shrinkedin'
            set(hLeg, 'Location', 'SouthEast');
    end

    saveEpsFigForPaper(hFigAreaOfInterest, curDirToSave);
end

%% Load Antenna Info

% Keep only the cell towers which can cover some part of the area of
% interest.
utmXYBoundaryToKeepCellTowers = extendUtmXYBoundOfInt( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ...
    simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);

boolsCellAntsToKeep ...
    = InPolygon(cellAntsXYH(:,1), cellAntsXYH(:,2), ...
    utmXYBoundaryToKeepCellTowers(:,1), ...
    utmXYBoundaryToKeepCellTowers(:,2));
% Effective cellular towers.
effeCellAntsXYH = cellAntsXYH(boolsCellAntsToKeep, :);
inEffeCellAntsXYH = cellAntsXYH(~boolsCellAntsToKeep, :);

% Save information for cells.
simState.CellAntsXyhAll = cellAntsXYH;
simState.CellAntsXyhEffective = effeCellAntsXYH;
simState.CellAntsEffectiveIds = find(boolsCellAntsToKeep);

curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider.png');
if ~exist(curDirToSave, 'file')
    % Overview plot in UTM (x, y).
    hFigCellOverview = figure; hold on; set(gca, 'fontWeight', 'bold');
    hEffeCells = plot(effeCellAntsXYH(:,1), ...
        effeCellAntsXYH(:,2), markerEffectiveTowers, ...
        'MarkerSize', markerSizeEffectiveTowers, ...
        'Color', colorEffectiveTowers);
    hExtendedArea = plot(polyshape([utmXYBoundaryToKeepCellTowers; ...
        nan nan; ...
        simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST]));
    hAreaOfInterest = plot( ...
        polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
        'FaceColor', areaOfInterestColor);
    axis equal; view(2); grid on; grid minor;
    axisXYToSet = [min(utmXYBoundaryToKeepCellTowers(:,1)), ...
        max(utmXYBoundaryToKeepCellTowers(:,1)), ...
        min(utmXYBoundaryToKeepCellTowers(:,2)), ...
        max(utmXYBoundaryToKeepCellTowers(:,2))];
    axisXYToSet = extendAxisByFactor(axisXYToSet, 0.2);
    hIneffeCells = plot(inEffeCellAntsXYH(:,1), ...
        inEffeCellAntsXYH(:,2), markerIneffectiveTowers, ...
        'Color', colorIneffectiveTowers, ...
        'LineWidth', lineWidthIneffectiveTowers);
    uistack(hIneffeCells, 'bottom'); uistack(hEffeCells, 'top');
    adjustFigSizeByContent(hFigCellOverview, axisXYToSet, 'height');
    xlabel('UTM x (m)'); ylabel('UTM y (m)'); box on;
    hLeg = legend([hAreaOfInterest, hExtendedArea, hEffeCells, hIneffeCells], ...
        'Area of interest', 'Extended area', 'Cell towers to consider', ...
        'Ineffective cell towers');
    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
        set(hLeg, 'Location', 'SouthEast');
    end
    saveas(hFigCellOverview, curDirToSave);
end

% Overview plot in GPS (lon, lat).
curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider_RoadMap');
if ~exist([curDirToSave, '.eps'], 'file')
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

    if simConfigs.RESIZE_FIG_FOR_PUBLICATION
        curCustomFigSize = customFigSize.*0.8;
    else
        curCustomFigSize = customFigSize;
    end
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
        [min(gpsLonsBoundaryToKeepCellTowers), ...
        max(gpsLonsBoundaryToKeepCellTowers), ...
        min(gpsLatsBoundaryToKeepCellTowers), ...
        max(gpsLatsBoundaryToKeepCellTowers)], extensionFactor, simConfigs);
    adjustFigSizeByContent(hFigCellOverview, axisLonLatToSet, ...
        'height', weightForWidth.*0.9);
    view(2); plot_google_map;
    grid on; grid minor;
    xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
    box on; makescale('sw', 'units', 'si');
    hLeg = legend( ...
        [hAreaOfInterest, hExtendedArea, hEffeCells, hIneffeCells], ...
        'Area of interest', 'Extended area', 'Cell towers to consider', ...
        'Ineffective cell towers', 'defaultLegendAutoUpdate','off');
    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'ExtendedTipp')
        % Manually adjust the figure for publication.
        set(hLeg, 'Position', [0.4043, 0.7667, 0.5091, 0.1680]);
        axis([-88.1448722958386, -85.5528198028410, ...
            39.3771676547221, 41.7336964261903]);
    elseif strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
        set(hLeg, 'Location', 'SouthEast');
        makescale('nw', 'units', 'si');
    end
    saveEpsFigForPaper(hFigCellOverview, curDirToSave);

    plot_google_map('MapType', 'satellite');
    curDirToSave = fullfile(pathToSaveResults, ...
        'Overview_CellularTowersToConsider_SatelliteMap');
    saveEpsFigForPaper(hFigCellOverview, curDirToSave);
end

%% Simulation

disp(' ')
disp('    Computing path losses ...')
disp(' ')
disp('        Closing figures to save RAM ...')
close all;
disp('        Done!')

disp(' ')
disp('    Initializing simulation ...')

% Cache file for recording the computation process.
[~, hostname] = system('hostname');
hostname = strtrim(hostname);
pathToCache = fullfile(pathToSaveResults, ...
    ['CovAnalysisCache_Task_', ...
    simConfigs.CURRENT_SIMULATION_TAG, ...
    '_LidarSet_', simConfigs.LIDAR_DATA_SET_TO_USE, ...
    '_NumPix_', ...
    num2str(simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE), ...
    '_TerProRes_', ...
    num2str(simConfigs ...
    .MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M), ...
    '_LidProRes_', ...
    num2str(simConfigs ...
    .MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M), ...
    '_HostName_', hostname, '.mat']);

if exist(pathToCache, 'file')
    disp(' ')
    disp('        Loading cached history path loss results ...')

    load(pathToCache); %#ok<LOAD>

    disp('        Done!')
end

% For initializing results.
[numOfEffeCellAnts, ~] = size(effeCellAntsXYH);
numOfRxHeightToInspect ...
    = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);
[numPixelsPerMap, ~] = size(simState.mapGridXYPts);
if ~isfield(simState, 'blockageMapsForEachCell')
    [simState.blockageMapsForEachCell, ...
        simState.coverageMapsForEachCell, ...
        simState.blockageDistMapsForEachCell, ...
        simState.blockageByTerrainDistMapsForEachCell, ...
        simState.pathLossWithVegMapsForEachCell, ...
        simState.TimeUsedInSForEachPixel] ...
        = deal(cell(numOfEffeCellAnts, 1));

    % Progress monitor (proMon): for reporting progress and estimating the
    % remaining simulation time.
    proMon.floatFomatter = '%.2f';
    proMon.pixCnt = 0;
    proMon.numOfPixToProcess = 0;

    % For each effective cellular tower, find grid pixels that are within
    % the coverage range and pre-assign them to workers. Also, initialize
    % result cells and count the total number of pixels to be processed for
    % progress monitoring.
    locIndicesForAllWorkersForAllCellsEff = cell(numOfEffeCellAnts,1);
    for idxEffeCellAnt = 1:numOfEffeCellAnts
        % Cellular location.
        curCellXYH = effeCellAntsXYH(idxEffeCellAnt, :);

        % Find drone locations that are within the current cellular tower's
        % coverage range in the UTM (x,y) system.
        curTxToRxDistsInM ...
            = vecnorm(simState.mapGridXYPts-curCellXYH(1:2), 2, 2);
        curIndicesRxLocsToConsider ...
            = find(curTxToRxDistsInM ...
            < simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);

        % We will use multiple workers to churn through the drone
        % locations. To avoid repeative environment set up and data
        % transfer, here we pre-assign the pixels to be processed by each
        % worker to avoid unnecessary data copying during worker
        % initialization.
        locIndicesForAllWorkers ...
            = preassignTaskIndicesToWorkers( ...
            length(curIndicesRxLocsToConsider));
        locIndicesForAllWorkersForAllCellsEff{idxEffeCellAnt} ...
            = cellfun(@(indices) ...
            curIndicesRxLocsToConsider(indices)', ...
            locIndicesForAllWorkers, 'UniformOutput', false);

        % Initialize results.
        for idxDroneHeightInM = 1:numOfRxHeightToInspect
            simState.blockageMapsForEachCell{idxEffeCellAnt} ...
                {idxDroneHeightInM} = nan(1, numPixelsPerMap);
            simState.coverageMapsForEachCell{idxEffeCellAnt} ...
                {idxDroneHeightInM} = nan(1, numPixelsPerMap);
            simState.blockageDistMapsForEachCell{idxEffeCellAnt} ...
                {idxDroneHeightInM} = nan(1, numPixelsPerMap);
            simState.blockageByTerrainDistMapsForEachCell{ ...
                idxEffeCellAnt}{idxDroneHeightInM} = ...
                nan(1, numPixelsPerMap);
            simState.pathLossWithVegMapsForEachCell{idxEffeCellAnt} ...
                {idxDroneHeightInM} = nan(1, numPixelsPerMap);
            simState.TimeUsedInSForEachPixel{idxEffeCellAnt} ...
                {idxDroneHeightInM} = zeros(1, numPixelsPerMap);
        end

        % Update the total number of pixels to be processed.
        proMon.numOfPixToProcess = proMon.numOfPixToProcess ...
            + length(curIndicesRxLocsToConsider);
    end

    % Cache the results.
    save(pathToCache, 'simState', 'locIndicesForAllWorkers', ...
        'locIndicesForAllWorkersForAllCellsEff', 'proMon', '-v7.3');

    disp('        Done!')
end

% Suppress this warning in the cluster to get clearer feedbacks from the
% program.
parfevalOnAll(gcp(), ...
    @warning, 0, 'off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

% Get the number of workers available for estimating progress.
curCluster = gcp;
numOfWorkersInCurCluster = curCluster.NumWorkers;

% The simulation progress is recorded in the integer variable
% nextIdxEffeCellAnt, the next idxEffeCellAnt to process.
if ~exist('nextIdxEffeCellAnt', 'var')
    nextIdxEffeCellAnt = 1;
end

mapGridXYPts = simState.mapGridXYPts;
for idxEffeCellAnt = nextIdxEffeCellAnt:numOfEffeCellAnts

    % Check the memory usage and restart the parallel pool if necessary.
    guardMemAvailOnLinux;

    disp(['        [', datestr(now, 'yyyy/mm/dd HH:MM:ss'), ...
        '] Effective cellular tower #', ...
        num2str(idxEffeCellAnt), '/', num2str(numOfEffeCellAnts), ' ...']);

    % Cellular location.
    curCellXYH = effeCellAntsXYH(idxEffeCellAnt, :);

    % Load the task assignment results.
    locIndicesForAllWorkers ...
        = locIndicesForAllWorkersForAllCellsEff{idxEffeCellAnt};
    numOfWorkers = length(locIndicesForAllWorkers);
    assert(numOfWorkers==numOfWorkersInCurCluster, ...
        'Preassigned tasks do not match current number of workers!');

    % Pre-allocate space. To make sure parfor works, we will store all
    % results from one worker into one matrix with each row being:
    % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
    % idxDroneHeightInM, blockageDistInM, blockageByTerrainDistInM].
    resultsFromWorkersCell = cell(numOfWorkers,1);
    for idxWorker = 1:numOfWorkers
        resultsFromWorkersCell{idxWorker} ...
            = nan(length(locIndicesForAllWorkers{idxWorker}) ...
            .*numOfRxHeightToInspect, 7);
    end

    % Utilize parfor only when some workers will get something to process.
    if (sum(cellfun(@(c) isempty(c), locIndicesForAllWorkers))...
            /numOfWorkers)>=MIN_RATIO_OF_EMPTY_WORKERS_TO_AVOID_PAR
        parforArg = 0;
    else
        parforArg = numOfWorkers;
    end

    % To make sure the overhead time is recorded by the first pixel
    % processed by the work.
    curExecTimeStartTic = tic;

    % For recording and estimating processing time in parfor.
    pathToSaveOverheadTimeMats = fullfile(pathToSaveResults, ...
        'ProcessingTimeCacheRecords');
    if exist(pathToSaveOverheadTimeMats, 'dir')
        % Specify 's' to also attempts to remove all subfolders and files,
        % regardless of their write permissions.
        rmdir(pathToSaveOverheadTimeMats, 's');
    end
    mkdir(pathToSaveOverheadTimeMats);
    if parforArg ~= 0
        pathsToOverheadTimeInSecStarts = cell(parforArg, 1);
        for curTaskId = 1:parforArg
            curPathToOverheadTimeInSecStart ...
                = fullfile(pathToSaveOverheadTimeMats, ...
                ['Worker_', num2str(curTaskId), '.mat']);
            save(curPathToOverheadTimeInSecStart, ...
                'curExecTimeStartTic');
            pathsToOverheadTimeInSecStarts{curTaskId} ...
                = curPathToOverheadTimeInSecStart;
        end
    else
        pathsToOverheadTimeInSecStarts = cell(1, 1);
        curPathToOverheadTimeInSecStart ...
            = fullfile(pathToSaveOverheadTimeMats, ...
            'Worker_0.mat');
        save(curPathToOverheadTimeInSecStart, ...
            'curExecTimeStartTic');
        pathsToOverheadTimeInSecStarts{1} ...
            = curPathToOverheadTimeInSecStart;
    end

    parfor (idxWorker = 1:numOfWorkers, parforArg)
        % Processing time considering the overhead.
        if parforArg ~= 0
            curTask = getCurrentTask();
            curTaskId = curTask.ID;
        else
            curTaskId = 1;
        end

        curExecTimeInSecStartMatContent ...
            = load(pathsToOverheadTimeInSecStarts{curTaskId}); %#ok<PFBNS>
        curExecTimeStartTic = curExecTimeInSecStartMatContent ...
            .curExecTimeStartTic;

        % Load the NTIA eHata library first, if necessary, to avoid the
        % "unable to find ehata" error. Note that a compiler is needed:
        %    https://www.mathworks.com/support/requirements/supported-compilers.html
        if ~libisloaded('ehata')
            loadlibrary('ehata');
        end

        % Load our Python module.
        py_addpath(fullfile(pwd, 'lib', 'python'));

        curWorkerPixCnt = 0;
        curDroneLocIndices = locIndicesForAllWorkers{idxWorker};
        curWorkerNumPixs = length(curDroneLocIndices);
        curWorkerNumPixsToReportProgress = ceil(curWorkerNumPixs ...
            .*simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT); %#ok<PFBNS>
        % Make sure curDroneLocIndices is a row vector.
        curDroneLocIndices = curDroneLocIndices(:)';
        for idxDroneLoc = curDroneLocIndices

            % Report progress when necessary.
            if mod(curWorkerPixCnt, ...
                    curWorkerNumPixsToReportProgress)==0
                disp(['        [', datestr(now, 'yyyy/mm/dd HH:MM:ss'), ...
                    '] Worker #', num2str(idxWorker), ...
                    '/', num2str(numOfWorkers), ' (', ...
                    num2str(curWorkerPixCnt/curWorkerNumPixs ...
                    /numOfRxHeightToInspect*100, ...
                    '%.2f'), '%) ...']);
            end

            % Drone location.
            curDroneXY = mapGridXYPts(idxDroneLoc, :); %#ok<PFBNS>

            % Generate terrain and LiDAR profiles. We will use an empty
            % cache path to skip saving the profiles.
            absPathToCacheMatFile = '';
            [curTerrainProfile, curLidarProfile] ...
                = fetchTerrainAndLidarProfiles( ...
                absPathToCacheMatFile, ...
                curCellXYH(:,1:2), ...
                curDroneXY, simConfigs, ...
                lidarFileXYCentroids, ...
                lidarFileXYCoveragePolyshapes, ...
                lidarMatFileAbsDirs); %#ok<PFBNS>

            % Fill NaN elements to the minimum value in that profile.
            numOfNanProfEles ...
                = sum(isnan([curTerrainProfile; curLidarProfile]));
            if numOfNanProfEles>0
                % warning([num2str(numOfNanProfEles), ...
                %     ' profile elements are NaN!']);
                curTerrainProfile(isnan(curTerrainProfile)) ...
                    = min(curTerrainProfile);
                curLidarProfile(isnan(curLidarProfile)) ...
                    = min(curLidarProfile);
            end

            for idxDroneHeightInM = 1:numOfRxHeightToInspect
                % Drone height.
                curDroneH ...
                    = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M( ...
                    idxDroneHeightInM);

                [curBlockagePL, curCoveragePL, ...
                    curBlockageDistInM, curBlockageByTerrainDistInM] ...
                    = computeBlockageAndCoveragePLs( ...
                    curLidarProfile, curTerrainProfile, ...
                    curCellXYH, [curDroneXY, curDroneH], ...
                    simConfigs);

                curExecTime = toc(curExecTimeStartTic);

                % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
                % idxDroneHeightInM, curBlockageDistInM,
                % curBlockageByTerrainDistInM].
                curWorkerPixCnt = curWorkerPixCnt+1;
                resultsFromWorkersCell{idxWorker} ...
                    (curWorkerPixCnt, :) ...
                    = [curBlockagePL, curCoveragePL, curExecTime, ...
                    idxDroneLoc, idxDroneHeightInM, ...
                    curBlockageDistInM, curBlockageByTerrainDistInM];

                % Reset timer.
                curExecTimeStartTic = tic;
                curPathToSaveOverheadTimeInSecStart ...
                    = pathsToOverheadTimeInSecStarts{curTaskId};
                parsave(curPathToSaveOverheadTimeInSecStart, ...
                    curExecTimeStartTic);
            end
        end

        % Final progress report which is expected to be 100% done.
        disp(['        [', datestr(now, 'yyyy/mm/dd HH:MM:ss'), ...
            '] Worker #', num2str(idxWorker), ...
            '/', num2str(numOfWorkers), ' (', ...
            num2str(curWorkerPixCnt/curWorkerNumPixs ...
            /numOfRxHeightToInspect*100, ...
            '%.2f'), '%) ...']);

        % Clean up the library.
        if libisloaded('ehata')
            unloadlibrary('ehata');
        end
    end

    % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
    % idxDroneHeightInM].
    resultsFromWorkersMat = vertcat(resultsFromWorkersCell{:});
    [numOfResultsFromWorkers, ~] = size(resultsFromWorkersMat);
    % Output the results.
    for idxResult = 1:numOfResultsFromWorkers
        curResult = resultsFromWorkersMat(idxResult, :);

        curBlockagePL = curResult(1);
        curCoveragePL = curResult(2);
        curPixelExecTime = curResult(3);
        curIdxDroneLoc = curResult(4);
        curIdxDroneHeightInM = curResult(5);
        curBlockageDistInM = curResult(6);
        curBlockageByTerrainDistInM = curResult(7);

        % Estimate the vegetable depth.
        curBlockageByVegDistInM = ...
            curBlockageDistInM - curBlockageByTerrainDistInM;
        % Estimate the overall path loss considering the vegetation. Note
        % that estimateExcessPL will set negative inputs to zero.
        curPathLossWithVeg = ...
            curCoveragePL + estimateExcessPL(curBlockageByVegDistInM);

        simState.blockageMapsForEachCell{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = curBlockagePL;
        simState.coverageMapsForEachCell{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = curCoveragePL;
        simState.blockageDistMapsForEachCell{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = curBlockageDistInM;
        simState.blockageByTerrainDistMapsForEachCell{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = ...
            curBlockageByTerrainDistInM;
        simState.pathLossWithVegMapsForEachCell{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = curPathLossWithVeg;
        simState.TimeUsedInSForEachPixel{idxEffeCellAnt} ...
            {curIdxDroneHeightInM}(curIdxDroneLoc) = curPixelExecTime;
    end

    latestEstiPixExecTimeInS = mean(resultsFromWorkersMat(:,3));
    latestEstiTowerExecTimeInS = sum(resultsFromWorkersMat(:,3)) ...
        ./numOfWorkersInCurCluster;
    % Take into account all the heights inspected.
    proMon.pixCnt = proMon.pixCnt ...
        + numOfResultsFromWorkers./numOfRxHeightToInspect;

    disp(' ')
    disp(['        [', datestr(now, 'yyyy/mm/dd HH:MM:ss'), ...
        '] Finished cellular tower #', ...
        num2str(idxEffeCellAnt), '/', num2str(numOfEffeCellAnts)]);
    disp(['            Estimated time used for this tower: ', ...
        num2str(latestEstiTowerExecTimeInS, proMon.floatFomatter), ...
        ' seconds']);
    if ~isnan(latestEstiTowerExecTimeInS)
        disp(['                (', ...
            seconds2human(latestEstiTowerExecTimeInS), ')']);
    end
    disp(['            Total progress: ', ...
        num2str(proMon.pixCnt./proMon.numOfPixToProcess.*100, ...
        proMon.floatFomatter), ...
        '% (Pixel #', num2str(proMon.pixCnt), ...
        '/', num2str(proMon.numOfPixToProcess), ')']);
    remainingTimeInS = latestEstiPixExecTimeInS ...
        .*(proMon.numOfPixToProcess-proMon.pixCnt) ...
        .*numOfRxHeightToInspect ...
        ./numOfWorkersInCurCluster;
    disp(['            Estimated remaining time: ', ...
        num2str(remainingTimeInS, proMon.floatFomatter), ...
        ' seconds']);
    if ~isnan(remainingTimeInS)
        disp(['                (', ...
            seconds2human(remainingTimeInS), ')']);
    end
    disp(' ');

    nextIdxEffeCellAnt = idxEffeCellAnt+1;
    % We will update simState for each tower.
    save(pathToCache, 'simState', 'locIndicesForAllWorkers', ...
        'locIndicesForAllWorkersForAllCellsEff', 'proMon', ...
        'nextIdxEffeCellAnt', '-v7.3');
end

disp('    Done!')

%% Combine Pathloss Results

disp(' ')
disp('    Generating aggregated path loss maps ...')

% Fetch all maps for each cell at this height, and aggregate them using
% pixel-wise min.
[simState.blockageMaps, simState.coverageMaps, ...
    simState.blockageDistMaps, simState.blockageByTerrainDistMaps, ...
    simState.blockageByVegDistMaps, simState.pathLossWithVegMaps] ...
    = deal(cell(1, numOfRxHeightToInspect));

for idxH = 1:numOfRxHeightToInspect
    % Load the results for the first effective tower to prepare for the min
    % comparison. Note that for simplicity we carry out the pixel-wise min
    % with the location coordinates, but they will be discarded eventually.
    simState.blockageMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'Blockage', simConfigs, 'utm');
    simState.coverageMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'Coverage', simConfigs, 'utm');
    simState.blockageDistMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'BlockageDist', simConfigs, 'utm');
    simState.blockageByTerrainDistMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'BlockageByTerrainDist', simConfigs, 'utm');
    simState.blockageByVegDistMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'BlockageByVegDist', simConfigs, 'utm');
    simState.pathLossWithVegMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'PathLossWithVeg', simConfigs, 'utm');

    for idxEffeCell = 2:numOfEffeCellAnts
        simState.blockageMaps{idxH} = min(simState.blockageMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Blockage', simConfigs, 'utm'));
        simState.coverageMaps{idxH} = min(simState.coverageMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Coverage', simConfigs, 'utm'));
        simState.blockageDistMaps{idxH} = min( ...
            simState.blockageDistMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'BlockageDist', simConfigs, 'utm'));
        simState.blockageByTerrainDistMaps{idxH} = min( ...
            simState.blockageByTerrainDistMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'BlockageByTerrainDist', simConfigs, 'utm'));
        simState.blockageByVegDistMaps{idxH} = min( ...
            simState.blockageByVegDistMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'BlockageByVegDist', simConfigs, 'utm'));
        simState.pathLossWithVegMaps{idxH} = min( ...
            simState.pathLossWithVegMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'PathLossWithVeg', simConfigs, 'utm'));
    end

    % Discard the location information in the aggregated path loss maps.
    simState.blockageMaps{idxH} = simState.blockageMaps{idxH}(:,3);
    simState.coverageMaps{idxH} = simState.coverageMaps{idxH}(:,3);
    simState.blockageDistMaps{idxH} = simState.blockageDistMaps{idxH}(:,3);
    simState.blockageByTerrainDistMaps{idxH} ...
        = simState.blockageByTerrainDistMaps{idxH}(:,3);
    simState.blockageByVegDistMaps{idxH} ...
        = simState.blockageByVegDistMaps{idxH}(:,3);
    simState.pathLossWithVegMaps{idxH} ...
        = simState.pathLossWithVegMaps{idxH}(:,3);
end

disp('    Done!')

end
% EOF