function [ simState, simConfigs ] ...
    = analyzeCoverageViaChannelSimulation(pathToSaveResults, ...
    lidarFileAbsDirs, lidarFileXYCoveragePolyshapes, ...
    cellAntsXYH, simConfigs)
%ANALYZECOVERAGEVIACHANNELSIMULATION Generate coverage and blockage maps
%via channel simulations.
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
%     data available as .img files, and their polygon boundries in the UTM
%     (x, y) system.
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
%         generating the coverage maps; this is default to the range
%         covered by the input LiDAR data set when it is empty.
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
%         blockage maps and the coverage maps, respectively, for all the RX
%         height that have been inspected.
%       - blockageMapsForEachCell, coverageMapsForEachCell
%         Two cells for the blockage and coverage maps (for all cell
%         antennas affecting the simulation for all drone heights),
%         respectively. For example,
%               simState.blockageMapsForEachCell{idxCell}{idxDroneHeight}
%         is the path loss vector (for locations in simState.mapGridXYPts)
%         for the idxCell-th cellular tower that has effect on the
%         simulation, at the idxDroneHeight-th drone height that needs to
%         be inspected.
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
% Yaguang Zhang, Purdue, 09/10/2019

%% Before Processing the Data

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% Overall boundries for the area covered by the LiDAR data set in UTM.
lidarFilesXYCoveragePolyshape ...
    = mergePolygonsForAreaOfInterest(lidarFileXYCoveragePolyshapes, 1);
% Centroids for the LiDAR files in UTM.
lidarFileXYCentroids ...
    = extractCentroidsFrom2DPolyCell(lidarFileXYCoveragePolyshapes);

% The area of interest is default to the area covered by the input LiDAR
% data.
if isempty(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST)
    [lidarFilesXYCoveragePolyshapeXs, lidarFilesXYCoveragePolyshapeYs] ...
        = boundary(lidarFilesXYCoveragePolyshape);
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST ...
        = [lidarFilesXYCoveragePolyshapeXs, ...
        lidarFilesXYCoveragePolyshapeYs];
end
assert(length(regions( ...
    polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST)))==1, ...
    'The area of interest should have only one region!');

% Default maximum radius around the area of interest to find cellular
% towers that are effective in the simulation.
if ~(simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M>0)
    if ~isfield(simConfigs, 'MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY')
        simConfigs.MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY ...
            = 'LowestAntennas';
    end

    switch lower(simConfigs.MAX_CELL_COVERAGE_RADIUS_GEN_STRATEGY)
        case 'lowestantennas'
            fctMaxCellCoverageRadiusGenStrategy = @(x) min(x);
        case 'highestantennas'
            fctMaxCellCoverageRadiusGenStrategy = @(x) max(x);
        otherwise
            error('');
    end

    % The maximum range for a tower. For example, it can be the distance
    % that one can see at the top of the highest/lowest cell tower to the
    % highest/lowest RX above the horizon (i.e. without blockage of the
    % earth):
    %       D_BL_IN_KM ~= 3.57(sqrt(h_TX_IN_M)+sqrt(h_RX_IN_M)).
    simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M ...
        = simConfigs.getCellCoverageRadiusInM( ...
        fctMaxCellCoverageRadiusGenStrategy(cellAntsXYH(:,3)), ...
        fctMaxCellCoverageRadiusGenStrategy( ...
        simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M));
end

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

%% Create the Map Grid

mapMinX = min(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
mapMaxX = max(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1));
mapMinY = min(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
mapMaxY = max(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

mapWidthInM = mapMaxX-mapMinX;
mapHeightInM = mapMaxY-mapMinY;

gridResolution = max([mapWidthInM, mapHeightInM]) ...
    ./simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE;

mapXLabels = constructAxisGrid( ...
    mean([mapMaxX, mapMinX]), ...
    floor((mapMaxX-mapMinX)./gridResolution), gridResolution);
mapYLabels = constructAxisGrid( ...
    mean([mapMaxY, mapMinY]), ...
    floor((mapMaxY-mapMinY)./gridResolution), gridResolution);
[mapXs,mapYs] = meshgrid(mapXLabels,mapYLabels);

% Discard map grid points out of the area of interest.
boolsMapGridPtsToKeep = InPolygon(mapXs(:), mapYs(:), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

simState.mapGridXYPts = [mapXs(boolsMapGridPtsToKeep), ...
    mapYs(boolsMapGridPtsToKeep)];

% Convert UTM (x, y) to (lat, lon).
[mapGridLats, mapGridLons] = simConfigs.utm2deg_speZone( ...
    simState.mapGridXYPts(:,1), simState.mapGridXYPts(:,2));
simState.mapGridLatLonPts = [mapGridLats, mapGridLons];

% Plot.
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
    simState.mapGridXYPts(:,2), '.', 'MarkerSize', 2.5, 'Color', darkBlue);
adjustFigSizeByContent(hFigAreaOfInterest, [], 'height', 0.9);
axis equal; view(2); %grid on; grid minor;
hLeg = legend([hAreaOfInterest, hGridPts], ...
    'Area of interest', 'UAV location grid');
xlabel('UTM x (m)'); ylabel('UTM y (m)');
% Adjust legend the exponent label for y axis.
switch lower(simConfigs.CURRENT_SIMULATION_TAG)
    case 'tipp'
        annotation(hFigAreaOfInterest, 'textbox',...
            [0.1409 0.8502 0.1785 0.0945],...
            'String', ['\times10^', num2str(hCurAxis.YAxis.Exponent)],...
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
            'String', ['\times10^', num2str(hCurAxis.YAxis.Exponent)],...
            'FontWeight', hCurAxis.YAxis.FontWeight,...
            'FontSize', hCurAxis.YAxis.FontSize,...
            'EdgeColor', 'none');
        yticks(yticks);
        yticklabels(yticklabels);
end

% tightfig(hFigAreaOfInterest);
curDirToSave = fullfile(pathToSaveResults, 'Overview_RxLocGrid');
saveEpsFigForPaper(hFigAreaOfInterest, curDirToSave);
% saveas(hFigAreaOfInterest,  [curDirToSave, '.eps'], 'epsc');
%  saveas(hFigAreaOfInterest,  [curDirToSave, '.png']);

%% Load Antenna Info

% Keep only the cell towers which can cover some part of the area of
% interest.
utmXYBoundaryToKeepCellTowers ...
    = extendPoly(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ...
    simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);
if length(utmXYBoundaryToKeepCellTowers)==1
    % Note that this polygon vector should already be closed.
    utmXYBoundaryToKeepCellTowers = utmXYBoundaryToKeepCellTowers{1};
else
    warning('The extended area of interest have more than one region!');
    utmXYBoundaryToKeepCellTowersAllPts  ...
        = vertcat(utmXYBoundaryToKeepCellTowers{:});
    % Combine all polygons when necessary.
    utmXYBoundaryToKeepCellTowersIndices = boundary( ...
        utmXYBoundaryToKeepCellTowersAllPts(:,1), ...
        utmXYBoundaryToKeepCellTowersAllPts(:,2));
    utmXYBoundaryToKeepCellTowers ...
        = utmXYBoundaryToKeepCellTowersAllPts( ...
        utmXYBoundaryToKeepCellTowersIndices, :);
end

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
legend([hAreaOfInterest, hExtendedArea, hEffeCells, hIneffeCells], ...
    'Area of interest', 'Extended area', 'Cell towers to consider', ...
    'Ineffective cell towers');

curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider.png');
saveas(hFigCellOverview, curDirToSave);

% Overview plot in GPS (lon, lat).
[effeCellAntsLats, effeCellAntsLons] ...
    = simConfigs.utm2deg_speZone(effeCellAntsXYH(:,1), ...
    effeCellAntsXYH(:,2));
[gpsLatsBoundaryToKeepCellTowers, gpsLonsBoundaryToKeepCellTowers] ...
    = simConfigs.utm2deg_speZone(utmXYBoundaryToKeepCellTowers(:,1), ...
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
end
curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider_RoadMap');
saveEpsFigForPaper(hFigCellOverview, curDirToSave);

plot_google_map('MapType', 'satellite');
curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider_SatelliteMap');
saveEpsFigForPaper(hFigCellOverview, curDirToSave);

%% Simulation

[~, hostname] = system('hostname');
hostname = strtrim(hostname);
pathToHistoryWorkspace = fullfile(pathToSaveResults, ...
    ['CovAnalysisHistory_Task_', ...
    simConfigs.CURRENT_SIMULATION_TAG, ...
    '_LidarSet_', simConfigs.LIDAR_DATA_SET_TO_USE, ...
    '_TerProRes_', ...
    num2str(simConfigs ...
    .MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M), ...
    '_NumPix_', ...
    num2str(simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE), ...
    '_LidProRes_', ...
    num2str(simConfigs ...
    .MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M), ...
    '_HostName_', hostname, '.mat']);
if exist(pathToHistoryWorkspace, 'file')
    disp(' ')
    disp('    Loading history path loss results ...')

    load(pathToHistoryWorkspace); %#ok<LOAD>
else
    disp(' ')
    disp('    Computing path losses ...')
    disp(' ')
    disp('        Closing figures to save RAM ...')
    close all;
    disp('        Done!')

    [numOfEffeCellAnts, ~] = size(effeCellAntsXYH);

    lidarMatFileAbsDirs = cellfun(@(d) regexprep(d, '\.img$', '.mat'), ...
        lidarFileAbsDirs, 'UniformOutput', false);
    [simState.blockageMapsForEachCell, simState.coverageMapsForEachCell, ...
        simState.TimeUsedInSForEachPixel] ...
        = deal(cell(numOfEffeCellAnts, 1));

    % Suppress this warning in the cluster to get clearer feedbacks from
    % the program.
    parfevalOnAll(gcp(), ...
        @warning, 0, 'off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

    % Progress monitor (proMon): for reporting progress and estimating the
    % remaining simulation time.
    proMonPixCnt = 0;
    proMonFloatFomatter = '%.2f';
    proMonNumOfPixToProcess = 0;

    % For each effective cellular tower, find grid pixels that are within
    % the coverage range and pre-assign them to workers. Also, initialize
    % result cells and count the total number of pixels to be processed for
    % progress monitoring.
    locIndicesForAllWorkersForAllCellsEff = cell(numOfEffeCellAnts,1);
    % For initializing results.
    numOfRxHeightToInspect ...
        = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);
    [numPixelsPerMap, ~] = size(simState.mapGridXYPts);
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
            simState.TimeUsedInSForEachPixel{idxEffeCellAnt} ...
                {idxDroneHeightInM} = zeros(1, numPixelsPerMap);
        end

        % Update the total number of pixels to be processed.
        proMonNumOfPixToProcess = proMonNumOfPixToProcess ...
            + length(curIndicesRxLocsToConsider);
    end

    % Path to cache computation results.
    absPathToCacheProfilesDir ...
        = fullfile(pathToSaveResults, 'CachedTerrainProfiles');
    if exist(absPathToCacheProfilesDir, 'dir')~=7
        mkdir(absPathToCacheProfilesDir);
    end

    % Get the number of workers available for estimating progress.
    localCluster = parcluster('local');
    numOfWorkersInLocalCluster = localCluster.NumWorkers;
    for idxEffeCellAnt = 1:numOfEffeCellAnts
        % Cellular location.
        curCellXYH = effeCellAntsXYH(idxEffeCellAnt, :);

        % Load the task assignment results.
        locIndicesForAllWorkers ...
            = locIndicesForAllWorkersForAllCellsEff{idxEffeCellAnt};
        numOfWorkers = length(locIndicesForAllWorkers);

        % Pre-allocate space. To make sure parfor works, we will store all
        % results from one worker into one matrix with each row being:
        % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
        % idxDroneHeightInM].
        resultsFromWorkersCell = cell(numOfWorkers,1);
        for idxWorker = 1:numOfWorkers
            resultsFromWorkersCell{idxWorker} ...
                = nan(length(locIndicesForAllWorkers{idxWorker}) ...
                .*numOfRxHeightToInspect, 5);
        end

        % Utilize parfor only when a majority of the workers will get
        % something to process.
        if (sum(cellfun(@(c) isempty(c), locIndicesForAllWorkers))...
                /numOfWorkers)>=0.5
            parforArg = 0;
        else
            parforArg = numOfWorkers;
        end

        % To make sure the overhead time is recorded by the first pixel
        % processed by the work.
        curOverheadTimeInSecStart = tic;
        curExecTimeInSecStart = curOverheadTimeInSecStart;

        % For recording and estimating processing time in parfor.
        pathToSaveOverheadTimeMats = fullfile(pathToSaveResults, ...
            'ProcessingTimeCacheRecords');
        if exist(pathToSaveOverheadTimeMats, 'dir')
            % Specify 's' to also attempts to remove all subfolders and
            % files, regardless of their write permissions.
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
                    'curExecTimeInSecStart');
                pathsToOverheadTimeInSecStarts{curTaskId} ...
                    = curPathToOverheadTimeInSecStart;
            end
        else
            pathsToOverheadTimeInSecStarts = cell(1, 1);
            curPathToOverheadTimeInSecStart ...
                = fullfile(pathToSaveOverheadTimeMats, ...
                'Worker_0.mat');
            save(curPathToOverheadTimeInSecStart, ...
                'curExecTimeInSecStart');
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
                = load(pathsToOverheadTimeInSecStarts{curTaskId});  ...
                %#ok<PFBNS>
            curExecTimeInSecStart = curExecTimeInSecStartMatContent ...
                .curExecTimeInSecStart;

            % Load the NTIA eHata library first, if necessary, to avoid the
            % "unable to find ehata" error.
            if ~libisloaded('ehata')
                loadlibrary('ehata');
            end

            % Load our Python module.
            py_addpath(fullfile(pwd, 'lib', 'python'));

            curWorkerPixCnt = 0;
            curWorkerNumPixs = length(locIndicesForAllWorkers{idxWorker});
            curWorkerNumPixsToReportProgress = ceil(curWorkerNumPixs ...
                .*simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT ...
                ); %#ok<PFBNS>
            curDroneLocIndices = locIndicesForAllWorkers{idxWorker};
            % Make sure curDroneLocIndices is a row vector.
            curDroneLocIndices = curDroneLocIndices(:)';
            for idxDroneLoc = curDroneLocIndices

                % Report progress when necessary.
                if mod(curWorkerPixCnt, ...
                        curWorkerNumPixsToReportProgress)==0
                    disp(['        Worker #', num2str(idxWorker), ...
                        '/', num2str(numOfWorkers), ' (', ...
                        num2str(curWorkerPixCnt/curWorkerNumPixs ...
                        /numOfRxHeightToInspect*100, ...
                        '%.2f'), '%) ...']);
                end

                % Drone location.
                curDroneXY ...
                    = simState.mapGridXYPts(idxDroneLoc, :); %#ok<PFBNS>

                % A unqiue .mat directory for caching the profiles. Note
                % that we are using the unique IDs for the cellular towers,
                % so there is no need to add labels for different max cell
                % coverage radii.
                uniqueMatFileName = ['Task_', ...
                    simConfigs.CURRENT_SIMULATION_TAG, ...
                    '_LidarSet_', simConfigs.LIDAR_DATA_SET_TO_USE, ...
                    '_TerProRes_', ...
                    num2str(simConfigs ...
                    .MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M), ...
                    '_NumPix_', ...
                    num2str(simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE), ...
                    '_LidProRes_', ...
                    num2str(simConfigs ...
                    .MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M), ...
                    '_CellId_', num2str( ...
                    simState.CellAntsEffectiveIds(idxEffeCellAnt)), ...
                    '_DroLocIdx_', num2str(idxDroneLoc)];
                absPathToCacheMatFile ...
                    = fullfile(absPathToCacheProfilesDir, ...
                    [uniqueMatFileName, '.mat']);

                % Generate terrain and LiDAR profiles.
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

                    [curBlockagePL, curCoveragePL] ...
                        = computeBlockageAndCoveragePLs( ...
                        curLidarProfile, curTerrainProfile, ...
                        curCellXYH, [curDroneXY, curDroneH], ...
                        simConfigs);

                    curExecTime = toc(curExecTimeInSecStart);

                    % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
                    % idxDroneHeightInM].
                    curWorkerPixCnt = curWorkerPixCnt+1;
                    resultsFromWorkersCell{idxWorker} ...
                        (curWorkerPixCnt, :) ...
                        = [curBlockagePL, curCoveragePL, curExecTime, ...
                        idxDroneLoc, idxDroneHeightInM];

                    % Reset timer.
                    curExecTimeInSecStart = tic;
                    curPathToSaveOverheadTimeInSecStart ...
                        = pathsToOverheadTimeInSecStarts{curTaskId};
                    parsave(curPathToSaveOverheadTimeInSecStart, ...
                        curExecTimeInSecStart);
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

            simState.blockageMapsForEachCell{idxEffeCellAnt} ...
                {curIdxDroneHeightInM}(curIdxDroneLoc) = curBlockagePL;
            simState.coverageMapsForEachCell{idxEffeCellAnt} ...
                {curIdxDroneHeightInM}(curIdxDroneLoc) = curCoveragePL;
            simState.TimeUsedInSForEachPixel{idxEffeCellAnt} ...
                {curIdxDroneHeightInM}(curIdxDroneLoc) = curPixelExecTime;
        end

        latestEstiPixExecTimeInS = mean(resultsFromWorkersMat(:,3));
        latestEstiTowerExecTimeInS = sum(resultsFromWorkersMat(:,3)) ...
            ./numOfWorkersInLocalCluster;
        % Take into account all the heights inspected.
        proMonPixCnt = proMonPixCnt ...
            + numOfResultsFromWorkers./numOfRxHeightToInspect;

        disp(' ')
        disp(['        [', datestr(now, 'yyyy/mm/dd HH:MM:ss'), ...
            '] Finished cellular tower #', ...
            num2str(idxEffeCellAnt), '/', num2str(numOfEffeCellAnts)]);
        disp(['            Estimated time used for this tower: ', ...
            num2str(latestEstiTowerExecTimeInS, proMonFloatFomatter), ...
            ' seconds']);
        if ~isnan(latestEstiTowerExecTimeInS)
            disp(['                (', ...
                seconds2human(latestEstiTowerExecTimeInS), ')']);
        end
        disp(['            Total progress: ', ...
            num2str(proMonPixCnt./proMonNumOfPixToProcess.*100, ...
            proMonFloatFomatter), ...
            '% (Pixel #', num2str(proMonPixCnt), ...
            '/', num2str(proMonNumOfPixToProcess), ')']);
        remainingTimeInS = latestEstiPixExecTimeInS ...
            .*(proMonNumOfPixToProcess-proMonPixCnt) ...
            .*numOfRxHeightToInspect ...
            ./numOfWorkersInLocalCluster;
        disp(['            Estimated remaining time: ', ...
            num2str(remainingTimeInS, proMonFloatFomatter), ...
            ' seconds']);
        if ~isnan(remainingTimeInS)
            disp(['                (', ...
                seconds2human(remainingTimeInS), ')']);
        end
        disp(' ');
    end

    % Save the whole workspace to a .mat file just in case.
    save(pathToHistoryWorkspace);
end

disp('    Done!')

%% Combine Pathloss Results

disp(' ')
disp('    Generating aggregated path loss maps ...')

% Fetch all maps for each cell at this height, and aggregate them using
% pixel-wise min.
[simState.blockageMaps, simState.coverageMaps] ...
    = deal(cell(1, numOfRxHeightToInspect));
for idxH = 1:numOfRxHeightToInspect
    simState.blockageMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'Blockage', simConfigs, 'utm');
    simState.coverageMaps{idxH} ...
        = fetchPathlossResultsFromSimState(simState, ...
        1, idxH, 'Coverage', simConfigs, 'utm');

    for idxEffeCell = 2:numOfEffeCellAnts
        simState.blockageMaps{idxH} = min(simState.blockageMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Blockage', simConfigs, 'utm'));
        simState.coverageMaps{idxH} = min(simState.coverageMaps{idxH}, ...
            fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Coverage', simConfigs, 'utm'));
    end

    % Discard the location information in the aggregated path loss maps.
    simState.blockageMaps{idxH} = simState.blockageMaps{idxH}(:,3);
    simState.coverageMaps{idxH} = simState.coverageMaps{idxH}(:,3);
end

disp('    Done!')

end

% EOF