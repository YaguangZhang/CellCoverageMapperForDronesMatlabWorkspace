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
%       - MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE
%         The guaranteed minimum number of LiDAR z (or possibly elevation)
%         elements in one terrain profile; this will ensure non-empty
%         terrain profiles.
%
% Output:
%   - simState
%     A structure for the state of the simulation, including fields
%       - mapGridXYPts
%         The UTM (x,y) locations for the map grid points. Note that the
%         shape of the map is determined by
%         simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST.
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
%
% Yaguang Zhang, Purdue, 09/10/2019

%% Before Processing the Data

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% We will save simState for recovery from interrupted simulation processes
% when necessary.
ABS_PATH_TO_SAVE_COMP_PROGRESS = fullfile(pathToSaveResults, ...
    [simConfigs.CURRENT_SIMULATION_TAG, '_SimState.mat']);

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

% Default maximum radius around the cellular tower to consider in the
% simulation.
if ~(simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M>0)
    % Set it to be the distance that one can see at the top of the highest
    % cell tower to the highest RX above the horizon (i.e. without blockage
    % of the earth): D_BL_IN_KM ~= 3.57(sqrt(h_TX_IN_M)+sqrt(h_RX_IN_M)).
    simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M ...
        = 3.57*(sqrt(max(cellAntsXYH(:,3)))...
        +sqrt(max(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M)))*1000;
end

% For plotting.
areaOfInterestColor = [0.9290 0.6940 0.1250];
lightBlue = [0.3010 0.7450 0.9330];
darkBlue = [0 0.4470 0.7410];

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
boolsMapGridPtsToKeep = inpolygon(mapXs(:), mapYs(:), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));

simState.mapGridXYPts = [mapXs(boolsMapGridPtsToKeep), ...
    mapYs(boolsMapGridPtsToKeep)];

% Plot.
hFigAreaOfInterest = figure; hold on;
hAreaOfInterest = plot( ...
    polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
    'FaceColor', areaOfInterestColor);
hGridPts = plot(simState.mapGridXYPts(:,1), ...
    simState.mapGridXYPts(:,2), '.b', 'MarkerSize', 3, 'Color', darkBlue);
axis equal; view(2); grid on; grid minor;
legend([hAreaOfInterest, hGridPts], ...
    'Area of interest', 'RX location grid points', ...
    'Location', 'SouthWest');
transparentizeCurLegends;

curDirToSave = fullfile(pathToSaveResults, 'Overview_RxLocGrid.png');
saveas(hFigAreaOfInterest, curDirToSave);

%% Load Antenna Info

% Keep only the cell towers which can cover some part of the area of
% interest.
utmXYBoundaryToKeepCellTowers ...
    = extendPoly(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST, ...
    simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M);
assert(length(utmXYBoundaryToKeepCellTowers)==1, ...
    'The area of interest should have only one region!');
utmXYBoundaryToKeepCellTowers = utmXYBoundaryToKeepCellTowers{1};

boolsCellAntsToKeep ...
    = inpolygon(cellAntsXYH(:,1), cellAntsXYH(:,2), ...
    utmXYBoundaryToKeepCellTowers(:,1), ...
    utmXYBoundaryToKeepCellTowers(:,2));
% Effective cellular towers.
effeCellAntsXYH = cellAntsXYH(boolsCellAntsToKeep, :);

% Plot.
hFigCellOverview = figure; hold on;
hEffeCells = plot(effeCellAntsXYH(:,1), ...
    effeCellAntsXYH(:,2), '.', 'Color', darkBlue);
plot(polyshape(utmXYBoundaryToKeepCellTowers));
hAreaOfInterest = plot( ...
    polyshape(simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST), ...
    'FaceColor', areaOfInterestColor);
axis equal; view(2); grid on; grid minor;
axisToSet = axis;
hAllCells = plot(cellAntsXYH(:,1), ...
    cellAntsXYH(:,2), '.', 'Color', 'yellow');
uistack(hAllCells, 'bottom'); uistack(hEffeCells, 'top');
axis(axisToSet);
legend([hAreaOfInterest, hEffeCells], ...
    'Area of interest', 'Cellular towers to consider', ...
    'Location', 'SouthWest');
transparentizeCurLegends;

curDirToSave = fullfile(pathToSaveResults, ...
    'Overview_CellularTowersToConsider.png');
saveas(hFigCellOverview, curDirToSave);

disp('    Done!')

%% Simulation

[numOfEffeCellAnts, ~] = size(effeCellAntsXYH);
[numOfDroneLocs, ~] = size(simState.mapGridXYPts);
numOfRxHeightToInspect = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);

lidarMatFileAbsDirs = cellfun(@(d) regexprep(d, '\.img$', '.mat'), ...
    lidarFileAbsDirs, 'UniformOutput', false);

[simState.blockageMapsForEachCell, simState.coverageMapsForEachCell, ...
    simState.TimeUsedInSForEachPixel] ...
    = deal(cell(numOfEffeCellAnts, 1));

% We will use multiple works to churn through the drone locations. To avoid
% repeative environment set up and data transfer, here we pre-assign the
% pixels to be processed by each worker. Pre-assign pixels to workers to
% avoid unnecessary data copying.
locIndicesForAllWorkers = preassignTaskIndicesToWorkers(numOfDroneLocs);
numOfWorkers = length(locIndicesForAllWorkers);

% Suppress this warning in the cluster to get clearer feedbacks from the
% program.
parfevalOnAll(gcp(), ...
    @warning, 0, 'off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

% Progress monitor (proMon): for reporting progress and estimating the
% remaining simulation time.
proMonPixCnt = 0;
proMonNumOfPixToProcess ...
    = numOfEffeCellAnts.*numOfDroneLocs.*numOfRxHeightToInspect;
proMonFloatFomatter = '%.2f';

disp(' ')
disp('    Computing path losses ...')

% Get the number of workers available for estimating progress.
localCluster = parcluster('local');
numOfWorkersInLocalCluster = localCluster.NumWorkers;
for idxEffeCellAnt = 1:numOfEffeCellAnts
    disp(' ')
    disp('        Closing figures to save RAM ...')
    close all;
    disp('        Done!')
    
    % Cellular location.
    curCellXYH = effeCellAntsXYH(idxEffeCellAnt, :);
    
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
    
    parfor idxWorker = 1:numOfWorkers
        % Load the NTIA eHata library first, if necessary, to avoid the
        % "unable to find ehata" error.
        if ~libisloaded('ehata')
            loadlibrary('ehata');
        end
        
        % Load our Python module.
        py_addpath(fullfile(pwd, 'lib', 'python'));
        
        % For recording and estimating processing time.
        curExecTimeInSecStart = tic;
        
        curWorkerPixCnt = 0;
        curWorkerNumPixs = length(locIndicesForAllWorkers{idxWorker});
        curWorkerNumPixsToReportProgress = ceil(curWorkerNumPixs ...
            .*simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT); %#ok<PFBNS>
        for idxDroneLoc = locIndicesForAllWorkers{idxWorker}
            
            % Report progress when necessary.
            if mod(curWorkerPixCnt, curWorkerNumPixsToReportProgress)==0
                disp(['        Worker #', num2str(idxWorker), ...
                    '/', num2str(numOfWorkers), ' (', ...
                    num2str(curWorkerPixCnt/curWorkerNumPixs ...
                    /numOfRxHeightToInspect*100, ...
                    '%.2f'), '%) ...']);
            end
            
            % Drone location.
            curDroneXY = simState.mapGridXYPts(idxDroneLoc, :); %#ok<PFBNS>
            
            % A unqiue .mat directory for caching the profiles.
            absPathToCacheProfilesDir ...
                = fullfile(pathToSaveResults, 'CachedTerrainProfiles');
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
                '_CelLocIdx_', num2str(idxEffeCellAnt), ...
                '_DroLocIdx_', num2str(idxDroneLoc), ...
                '.mat'];
            if ~exist(absPathToCacheProfilesDir, 'dir')
                mkdir(absPathToCacheProfilesDir);
            end
            absPathToCacheMatFile ...
                = fullfile(absPathToCacheProfilesDir, uniqueMatFileName);
            
            % Generate terrain and LiDAR profiles.
            [curTerrainProfile, curLidarProfile] ...
                = fetchTerrainAndLidarProfiles(absPathToCacheMatFile, ...
                curCellXYH(:,1:2), ...
                curDroneXY, simConfigs, ...
                lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
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
                
                curWorkerPixCnt = curWorkerPixCnt+1;
                
                % [blockagePl, coveragePl, pixelExecTime, idxDroneLoc,
                % idxDroneHeightInM].
                resultsFromWorkersCell{idxWorker}(curWorkerPixCnt, :) ...
                    = [curBlockagePL, curCoveragePL, curExecTime, ...
                    idxDroneLoc, idxDroneHeightInM];
                
                % Reset timer.
                curExecTimeInSecStart = tic;
            end
        end
        
        % Final progress report which is expected to be 100% done.
        disp(['        Worker #', num2str(idxWorker), ...
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
    latestEstiTowerExecTimeInS = sum(resultsFromWorkersMat(:,3));
    proMonPixCnt = proMonPixCnt+numOfResultsFromWorkers;
    
    disp(' ')
    disp(['        Finished cellular tower #', ...
        num2str(idxEffeCellAnt), '/', num2str(numOfEffeCellAnts)]);
    disp(['            Times used for this tower: ', ...
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
        ./numOfWorkersInLocalCluster;
    disp(['            Estimated remaining time: ', ...
        num2str(remainingTimeInS, proMonFloatFomatter), ...
        ' seconds']);
    if ~isnan(remainingTimeInS)
        disp(['                (', seconds2human(remainingTimeInS), ')']);
    end
end
disp('    Done!')

%% Plot Maps for Each Cellular Tower
% We will plot the path losses to generate both the blockage maps and the
% coverage maps for each cellular tower and each inspected height.

disp(' ')
disp('    Plotting maps for each cellular tower ...')

parfor idxWorker = 1:numOfWorkers
    curFigFileName = [ ...
        'CellTowerPathLoss_RxHeight_', ...
        num2str(curRxAntH), ...
        '_TxCell_', num2str(idxCellAntenna), ...
        '_eHataLib_', LIBRARY_TO_USE, ...
        '_Terrain_', REGION];
    
    if FLAG_GEN_FIGS
        % Generate a figure on Google map to show the path loss map.
        [curXs, curYs] ...
            = meshgrid(pathLossMapXLabels{idxCellAntenna}, ...
            pathLossMapYLabels{idxCellAntenna});
        curBoolsToShow ...
            = pathLossMaps{idxCellAntenna}(:) ...
            >=EXPECTED_PL_RANGE_IN_DB(1) ...
            & pathLossMaps{idxCellAntenna}(:) ...
            <=EXPECTED_PL_RANGE_IN_DB(2);
        [curLatsToShow, curLonsToShow] ...
            = utm2deg(curXs(curBoolsToShow), ...
            curYs(curBoolsToShow), ...
            repmat(UTM_ZONE, sum(curBoolsToShow), 1));
        
        hCurPLMap = figure;
        plot3k([curLonsToShow curLatsToShow ...
            pathLossMaps{idxCellAntenna}(curBoolsToShow)], ...
            'Labels', {'', ...
            'Longitude (degrees)', 'Latitude (degrees)', ...
            '', 'Path Loss (dB)'}, ...
            'ColorRange', EXPECTED_PL_RANGE_IN_DB);
        grid on; view(2); axis tight;
        plotGoogleMapAfterPlot3k(gcf, 'satellite');
        
        pathToSaveFig = fullfile(pathToSaveResults, ...
            [curFigFileName, '.png']);
        saveas(hCurPLMap, pathToSaveFig);
        
        hCurTransPLMap = figure;
        curAlpha = 0.95;
        curDotSize = 15;
        scatter(curLonsToShow, curLatsToShow, ...
            curDotSize.*ones(length(curLatsToShow),1), ...
            pathLossMaps{idxCellAntenna}(curBoolsToShow), ...
            'Marker', '.', ...
            'MarkerEdgeAlpha', curAlpha);
        curCB = colorbar; curCB.Label.String = 'Path Loss (dB)';
        xlabel('Longitude (degrees)');
        ylabel('Latitude (degrees)')
        grid on; view(2); axis tight;
        plot_google_map('MapType', 'satellite');
        
        pathToSaveFig = fullfile(pathToSaveResults, ...
            [curFigFileName, '_Transparent.png']);
        saveas(hCurTransPLMap, pathToSaveFig);
    end
end

disp('    Done!')

%% Combine Pathloss Results


%% Clear Things Up


end
% EOF