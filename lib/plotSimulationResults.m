function plotSimulationResults( ...
    pathToSaveResults, simState, simConfigs)
%PLOTSIMULATIONRESULTS Generate plots for the simulation results from
%analyzeCellularCoverage.m.
%
% Inputs:
%   - pathToSaveResults
%     The full path to the directory for saving plots.
%   - simState, simConfigs
%     The outputs from simulateChannelForExtendedTipp.m representing the
%     simulation results and simulaiton configurations, respectively.
%
% This script is based on:
%   5_SimulationForExtendedTipp/plotCoverageSimulationResults.m.
%
% Yaguang Zhang, Purdue, 05/07/2021

% Parameters.
numOfFigsToGenForSingleCellInspection = 50;
[numOfEffeCellAnts, ~] = size(simState.CellAntsXyhEffective);
numOfRxHeightToInspect = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);
figSizeRatioForIn = 1.5;

%% Plot Some Maps for Each Cellular Tower
% We will plot the path losses to generate both the blockage maps and the
% coverage maps for each cellular tower and each inspected height. We will
% only plot a randomly chosen subset of all the figures to save Google Maps
% quota.

disp(' ')
disp('    Plotting maps for each cellular tower ...')

% Path to save maps for single celluar towers.
dirToSaveMapsForSingleCellTowers = fullfile(pathToSaveResults, ...
    'MapsForSingleCellularTowers');
% Create directories if necessary.
if exist(dirToSaveMapsForSingleCellTowers, 'dir')~=7
    mkdir(dirToSaveMapsForSingleCellTowers);
end

totalNumOfFigs = numOfEffeCellAnts.*numOfRxHeightToInspect;
% Randomly choose a subset of all the figures to plot if necessary.
if totalNumOfFigs>numOfFigsToGenForSingleCellInspection %#ok<BDSCI>
    s = RandStream('mlfg6331_64');
    indicesOfFigTasks = randsample(s,totalNumOfFigs, ...
        numOfFigsToGenForSingleCellInspection);
    indicesOfFigTasks = sort(indicesOfFigTasks);
else
    indicesOfFigTasks = 1:totalNumOfFigs;
end

[effeCellIndices, rxHeightIndices] ...
    = meshgrid(1:numOfEffeCellAnts, 1:numOfRxHeightToInspect);
matEffeCellIdRxHId = [effeCellIndices(:), rxHeightIndices(:)];

matEffeCellIdRxHIdToDraw = matEffeCellIdRxHId(indicesOfFigTasks,:);
[totalNumOfFigsToDraw,~] = size(matEffeCellIdRxHIdToDraw);

% Pre-assign the plot tasks to workers.
plotIndicesForAllWorkers ...
    = preassignTaskIndicesToWorkers(totalNumOfFigsToDraw, 1);
% plotIndicesForAllWorkers ...
%     = preassignTaskIndicesToWorkers(totalNumOfFigsToDraw);
numOfWorkers = length(plotIndicesForAllWorkers);

matEffeCellIdRxHIdForAllWorders = cell(numOfWorkers,1);
for idxWorker = 1:numOfWorkers
    matEffeCellIdRxHIdForAllWorders{idxWorker} ...
        = matEffeCellIdRxHIdToDraw(plotIndicesForAllWorkers{idxWorker}, :);
end

% Blockage and coverage maps.
%   Bug: plot_google_maps not working in parfor.
for idxWorker = 1:numOfWorkers
    curEffeCellIds = matEffeCellIdRxHIdForAllWorders{idxWorker}(:, 1);
    curRxHIds = matEffeCellIdRxHIdForAllWorders{idxWorker}(:, 2);

    [curWorkerNumFigSets, ~] ...
        = size(matEffeCellIdRxHIdForAllWorders{idxWorker});
    curWorkerNumFigSetsToReportProgress = ceil(curWorkerNumFigSets ...
        .*simConfigs.WORKER_MIN_PROGRESS_RATIO_TO_REPORT);

    for idxFigSet = 1:curWorkerNumFigSets
        % Report progress when necessary.
        if mod(idxFigSet-1, curWorkerNumFigSetsToReportProgress)==0
            disp(['        Worker #', num2str(idxWorker), '/', ...
                num2str(numOfWorkers), ': ', ...
                num2str(idxFigSet),'/', ...
                num2str(curWorkerNumFigSets), ' (', ...
                num2str((idxFigSet-1)/curWorkerNumFigSets*100, ...
                '%.2f'), '%) ...']);
        end

        idxEffeCell = curEffeCellIds(idxFigSet);
        idxH = curRxHIds(idxFigSet);

        % Fetch simulation results for the current cellular tower at the
        % current inspected RX height.
        locType = 'GPS';
        curFlagGenFigsSilently = true;

        % Blockage map.
        [matRxLonLatWithPathloss, rxAntH, cellAntLonLat, ~] ...
            = fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Blockage', simConfigs, locType);

        % Use plot3k for the blockage maps to avoid fitting data for spots
        % with NaN path loss values.
        pathToSaveBlockFig = fullfile(dirToSaveMapsForSingleCellTowers, ...
            ['CellTowerPathLoss_Type_', 'Blockage', ...
            '_RxHeight_', num2str(rxAntH), ...
            '_EffeCellId_', num2str(idxEffeCell), '.png']);
        if ~exist(pathToSaveBlockFig, 'file')
            [ hCurPLMap ] ...
                = plotPathLossMap(matRxLonLatWithPathloss, ...
                cellAntLonLat, simConfigs, ~curFlagGenFigsSilently, ...
                false, 'plot3k');

            saveas(hCurPLMap,  pathToSaveBlockFig);
            close(hCurPLMap);
        end

        % Blockage distance map.
        [matRxLonLatWithDist, rxAntH, ~, ~] ...
            = fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'BlockageDist', simConfigs, locType);

        pathToSaveBlockDistFig = fullfile( ...
            dirToSaveMapsForSingleCellTowers, ...
            ['CellTowerPathLoss_Type_', 'BlockageDist', ...
            '_RxHeight_', num2str(rxAntH), ...
            '_EffeCellId_', num2str(idxEffeCell), '.png']);
        if ~exist(pathToSaveBlockDistFig, 'file')
            [ hCurDistMap ] = plotPathLossMap(matRxLonLatWithDist, ...
                cellAntLonLat, simConfigs, ~curFlagGenFigsSilently, ...
                false, 'griddatasurf');

            saveas(hCurDistMap,  pathToSaveBlockDistFig);
            close(hCurDistMap);
        end

        % Coverage map.
        [matRxLonLatWithPathloss, rxAntH, cellAntLonLat, ~] ...
            = fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'Coverage', simConfigs, locType);
        pathToSaveCovFig = fullfile(dirToSaveMapsForSingleCellTowers, ...
            ['CellTowerPathLoss_Type_', 'Coverage', ...
            '_RxHeight_', num2str(rxAntH), ...
            '_EffeCellId_', num2str(idxEffeCell), '.png']);
        if ~exist(pathToSaveCovFig, 'file')
            try
                [ hCurPLMap ] ...
                    = plotPathLossMap(matRxLonLatWithPathloss, ...
                    cellAntLonLat, simConfigs, ~curFlagGenFigsSilently);
            catch
                if isvalid(hCurPLMap)
                    close(hCurPLMap);
                end
                [ hCurPLMap ] ...
                    = plotPathLossMap(matRxLonLatWithPathloss, ...
                    cellAntLonLat, simConfigs, ~curFlagGenFigsSilently, ...
                    false, 'plot3k');
            end

            saveas(hCurPLMap,  pathToSaveCovFig);
            close(hCurPLMap);
        end

        % Coverage path loss map considering vegetation.
        [matRxLonLatWithPathloss, rxAntH, cellAntLonLat, ~] ...
            = fetchPathlossResultsFromSimState(simState, ...
            idxEffeCell, idxH, 'PathLossWithVeg', simConfigs, locType);
        pathToSaveCovFig = fullfile(dirToSaveMapsForSingleCellTowers, ...
            ['CellTowerPathLoss_Type_', 'PathLossWithVeg', ...
            '_RxHeight_', num2str(rxAntH), ...
            '_EffeCellId_', num2str(idxEffeCell), '.png']);
        if ~exist(pathToSaveCovFig, 'file')
            try
                [ hCurPLMap ] ...
                    = plotPathLossMap(matRxLonLatWithPathloss, ...
                    cellAntLonLat, simConfigs, ~curFlagGenFigsSilently);
            catch
                if isvalid(hCurPLMap)
                    close(hCurPLMap);
                end
                [ hCurPLMap ] ...
                    = plotPathLossMap(matRxLonLatWithPathloss, ...
                    cellAntLonLat, simConfigs, ~curFlagGenFigsSilently, ...
                    false, 'plot3k');
            end

            saveas(hCurPLMap,  pathToSaveCovFig);
            close(hCurPLMap);
        end
    end

    disp(['        Worker #', num2str(idxWorker), ': ', ...
        num2str(idxFigSet),'/', ...
        num2str(curWorkerNumFigSets), ' (', ...
        num2str(idxFigSet/curWorkerNumFigSets*100, ...
        '%.2f'), '%) ...']);
end

disp('    Done!')

%% Plots for Aggregated Maps

disp(' ')
disp('    Plotting aggregated path loss maps ...')

curFlagGenFigsSilently = true;
defaultCmdToPlotPLMaps = 'surf';

[effeCellAntLats, effeCellAntLons] = simConfigs.utm2deg_speZone( ...
    simState.CellAntsXyhEffective(:, 1), ...
    simState.CellAntsXyhEffective(:, 2));
effeCellAntLonLats = [effeCellAntLons, effeCellAntLats];
mapGridLonLats = simState.mapGridLatLonPts(:, [2,1]);

% Overview for settings.
[numOfPixels, ~] = size(simState.blockageMaps{1});
[ hCurPLMap ] = plotPathLossMap( ...
    [mapGridLonLats, nan(numOfPixels, 1)], ...
    effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, false, ...
    defaultCmdToPlotPLMaps);
legend off;

pathToSaveFig = fullfile(pathToSaveResults, ...
    'Overview_SimulationSetup.png');
saveas(hCurPLMap,  pathToSaveFig);
close(hCurPLMap);

% Smaller maps for publication.
if ~isfield(simState,'flagResizeFigForPublication')
    % By default, we do not need to resize figures.
    evalin('base', 'flagResizeFigForPublication = false;');
else
    evalin('base', ...
        'flagResizeFigForPublication=simState.RESIZE_FIG_FOR_PUBLICATION');
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

% Adjustment for IN.
if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
    customFigSize = customFigSize.*figSizeRatioForIn;
end

% Make sure the color bar for all blockage distance maps (with different
% drone heights) is of the same limit.
maxBlockDistInM = max(vertcat(simState.blockageDistMaps{:}));
for idxH = 1:numOfRxHeightToInspect
    rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);

    %------ 1 ------
    % For blockage path loss maps.
    %---------------
    curFlagZoomIn = true;
    % Use plot3k to avoid fitting data for NaN spots.
    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, simState.blockageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, 'plot3k', customFigSize);
    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
        hLeg = findobj(hCurPLMap, 'Type', 'Legend');
        set(hLeg, 'Location', 'best');
    end

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_Blockage_RxHeight_', num2str(rxAntH), '.png']);
    export_fig(hCurPLMap, pathToSaveFig, '-m1');
    close(hCurPLMap);

    %------ 2 ------
    % For blockage status maps.
    %---------------
    [ hCurBlMap ] ...
        = plotBlockageMap( ...
        [mapGridLonLats, simState.blockageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, customFigSize);
    hLeg = findobj(hCurBlMap, 'Type', 'Legend');
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'extendedtipp'
            set(hLeg, 'Position', [0.1369, 0.7687, 0.4374, 0.1538]);
        case 'tipp'
            set(hLeg, 'Position', [0.5201, 0.1132, 0.3836, 0.1538]);
        case 'shrinkedin'
            set(hLeg, 'Location', 'best');
    end
    % Further tighten the figure.
    xlabel(''); ylabel('');
    tightfig(hCurBlMap);

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['BlockageStatusMap_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurBlMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurBlMap, [pathToSaveFig, '.png']);
    saveas(hCurBlMap, [pathToSaveFig, '.fig']);
    close(hCurBlMap);

    %------ 3 ------
    % For blockage distance maps.
    %---------------
    [ hCurDistMap ] = plotPathLossMap( ...
        [mapGridLonLats, simState.blockageDistMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, 'griddatasurf', customFigSize);

    if strcmpi(simConfigs.CURRENT_SIMULATION_TAG, 'shrinkedin')
        hLeg = findobj(gcf, 'Type', 'Legend');
        set(hLeg, 'Location', 'best');
    end

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['BlockageDist_RxHeight_', num2str(rxAntH)]);
    saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurDistMap, [pathToSaveFig, '.png']);
    saveas(hCurDistMap, [pathToSaveFig, '.fig']);

    %------ 3_extra ------
    % The version with colorbar adjusted to be the same one.
    %---------------------
    caxis([0, maxBlockDistInM]);
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['BlockageDist_SameColorBar_RxHeight_', num2str(rxAntH)]);
    saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurDistMap, [pathToSaveFig, '.png']);
    saveas(hCurDistMap, [pathToSaveFig, '.fig']);

    %------ 3_extra_extra ------
    % The version with colorbar adjusted to be the same one but hidden.
    %---------------------------
    hCb = findobj(gcf, 'Type', 'Colorbar');
    set(hCb, 'Visible', 'off');
    tightfig(hCurDistMap);

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['BlockageDist_SameColorBarHidden_RxHeight_', num2str(rxAntH)]);
    saveas(hCurDistMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurDistMap, [pathToSaveFig, '.png']);
    saveas(hCurDistMap, [pathToSaveFig, '.fig']);

    close(hCurDistMap);

    %------ 4 ------
    % For coverage path loss maps.
    %---------------
    curFlagZoomIn = true;

    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, simState.coverageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, defaultCmdToPlotPLMaps, customFigSize);
    hLeg = findobj(gcf, 'Type', 'Legend');
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'extendedtipp'
            set(hLeg, 'Position', [0.1075, 0.8640, 0.3745, 0.0590]);
        case 'tipp'
            set(hLeg, 'Position', [0.4053, 0.1144, 0.3389, 0.0629]);
        case 'shrinkedin'
            set(hLeg, 'Location', 'best');
    end

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_Coverage_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    %------ 4_extra ------
    % The version with colorbar hidden.
    %---------------------
    hCb = findobj(gcf, 'Type', 'Colorbar');
    set(hCb, 'Visible', 'off');
    tightfig(hCurPLMap);

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_Coverage_ColorBarHidden_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    close(hCurPLMap);

    %------ 5 ------
    % For coverage path loss maps considering propogation loss through
    % vegetation.
    %---------------
    curFlagZoomIn = true;

    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, simState.pathLossWithVegMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, defaultCmdToPlotPLMaps, customFigSize);
    hLeg = findobj(gcf, 'Type', 'Legend');
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'extendedtipp'
            set(hLeg, 'Position', [0.1075, 0.8640, 0.3745, 0.0590]);
        case 'tipp'
            set(hLeg, 'Position', [0.4053, 0.1144, 0.3389, 0.0629]);
        case 'shrinkedin'
            set(hLeg, 'Location', 'best');
    end

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_WithVeg_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    %------ 5_extra ------
    % The version with colorbar hidden.
    %---------------------
    hCb = findobj(gcf, 'Type', 'Colorbar');
    set(hCb, 'Visible', 'off');
    tightfig(hCurPLMap);

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_WithVeg_ColorBarHidden_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    %------ 6 ------
    % For path loss from propogation loss through vegetation only.
    %---------------
    curFlagZoomIn = true;

    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, ...
        simState.pathLossWithVegMaps{idxH} ...
        - simState.coverageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, defaultCmdToPlotPLMaps, customFigSize);
    hLeg = findobj(gcf, 'Type', 'Legend');
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'extendedtipp'
            set(hLeg, 'Position', [0.1075, 0.8640, 0.3745, 0.0590]);
        case 'tipp'
            set(hLeg, 'Position', [0.4053, 0.1144, 0.3389, 0.0629]);
        case 'shrinkedin'
            set(hLeg, 'Location', 'best');
    end

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_ByVeg_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    %------ 6_extra ------
    % The version with colorbar hidden.
    %---------------------
    hCb = findobj(gcf, 'Type', 'Colorbar');
    set(hCb, 'Visible', 'off');
    tightfig(hCurPLMap);

    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_ByVeg_ColorBarHidden_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap, [pathToSaveFig, '.png']);
    saveas(hCurPLMap, [pathToSaveFig, '.fig']);

    close(hCurPLMap);
end

disp('    Done!')

%% Figures for Computation Time Statistics

disp(' ')
disp('    Plotting processing time statistics ...')

% Processing time for all pixels.
curFlagGenFigsSilently = true;
curPlotType = 'Pixels';

[curProcTimeFig, curProcTimeHistFig] ...
    = plotProcessingTime(simState, curPlotType, ...
    ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '.png']);
saveas(curProcTimeFig,  pathToSaveFig);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '_Hist.png']);
saveas(curProcTimeHistFig,  pathToSaveFig);
if curFlagGenFigsSilently
    close(curProcTimeFig);
    close(curProcTimeHistFig);
end

% Processing time for all cellular towers.
curPlotType = 'MapSets';

[curProcTimeFig, curProcTimeHistFig] ...
    = plotProcessingTime(simState, curPlotType, ...
    ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '.png']);
saveas(curProcTimeFig,  pathToSaveFig);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '_Hist.png']);
saveas(curProcTimeHistFig,  pathToSaveFig);
if curFlagGenFigsSilently
    close(curProcTimeFig);
    close(curProcTimeHistFig);
end

% Processing time for all cellular towers.
curPlotType = 'CellTowers';

[curProcTimeFig, curProcTimeHistFig] ...
    = plotProcessingTime(simState, curPlotType, ...
    ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '.png']);
saveas(curProcTimeFig,  pathToSaveFig);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['ProcessingTime_', curPlotType, '_Hist.png']);
saveas(curProcTimeHistFig,  pathToSaveFig);
if curFlagGenFigsSilently
    close(curProcTimeFig);
    close(curProcTimeHistFig);
end

disp('    Done!')

%% Empirical CDFs

disp(' ')
disp('    Plotting empirical CDFs ...')

% For blockage maps.
curFlagGenFigsSilently = true;
mapType = 'Blockage';

[curEmpCdfFig, simState.blockageMapsCovRatioMeta] ...
    = plotEmpiricalCdfForCoverage(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType, '.png']);
saveas(curEmpCdfFig, pathToSaveFig);
if curFlagGenFigsSilently
    close(curEmpCdfFig);
end

% For coverage maps.
mapType = 'Coverage';

[curEmpCdfFig, simState.coverageMapsCovRatioMeta] ...
    = plotEmpiricalCdfForCoverage(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType]);
saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
saveas(curEmpCdfFig, [pathToSaveFig, '.fig']);
% Coverage ratio gain.
[ curCovRatioGainFig ] ...
    = plotCoverageRatioGain(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType, 'CovRatGain']);
saveEpsFigForPaper(curCovRatioGainFig, pathToSaveFig);
saveas(curCovRatioGainFig, [pathToSaveFig, '.fig']);

if curFlagGenFigsSilently
    close(curEmpCdfFig);
    close(curCovRatioGainFig);
end

% For blockage distance maps.
mapType = 'BlockageDist';

[curEmpCdfFig, simState.blockageDistMapsCovRatioMeta] ...
    = plotEmpiricalCdfForCoverage(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType]);
saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
saveas(curEmpCdfFig, [pathToSaveFig, '.fig']);
% Coverage ratio gain.
[ curCovRatioGainFig ] ...
    = plotCoverageRatioGain(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType, 'CovRatGain']);
saveEpsFigForPaper(curCovRatioGainFig, pathToSaveFig);
saveas(curCovRatioGainFig, [pathToSaveFig, '.fig']);

if curFlagGenFigsSilently
    close(curEmpCdfFig);
    close(curCovRatioGainFig);
end

% This is not very helpful.
if false
    % For excess path loss maps by vegetation.
    mapType = 'PathLossByVeg';

    [curEmpCdfFig, simState.coverageMapsCovRatioMeta] ...
        = plotEmpiricalCdfForCoverage(simState, simConfigs, ...
        mapType, ~curFlagGenFigsSilently);
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType]);
    saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
    saveas(curEmpCdfFig, [pathToSaveFig, '.fig']);
    % Coverage ratio gain.
    [ curCovRatioGainFig ] ...
        = plotCoverageRatioGain(simState, simConfigs, ...
        mapType, ~curFlagGenFigsSilently);
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['EmpiricalCdf_', mapType, 'CovRatGain']);
    saveEpsFigForPaper(curCovRatioGainFig, pathToSaveFig);
    saveas(curCovRatioGainFig, [pathToSaveFig, '.fig']);

    if curFlagGenFigsSilently
        close(curEmpCdfFig);
        close(curCovRatioGainFig);
    end
end

% For coverage path loss maps with vegetation.
mapType = 'PathLossWithVeg';

[curEmpCdfFig, simState.coverageMapsCovRatioMeta] ...
    = plotEmpiricalCdfForCoverage(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType]);
saveEpsFigForPaper(curEmpCdfFig, pathToSaveFig);
saveas(curEmpCdfFig, [pathToSaveFig, '.fig']);
% Coverage ratio gain.
[ curCovRatioGainFig ] ...
    = plotCoverageRatioGain(simState, simConfigs, ...
    mapType, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    ['EmpiricalCdf_', mapType, 'CovRatGain']);
saveEpsFigForPaper(curCovRatioGainFig, pathToSaveFig);
saveas(curCovRatioGainFig, [pathToSaveFig, '.fig']);

if curFlagGenFigsSilently
    close(curEmpCdfFig);
    close(curCovRatioGainFig);
end

disp('    Done!')

%% LoS Coverage Ratio vs Drone Height

disp(' ')
disp(['    Plotting LoS coverage ratio vs drone height ', ...
    'based on blockage maps...'])

[curCovRatioVsHFig, curCovRatios] ...
    = plotCovRatioVsInspectedHeight(simState.blockageMaps, ...
    simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M, ~curFlagGenFigsSilently);
pathToSaveFig = fullfile(pathToSaveResults, ...
    'CoverageRatioVsDroneHeight_Blockage');
saveEpsFigForPaper(curCovRatioVsHFig,  pathToSaveFig);
if curFlagGenFigsSilently
    close(curCovRatioVsHFig);
end

CoverageRatio = curCovRatios;
UavHeight = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;

curCovRatioTable = table(UavHeight, CoverageRatio);
disp(curCovRatioTable);

pathToSaveTab = fullfile(pathToSaveResults, ...
    'CoverageRatioVsDroneHeight_Blockage.csv');
writetable(curCovRatioTable, pathToSaveTab);

disp('    Done!')

end
% EOF