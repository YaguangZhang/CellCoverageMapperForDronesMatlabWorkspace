function plotCoverageSimulationResults( ...
    pathToSaveResults, simState, simConfigs)
%PLOTCOVERAGESIMULATIONRESULTS Generate plots for the simulation results
%from simulateChannelForExtendedTipp.m.
%
% Inputs:
%   - pathToSaveResults
%     The full path to the directory for saving plots.
%   - simState, simConfigs
%     The outputs from simulateChannelForExtendedTipp.m representing the
%     simulation results and simulaiton configurations, respectively.
%
% Yaguang Zhang, Purdue, 10/14/2019

% Parameters.
numOfFigsToGenForSingleCellInspection = 50;
[numOfEffeCellAnts, ~] = size(simState.CellAntsXyhEffective);
numOfRxHeightToInspect = length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M);

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
    indicesOfFigTasks = 1:numOfFigsToGenForSingleCellInspection;
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

pathToSaveFig = fullfile(pathToSaveResults, ...
    'Overview_SimulationSetup.png');
saveas(hCurPLMap,  pathToSaveFig);
close(hCurPLMap);

% Smaller maps for publication.
if ~isfield(simState,'flagResizeFigForPublication')
    % By default, we do not need to resize figures.
    evalin('base', 'flagResizeFigForPublication = false');
else
    evalin('base', ...
        'flagResizeFigForPublication=simState.RESIZE_FIG_FOR_PUBLICATION');
end

flagResizeFigForPublication ...
    = evalin('base', 'flagResizeFigForPublication');
if flagResizeFigForPublication
    customFigSize = [500, 500].*0.6;
else
    defaultFigPos = get(0, 'defaultfigureposition');
    customFigSize = defaultFigPos(3:4);
end

for idxH = 1:numOfRxHeightToInspect
    rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
    
    % For blockage maps.
    curFlagZoomIn = true;
    % Use plot3k to avoid fitting data for NaN spots.
    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, simState.blockageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, 'plot3k');
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_Blockage_RxHeight_', num2str(rxAntH), '.png']);
    saveas(hCurPLMap,  pathToSaveFig);
    close(hCurPLMap);
    
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
    end
    % Further tighten the figure.
    xlabel(''); ylabel('');
    tightfig(hCurBlMap);
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['BlockageMap_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurBlMap, [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurBlMap, [pathToSaveFig, '.png']);
    close(hCurBlMap);
    
    % For coverage maps.
    curFlagZoomIn = true;
    
    [ hCurPLMap ] ...
        = plotPathLossMap( ...
        [mapGridLonLats, simState.coverageMaps{idxH}], ...
        effeCellAntLonLats, simConfigs, ~curFlagGenFigsSilently, ...
        curFlagZoomIn, defaultCmdToPlotPLMaps, customFigSize);
    hCb = findobj(hCurPLMap, 'Type', 'Colorbar');
    hLeg = findobj(hCurPLMap, 'Type', 'Legend');
    switch lower(simConfigs.CURRENT_SIMULATION_TAG)
        case 'extendedtipp'
            set(hLeg, 'Position', [0.1075, 0.8640, 0.3745, 0.0590]);
        case 'tipp'
            set(hLeg, 'Position', [0.4053, 0.1144, 0.3389, 0.0629]);
            set(hCb, 'Visible', 'off');
    end
    % Further tighten the figure.
    xlabel(''); ylabel('');
    tightfig(hCurPLMap);
    % Move the colorbar title a bit to fully show it.
    hCbTitle = get(hCb, 'Title');
    set(hCbTitle, 'Unit', 'pixel');
    curCbTitlePos = hCbTitle.Position;
    set(hCbTitle, 'Position', [5 curCbTitlePos(2)+3 curCbTitlePos(3)]);
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['PathLossMap_Coverage_RxHeight_', ...
        strrep(num2str(rxAntH), '.', '_')]);
    saveas(hCurPLMap,  [pathToSaveFig, '.eps'], 'epsc');
    saveas(hCurPLMap,  [pathToSaveFig, '.png']);
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

disp('    Done!')

%% Coverage Ratio vs Drone Height

disp(' ')
disp(['    Plotting coverage ratio vs drone height ', ...
    'for the blocakge maps...'])

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