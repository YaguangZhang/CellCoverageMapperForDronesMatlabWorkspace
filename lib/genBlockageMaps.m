% GENBLOCKAGEMAPS The snippet for locating blocked locations on maps.
%
% Yaguang Zhang, Purdue, 07/15/2019

disp(' ')
disp('    Generating coverage maps ...')

% Note that although we do not need the eHata model for blockage map
% generation, we still need to initialize the inputs properly.
REGION = 'LoS';
LIBRARY_TO_USE = 'FSPL';

% Generate a blockage map for each height to inspect.
numOfHs = length(RX_ANT_HEIGHTS_TO_INSPECT_IN_M);

if exist(ABS_PATH_TO_SAVE_COMP_PROGRESS, 'file')==2
    disp('        Loading previous results ...')
    load(ABS_PATH_TO_SAVE_COMP_PROGRESS);
    boolResumeCompProg = true;
    flagProcessInterrupted = true;
else
    pathLossMapCompProgress = cell(numOfHs,1);
    % Initializing the labels.
    for idxCompProg = 1:numOfHs
        pathLossMapCompProgress{idxCompProg} = false(numOfCellAnts, 1);
    end
    
    [towerPathLossMapsEHata, towerPathLossMapsEHataXLabels, ...
        towerPathLossMapsEHataYLabels] ...
        = deal(cell(numOfHs, 1));
    execTimeInSecForAllHs = nan(numOfHs,1);
    boolResumeCompProg = false;
end

disp(' ')
disp('        Ray tracing via terrian profiles ...')

for idxH = 1:numOfHs
    disp(' ')
    disp('            Closing figures to save RAM ...')
    close all;
    disp('            Done!')
    disp(' ')
    
    curExecTimeInSecStart = tic;
    
    curRxAntH = RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
    
    if ~boolResumeCompProg
        % Generate a propogation map for each cell tower. We will first
        % generate the path loss maps, one for each cell tower location.
        [pathLossMaps, pathLossMapXLabels, pathLossMapYLabels] ...
            = deal(cell(numOfCellAnts, 1));
    end
    
    for idxCellAntenna = 1:numOfCellAnts
        disp(' ');
        disp(['            Height #', num2str(idxH), '/', ...
            num2str(numOfHs)]);
        disp(['                Cell antenna #', ...
            num2str(idxCellAntenna), '/', num2str(numOfCellAnts)]);
        
        if pathLossMapCompProgress{idxH}(idxCellAntenna)
            % This cell tower is already processed.
            disp(['                    Skipped ', ...
                'because it is already processed.']);
        else
            curBaseAntXY = cellAntsLatLonXYAlt(idxCellAntenna, 3:4);
            % Get the altitude (ground elevation + antenna height) of the
            % base station antenna.
            curBaseAntAltInM = cellAntsLatLonXYAlt(idxCellAntenna, 5);
            % Get the ground elevation for the base station antenna.
            curBaseAntEleInM ...
                = getEleFromXYFct(curBaseAntXY(1), curBaseAntXY(2));
            
            curBaseAntHeightInM = curBaseAntAltInM - curBaseAntEleInM;
            
            [pathLossMaps{idxCellAntenna}, ...
                pathLossMapXLabels{idxCellAntenna}, ...
                pathLossMapYLabels{idxCellAntenna}, ...
                numOfAvailableWorkers] ...
                = genMedianBasicPropLossMaps( ...
                CARRIER_FREQUENCY_IN_MHZ, ...
                curBaseAntXY, curBaseAntHeightInM, ...
                areaOfInterestXs, areaOfInterestYs, curRxAntH, ...
                getEleFromXYFct, getLiDarZFromXYFct, REGION, ...
                TERRAIN_RES_IN_M, LIBRARY_TO_USE, [], inf, [-inf, inf]);
            
            curFigFileName = [ ...
                'CellTowerPathLoss_RxHeight_', ...
                num2str(curRxAntH), ...
                '_TxCell_', num2str(idxCellAntenna), ...
                '_eHataLib_', LIBRARY_TO_USE, ...
                '_Terrain_', REGION];
            
            if FLAG_GEN_FIGS
                % Generate a figure on Google map to show the path loss
                % map.
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
            
            disp(' ');
            disp('                    Creating restore point...');
            pathLossMapCompProgress{idxH}(idxCellAntenna) = true;
            boolResumeCompProg = false;
            
            if exist(ABS_PATH_TO_SAVE_COMP_PROGRESS, 'file')
                save(ABS_PATH_TO_SAVE_COMP_PROGRESS, ...
                    'pathLossMaps', ...
                    'pathLossMapXLabels', ...
                    'pathLossMapYLabels', ...
                    'pathLossMapCompProgress', '-append');
            else
                save(ABS_PATH_TO_SAVE_COMP_PROGRESS, ...
                    'pathLossMaps', ...
                    'pathLossMapXLabels', ...
                    'pathLossMapYLabels', ...
                    'pathLossMapCompProgress', ...
                    'dataTimeStrStart', 'timerValueStart');
            end
            
            disp('                    Done!');
        end
    end
    
    if ~boolResumeCompProg
        % Save the results.
        towerPathLossMapsEHata{idxH} = pathLossMaps;
        towerPathLossMapsEHataXLabels{idxH} = pathLossMapXLabels;
        towerPathLossMapsEHataYLabels{idxH} = pathLossMapYLabels;
        
        save(ABS_PATH_TO_SAVE_COMP_PROGRESS, ...
            'towerPathLossMapsEHata', ...
            'towerPathLossMapsEHataXLabels', ...
            'towerPathLossMapsEHataYLabels', '-append');
    end
    
    if isnan(execTimeInSecForAllHs(idxH))
        execTimeInSecForAllHs(idxH) = toc(curExecTimeInSecStart);
        save(ABS_PATH_TO_SAVE_COMP_PROGRESS, ...
            'execTimeInSecForAllHs', '-append');
    end
end

disp('    Done!')

% EOF