% COMBINEPATHLOSSMAPS The snippet for combining path loss maps.
%
% Yaguang Zhang, Purdue, 07/03/2019

disp(' ')
disp('        Combinning path loss maps ...')

[coverageMapsEHata, coverageMapsEHataXLabels, coverageMapsEHataYLabels] ...
    = deal(cell(numOfHs, 1));

% Get the range of path loss values to be shown.
minPathLossValueToShow = inf;
maxPathLossValueToShow = -inf;
for idxH = 1:numOfHs
    curRxAntH = RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
    
    [coverageMapsEHata{idxH}, coverageMapsEHataXLabels{idxH}, ...
        coverageMapsEHataYLabels{idxH}] = combineTowerPathLossMaps( ...
        towerPathLossMapsEHata{idxH}, ...
        towerPathLossMapsEHataXLabels{idxH}, ...
        towerPathLossMapsEHataYLabels{idxH}, ...
        EXPECTED_PL_RANGE_IN_DB);
    
    minPathLossValueToShow = min([minPathLossValueToShow; ...
        coverageMapsEHata{idxH}(:)]);
    maxPathLossValueToShow = max([maxPathLossValueToShow; ...
        coverageMapsEHata{idxH}(:)]);
end

if FLAG_GEN_FIGS
    colorRangeForFinalMaps = [floor(minPathLossValueToShow), ...
        ceil(maxPathLossValueToShow)];
    
    for idxH = 1:numOfHs
        curRxAntH = RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);
        
        % Generate a figure on Google map to show the path loss map.
        [curXs, curYs] ...
            = meshgrid(coverageMapsEHataXLabels{idxH}, ...
            coverageMapsEHataYLabels{idxH});
        [curLats, curLons] ...
            = utm2deg(curXs(:), curYs(:), ...
            repmat(UTM_ZONE, length(curXs(:)), 1));
        
        hCurPLMap = figure;
        hold on;
        hCurTxs = plot3(cellAntsLatLon(:,2), cellAntsLatLon(:,1), ...
            ones(length(cellAntsLatLon(:,1)), 1)...
            .*(max(coverageMapsEHata{idxH}(:))+1), ...
            'xr', 'LineWidth', 1.5);
        plot3k([curLons curLats ...
            coverageMapsEHata{idxH}(:)], ...
            'Labels', {'', ...
            'Longitude (degrees)', 'Latitude (degrees)', ...
            '', 'Path Loss (dB)'}, ...
            'ColorRange', colorRangeForFinalMaps);
        grid on; view(2); axis tight;
        legend(hCurTxs, 'TXs', 'Location', 'SouthEast');
        plotGoogleMapAfterPlot3k(gcf, 'satellite');
        
        pathToSaveFig = fullfile(pathToSaveResults, [ ...
            'Coverage_RxHeight_', num2str(curRxAntH), ...
            '_eHataLib_', LIBRARY_TO_USE, '.png']);
        saveas(hCurPLMap, pathToSaveFig);
    end
end

save(ABS_PATH_TO_SAVE_COVERAGE_MAPS, 'RX_ANT_HEIGHTS_TO_INSPECT_IN_M', ...
    'coverageMapsEHata', 'towerPathLossMapsEHata', ...
    'coverageMapsEHataXLabels', 'coverageMapsEHataYLabels');

disp('    Done!')

% EOF