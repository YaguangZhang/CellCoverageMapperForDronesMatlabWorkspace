%GENEXTRAPLOTSFORWHINWEATHERSTATIONS Generate coverage and blockage maps
%via channel simulations.
%
% Yaguang Zhang, Purdue, 10/14/2021

% We will highlight the best some percent of stations.
percentOfStationsToSelect = 10;

% And we will export their information into .csv files.
[stationLats, stationLons, stationIntIds, stationNames] ...
    = loadWhinWeatherStationInfo();

% Find the best stations via empirical CDF.
fThreshold = percentOfStationsToSelect./100;
% For plotting.
areaOfInterestColor = ones(1,3).*0.5;
[gpsLatsBoundaryOfInterest, gpsLonsBoundaryOfInterest] ...
    = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
markerSelected = 'ro';
for idxH = 1:length(simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M)
    curRxAntHInM = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxH);

    % For both blockage distance and path loss maps.
    valueMaps = {simState.blockageDistMaps{idxH}; ...
        simState.coverageMaps{idxH}};
    mapLabels = {'BlockDistInM', 'PathlossInDb'};
    valueLabels = {'Accumulate Blockage Distance (m)', 'Path loss (dB)'};
    for idxM = 1:length(valueMaps)
        curValueMap = valueMaps{idxM};
        curMapLabel = mapLabels{idxM};
        curValueLabel = valueLabels{idxM};

        [f, x] = ecdf(curValueMap);
        idxThreshold = find(f>=fThreshold, 1, 'First');
        xThreshold = interp1(f([idxThreshold-1, idxThreshold]), ...
            x([idxThreshold-1, idxThreshold]), fThreshold);
        indicesBestStations = find(curValueMap<=xThreshold);
        bestStationIntIds = stationIntIds(indicesBestStations);
        bestStationNames = stationNames(indicesBestStations);
        bestStationLats = stationLats(indicesBestStations);
        bestStationLons = stationLats(indicesBestStations);
        bestStationVs = curValueMap(indicesBestStations);

        hFigCdf = figure; hold on; plot(x, f, '.b-');
        axis([x(1), x(end), 0, 1]);
        plot([0, x(end)+1], [fThreshold, fThreshold], 'r--');
        grid on; grid minor; ylabel('Empirical CDF');
        xlabel(curValueLabel);
        saveas(hFigCdf, fullfile(pathToSaveResults, [ 'best', ...
            num2str(percentOfStationsToSelect), 'pct-ecdf-', ...
            strrep(num2str(curRxAntHInM), '.', '_'), 'm-', ...
            curMapLabel, '.jpg']));

        hMap = figure; hold on;
        plot(polyshape( ...
            [gpsLonsBoundaryOfInterest, gpsLatsBoundaryOfInterest]), ...
            'FaceColor', areaOfInterestColor);
        plot3k([simState.mapGridLatLonPts(:,2), ...
            simState.mapGridLatLonPts(:,1), curValueMap]);
        hSelectedStations ...
            = plot(simState.mapGridLatLonPts(indicesBestStations, 2), ...
            simState.mapGridLatLonPts(indicesBestStations, 1), ...
            markerSelected);
        legend(hSelectedStations, ['Best ', ...
            num2str(percentOfStationsToSelect), '% Stations'], ...
            'Location', 'northwest');
        view(2); plot_google_map('MapType', 'roadmap');
        xlabel('Longitude'); ylabel('Latitude');
        title(curValueLabel);
        saveas(hMap, fullfile(pathToSaveResults, [ 'best', ...
            num2str(percentOfStationsToSelect), 'pct-map_view-', ...
            strrep(num2str(curRxAntHInM), '.', '_'), 'm-', ...
            curMapLabel, '.jpg']));

        % Sort the records by their integer IDs.
        [bestStationIntIds, indicesNewOrder] = sort(bestStationIntIds);
        bestStationNames = bestStationNames(indicesNewOrder);
        bestStationLats = bestStationLats(indicesNewOrder);
        bestStationLons = bestStationLons(indicesNewOrder);
        eval(['bestStation', curMapLabel, ...
            ' = bestStationVs(indicesNewOrder);']);

        writetable(eval(['table(bestStationIntIds, bestStationNames, ', ...
            'bestStationLats, bestStationLons, ', ...
            'bestStation', curMapLabel, ')']), ...
            fullfile(pathToSaveResults, [ 'best', ...
            num2str(percentOfStationsToSelect), 'pct-station_info-', ...
            strrep(num2str(curRxAntHInM), '.', '_'), 'm-', ...
            curMapLabel, '.csv']));
    end
end
% EOF