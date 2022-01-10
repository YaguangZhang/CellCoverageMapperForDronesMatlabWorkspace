function [matRxLocWithPathloss, rxAntH, cellAntLoc, cellAntH] ...
    = fetchPathlossResultsFromSimState(simState, idxEffeCell, ...
    idxRxHeight, mapType, simConfigs, locType)
%FETCHPATHLOSSRESULTSFROMSIMSTATE A helper function to fetch pathloss
%results from the simulation state struct variable simState.
%
% Note that when GPS system location is requested, we use (lon, lat) (note
% the order) as the output because it matches the UTM (x, y) in plots.
%
% Inputs:
%   - simState
%     The simulation result struct. We need fields:
%       - CellAntsXyhEffective
%         The UTM (x, y) coordinates and heights for the effective cellular
%         tower antennas in the simulation.
%   - idxEffeCell, idxRxHeight
%     The indices for the effective cellular antenna and inspected RX
%     height of interest.
%   - mapType
%     String. 'Blockage', 'Coverage', or 'BlockageDist'. The first two
%     correspond to the pathloss results stored in fields
%     blockageMapsForEachCell and coverageMapsForEachCell of simState,
%     respectively. And the third one corresponds to the blockage distance
%     results stored in the field blockageDistMapsForEachCell of simState.
%   - simConfigs
%     The configuration struct for the simulation. We need fields:
%       - RX_ANT_HEIGHTS_TO_INSPECT_IN_M
%         Rx heights for the simulator to inspect.
%   - locType
%     Either 'UTM' or 'GPS', corresponding to the UTM (x, y) and GPS (lon,
%     lat) systems, respectively, for the output locations.
%
% Ouputs:
%   - matRxLocWithPathloss
%     A matrix with each row being (x, y, pathloss) or (lon, lat,
%     pathloss), depending on the value of input locType, for one pixel of
%     the specificed type of map via input mapType.
%   - rxAntH
%     The height of the RX antenna.
%   - cellAntLoc, cellAntH
%     The location of the cellular tower antenna and height, respectively.
%     Similarly to the RX location, cellAntLoc is (x, y) or (lon, lat)
%     depending on locType.
%
% Yaguang Zhang, Purdue, 05/07/2021

rxAntH = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M(idxRxHeight);
cellAntXY = simState.CellAntsXyhEffective(idxEffeCell, 1:2);
cellAntH = simState.CellAntsXyhEffective(idxEffeCell, 3);

switch lower(mapType)
    case('blockage')
        values = simState.blockageMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}';
    case('coverage')
        values = simState.coverageMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}';
    case('blockagedist')
        values = simState.blockageDistMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}';
    case('blockagebyterraindist')
        values = simState.blockageByTerrainDistMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}';
    case('blockagebyvegdist')
        values = simState.blockageDistMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}' ...
            - simState.blockageByTerrainDistMapsForEachCell ...
            {idxEffeCell}{idxRxHeight}';
    case('pathlosswithveg')
        values = simState.pathLossWithVegMapsForEachCell...
            {idxEffeCell}{idxRxHeight}';
    otherwise
        error(['Unsupported mapType ', mapType, '!']);
end

switch lower(locType)
    case('utm')
        rxLocs = simState.mapGridXYPts;
        cellAntLoc = cellAntXY;
    case('gps')
        rxLocs = simState.mapGridLatLonPts(:, [2,1]);
        [cellAntLat, cellAntLon] ...
            = simConfigs.utm2deg_speZone(cellAntXY(1), cellAntXY(2));
        cellAntLoc = [cellAntLon, cellAntLat];
    otherwise
        error(['Unsupported locType ', locType, '!']);
end

matRxLocWithPathloss = [rxLocs values];

end
% EOF