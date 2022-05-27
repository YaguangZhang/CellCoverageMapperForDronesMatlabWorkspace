function [hFigRelayGain, meta] = plotRelayGain( ...
    simState, simConfigs, mapType, flagVisible, ...
    flagShowFSPL, refFsplVs)
%PLOTCOVERAGERATIOGAIN Plot the coverage ratio gain according to the the
%empirical CDF for the path loss map.
%
% Inputs:
%   - simState
%     The struct for the simulation results. We need fields:
%       - blockageMaps, coverageMaps
%         Two cells holding the final aggregated path loss maps, for the
%         blockage maps and the coverage maps, respectively, for all the RX
%         height that have been inspected.
%   - simConfigs
%     The struct for the simulation configurations. We need field:
%       - RX_ANT_HEIGHTS_TO_INSPECT_IN_M
%         Rx heights for the simulator to inspect.
%   - mapType
%     A case-insensitive string, 'Blockage', 'Coverage', or 'BlockageDist',
%     controling what type of map to consider.
%   - flagVisible
%     An optional boolean for whether to show the resultant figure or not.
%     Default to true.
%   - flagShowFSPL
%     An optional boolean for whether to show the ideal FSPL line as
%     reference. Default to false. If true, the best scenario FSPL will by
%     default be evaluated in 2D based on simState.mapGridXYPts and
%     simState.CellAntsXyhEffective, unless refFsplVs is present.
%   - refFsplVs
%     Optional. If refFsplVs is present and flagShowFSPL is true, this
%     function will skip the FSPL evaluation and use directly refFsplVs.
%
% Outputs:
%   - hFigRelayGain
%     The handle to the resultant figure.
%   - meta
%     For debugging. The values for the plot and some relavant parameters.
%
% Yaguang Zhang, Purdue, 05/26/2022

% Set an appropriate figure size for publication.
try
    flagResizeFigForPublication = evalin('base', ...
        'flagResizeFigForPublication');
catch
    flagResizeFigForPublication = false;
end

if flagResizeFigForPublication
    desiredFigSizeInPixel = [500, 300];
else
    defaultFigPos =get(0,'defaultfigureposition');
    desiredFigSizeInPixel = defaultFigPos(3:4);
end

% For plotting.
markers = {'-', '-.', '--', ':'};
lineWidth = 1;
numOfMarkers = length(markers);

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

if ~exist('flagShowFSPL', 'var')
    flagShowFSPL = false;
end

switch lower(mapType)
    case 'blockage'
        maps = simState.blockageMaps;
    case 'coverage'
        maps = simState.coverageMaps;
    case 'blockagedist'
        maps = simState.blockageDistMaps;
    case 'pathlosswithveg'
        maps = simState.pathLossWithVegMaps;
    case 'pathlossbyveg'
        maps = simState.pathLossWithVegMaps - simState.coverageMaps;
    otherwise
        error(['Unsupported mapType ', mapType, '!']);
end

meta.totalNumOfPixelsOnMap = length(maps{1});
meta.rxHeightsInM = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;

% Get information for empirical CDFs. We need a cell, cdfVs, holding the
% empirical cumulative distribution function (cdf) values, evaluated at the
% points in the cooresponding element of the cell cdfXs.
numMaps = length(maps);
numsOfCovPixels = deal(nan(numMaps, 1));
[cdfXs, cdfVs] = deal(cell(numMaps, 1));

if flagShowFSPL
    if ~exist('refFsplVs', 'var')
        refFsplVs = genOptimisticFsplMap(simState);
    end
    maps = [maps, refFsplVs];
    numMaps = numMaps+1;
end

for idxMap = 1:numMaps
    curM = maps{idxMap};
    assert(length(curM(:))==meta.totalNumOfPixelsOnMap, ...
        'All maps should have the same number of pixels!');

    boolsCovPs = ~isnan(curM);
    numsOfCovPixels(idxMap) = sum(boolsCovPs);

    curMNanToInf = curM; curMNanToInf(~boolsCovPs) = inf;
    [cdfVs{idxMap}, cdfXs{idxMap}] = ecdf(curMNanToInf(:));
end

% Get the relay gain compared to the first map. For the comparison, we need
% to construct a shared grid for the target coverage values and fit the
% empirical data for each map over the grid.
numOfGridPts = 10^3;
gridValues = linspace(0, 1, numOfGridPts)';

[ys, relayGains] = deal(cell(numMaps, 1));
for idxMap = 1:numMaps
    % Interpolate the CDF values on the grid. Note that we need to flip the
    % axis.
    ys{idxMap} = interp1(cdfVs{idxMap}, cdfXs{idxMap}, ...
        gridValues, 'linear', nan);

    relayGains{idxMap} = ys{1} - ys{idxMap};
    % Always end at the x axis.
    idxFirstInvalidGain = find(~isfinite(relayGains{idxMap}), 1, 'first');
    assert(all(~isfinite(relayGains{idxMap}(idxFirstInvalidGain:end))), ...
        'NaN values found in the middle of relay gain values!');

    if ~isempty(idxFirstInvalidGain)
        relayGains{idxMap}(idxFirstInvalidGain) = 0;
    end
end

meta.gridValues = gridValues;
meta.ys = ys;
meta.relayGains = relayGains;

hFigRelayGain = figure('Visible', flagVisible, ...
    'Position', [0, 0, desiredFigSizeInPixel]);
hold on; set(gca, 'fontWeight', 'bold');
% Extend the final finite result to show a flat line at the end of each
% track.
for idxMap = 1:numMaps
    if flagShowFSPL && (idxMap==numMaps)
        curMarker = 'k-';
    else
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
    end

    hsRelayGains(idxMap) = plot(gridValues.*100, relayGains{idxMap}, ...
        curMarker, 'LineWidth', lineWidth); %#ok<AGROW>
end
axis tight;

xlabel('Target Coverage Ratio {\it{C}} (%)');
curLegLoc = 'best';
switch lower(mapType)
    case {'blockage', 'coverage', 'pathlossbyveg', 'pathlosswithveg'}
        ylabel('Relay Gain {\it{G_r}} (dB)');
    case 'blockagedist'
        ylabel('Relay Gain in Block Distance (m)');

        set(gca, 'YScale', 'log');
end
grid on; grid minor; box on;

if flagShowFSPL
    hLeg = legend(hsRelayGains, ...
        [arrayfun(@(n) [num2str(n), ' m'], ...
        meta.rxHeightsInM, 'UniformOutput', false); ...
        'Reference FSPL'], ...
        'Location', curLegLoc);
else
    hLeg = legend(hsRelayGains, ...
        arrayfun(@(n) [num2str(n), ' m'], ...
        meta.rxHeightsInM, 'UniformOutput', false), ...
        'Location', curLegLoc);
end
hLegTitle = get(hLeg, 'Title');
set(hLegTitle, 'String', 'Relay height {\it{h_R}}');

delete(hsRelayGains(1));
set(hLeg, 'AutoUpdate', 'off');

end
% EOF