function [ hFigCdf, cdfMeta ] = plotEmpiricalCdfForCoverage( ...
    simState, simConfigs, mapType, flagVisible)
%PLOTEMPIRICALCDFFORCOVERAGE Plot the empirical CDFs of the aggregated path
%loss maps of all RX heights.
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
%     A case-insensitive string, 'Blockage' or 'Coverage', controling what
%     type of map to consider.
%   - flagVisible
%     An optional boolean for whether to show the resultant figure or not.
%     Default to true.
%
% Outputs:
%   - hFigCdf
%     The handle to the resultant CDF figure.
%   - cdfMeta
%     A struct holding the meta information for generating the CDF plot. It
%     has fields:
%       - totalNumOfPixelsOnMap
%         A scalar indicating the number of pixels in the whole path loss
%         map.
%       - rxHeightsInM
%         A column vector for the RX heights in meter.
%       - numsOfCovPixels
%         A column vector for the total numbers of covered pixels for each
%         RX height. For the blockage map, we consider any non-NaN path
%         loss value as covered; for the coverage map, we consider any
%         non-NaN path loss value equal to or below
%         simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2) as covered.
%       - coverageRatios
%         A column vector for the coverage ratios
%         (numsOfCovPixels./totalNumOfPixelsOnMap) for each RX height.
%
% Yaguang Zhang, Purdue, 10/05/2019

% Set an appropriate figure size for publication.
flagResizeFigForPublication = evalin('base', 'RESIZE_FIG_FOR_PUBLICATION');
if flagResizeFigForPublication
    desiredFigSizeInPixel = [500, 300];
else
    defaultFigPos =get(0,'defaultfigureposition');
    desiredFigSizeInPixel = defaultFigPos(3:4);
end

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

switch lower(mapType)
    case 'blockage'
        maps = simState.blockageMaps;
    case 'coverage'
        maps = simState.coverageMaps;
    otherwise
        error(['Unsupported mapType ', mapType, '!']);
end

cdfMeta.totalNumOfPixelsOnMap = length(maps{1});
cdfMeta.rxHeightsInM = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;

% Get information for plotting empirical CDFs. We need a cell, cdfVs,
% holding the empirical cumulative distribution function (cdf) values,
% evaluated at the points in the cooresponding element of the cell cdfXs.
numMaps = length(maps);
numsOfCovPixels = deal(nan(numMaps, 1));
[cdfXs, cdfVs] = deal(cell(numMaps, 1));

for idxMap = 1:numMaps
    curM = maps{idxMap};
    assert(length(curM(:))==cdfMeta.totalNumOfPixelsOnMap, ...
        'All maps should have the same number of pixels!');
    
    switch lower(mapType)
        case 'blockage'
            boolsCovPs = ~isnan(curM);
        case 'coverage'
            boolsCovPs = (~isnan(curM)) ...
                &(curM<=simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2));
    end
    numsOfCovPixels(idxMap) = sum(boolsCovPs);
    
    curMNanToInf = curM; curMNanToInf(~boolsCovPs) = inf;
    [cdfVs{idxMap}, cdfXs{idxMap}] = ecdf(curMNanToInf(:));
end

% Output the meta data for the coverage ratio vs RX height plot.
cdfMeta.numsOfCovPixels = numsOfCovPixels;
cdfMeta.coverageRatio ...
    = cdfMeta.numsOfCovPixels./cdfMeta.totalNumOfPixelsOnMap;

% For plotting.
markers = {'-', '-.', '--', ':'};
lineWidth = 1;
numOfMarkers = length(markers);
infinitePathLossInDb = 1000;

hFigCdf = figure('Visible', flagVisible, ...
    'Position', [0, 0, desiredFigSizeInPixel]);
hold on; set(gca, 'fontWeight', 'bold');
xMax = -inf; xMin = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2);
for idxMap = 1:numMaps
    curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
    maxPtIdxToShow = find(~isfinite(cdfXs{idxMap}), 1)-1;
    if isempty(maxPtIdxToShow)
        maxPtIdxToShow = length(cdfXs{idxMap});
    end
    xs = [cdfXs{idxMap}(1:maxPtIdxToShow); ...
        infinitePathLossInDb];
    ys = [cdfVs{idxMap}(1:maxPtIdxToShow); ...
        cdfVs{idxMap}(maxPtIdxToShow)];
    hsCdf(idxMap) = plot(xs, ys, ...
        curMarker, 'LineWidth', lineWidth); %#ok<AGROW>
    xMax = max(cdfXs{idxMap}(maxPtIdxToShow), xMax);
    
    % Find the last point with y less than 5% and keep a record of the
    % corresponding x. We will use the minimus x for the axis xMin.
    curXMinIdx = find(ys<0.05, 1, 'last');
    if ~isempty(curXMinIdx)
        xMin = min(xMin, xs(curXMinIdx));
    end
end
axis tight; curAxis = axis; axis([xMin xMax curAxis(3:4)]);
grid on; grid minor; box on;
xlabel('Path Loss (dB)'); ylabel('Empirical CDF');
legend(hsCdf, ...
    arrayfun(@(n) ['UAV at ', num2str(n), ' m'], cdfMeta.rxHeightsInM, ...
    'UniformOutput', false), ...
    'Location', 'NorthWest');
end
% EOF