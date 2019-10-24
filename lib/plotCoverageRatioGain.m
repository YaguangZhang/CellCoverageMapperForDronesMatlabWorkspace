function [hFigCovRatGain] = plotCoverageRatioGain( ...
    simState, simConfigs, mapType, flagVisible)
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
%     A case-insensitive string, 'Blockage' or 'Coverage', controling what
%     type of map to consider.
%   - flagVisible
%     An optional boolean for whether to show the resultant figure or not.
%     Default to true.
%
% Outputs:
%   - hFigCovRatGain
%     The handle to the resultant CDF figure.
%
% Yaguang Zhang, Purdue, 10/05/2019

% Set an appropriate figure size for publication.
desiredFigSizeInPixel = [500, 300];

% For plotting.
markers = {'-', '-.', '--', ':'};
lineWidth = 1;
numOfMarkers = length(markers);
infinitePathLossInDb = 1000;

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

% Get information for empirical CDFs. We need a cell, cdfVs, holding the
% empirical cumulative distribution function (cdf) values, evaluated at the
% points in the cooresponding element of the cell cdfXs.
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

% Get the coverage ratio gain compared to the first map. For the
% comparison, we need to construct a shared grid for the path loss values
% and fit the empirical data for each map over the grid.
numOfGridPts = 1000;
allCdfXs = vertcat(cdfXs{:});
% Remove invalid grid points.
allCdfXs = allCdfXs(isfinite(allCdfXs));
gridValues = linspace(min(allCdfXs), max(allCdfXs), numOfGridPts)';

covRatios = cell(numMaps, 1);
% Find the axis following the same procedure as that for the empirical CDF
% plots.
xMax = -inf; xMin = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2);
for idxMap = 1:numMaps
    % Pre-compute the axis ranges.    
    maxPtIdxToShow = find(~isfinite(cdfXs{idxMap}), 1)-1;
    if isempty(maxPtIdxToShow)
        maxPtIdxToShow = length(cdfXs{idxMap});
    end
    xs = [cdfXs{idxMap}(1:maxPtIdxToShow); ...
        infinitePathLossInDb];
    ys = [cdfVs{idxMap}(1:maxPtIdxToShow); ...
        cdfVs{idxMap}(maxPtIdxToShow)];
    
    % Find the last point with y less than 5% and keep a record of the
    % corresponding x. We will use the minimus x for the axis xMin.
    xMax = max(cdfXs{idxMap}(maxPtIdxToShow), xMax);
    curXMinIdx = find(ys<0.05, 1, 'last');
    if ~isempty(curXMinIdx)
        xMin = min(xMin, xs(curXMinIdx));
    end
    
    % Get the coverage ratios.
    boolsValidData = isfinite(cdfXs{idxMap}) & isfinite(cdfVs{idxMap});
    validCdfXs = cdfXs{idxMap}(boolsValidData);
    validCdfVs = cdfVs{idxMap}(boolsValidData);
    [~, indicesUniqueXs] = unique(validCdfXs);
    covRatios{idxMap} ...
        = interp1(validCdfXs(indicesUniqueXs), ...
        validCdfVs(indicesUniqueXs), gridValues);
end

% We will plot the reference track and delete it later to get a consistant
% legend with the empirical CDF plots.
refCdfVs = covRatios{1};

hFigCovRatGain = figure('Visible', flagVisible, ...
    'Position', [0, 0, desiredFigSizeInPixel]);
hold on; set(gca, 'fontWeight', 'bold');
% Extend the final finite result to show a flat line at the end of each
% track.
extendedGridVs = [gridValues; infinitePathLossInDb];
for idxMap = 1:numMaps
    curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
    curCdfVs = covRatios{idxMap}-refCdfVs;
    
    curExtendedCdfVs = [curCdfVs; ...
        curCdfVs(find(isfinite(curCdfVs), 1, 'last'))];
    curBoolsValidPts = isfinite(curExtendedCdfVs);
    
    hsCdf(idxMap) = plot( ...
        extendedGridVs(curBoolsValidPts), ...
        curExtendedCdfVs(curBoolsValidPts).*100, ...
        curMarker, 'LineWidth', lineWidth); %#ok<AGROW>
end
axis auto; curAxisAuto = axis; axis tight; curAxisTight = axis;
curYLimToSet = curAxisTight(3:4) ...
    + (curAxisAuto(3:4)-curAxisTight(3:4)).*0.5;
axis([xMin xMax curYLimToSet]);
grid on; grid minor; box on;
xlabel('Maximum Allowed Path Loss (dB)'); 
ylabel('Coverage Ratio Gain (%)');
hLeg = legend(hsCdf, ...
    arrayfun(@(n) ['UAV at ', num2str(n), ' m'], cdfMeta.rxHeightsInM, ...
    'UniformOutput', false), ...
    'Location', 'NorthWest');
delete(hsCdf(1));
set(hLeg, 'AutoUpdate', 'off');
line(xlim, [0 0], 'Color', 'k');
end
% EOF