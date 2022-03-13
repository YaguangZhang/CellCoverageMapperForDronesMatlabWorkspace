function [hFigCovRatGain] = plotCoverageRatioGain( ...
    simState, simConfigs, mapType, flagVisible, ...
    flagManualXlim, flagShowFSPL, refFsplVs)
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
%   - flagManualXlim
%     An optional boolean for whether to show the resultant figure or not.
%     Default to false. If true, the xlim will be set to
%     simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB.
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
%   - hFigCovRatGain
%     The handle to the resultant CDF figure.
%
% Update 20191112: If flagResizeFigForPublication is set to true in the
% base workspace, the figure will be resized for publication.
%
% Yaguang Zhang, Purdue, 05/07/2021

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
infinitePathLossInDb = 10^6;

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

if ~exist('flagManualXlim', 'var')
    flagManualXlim = false;
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

cdfMeta.totalNumOfPixelsOnMap = length(maps{1});
cdfMeta.rxHeightsInM = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;

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
    assert(length(curM(:))==cdfMeta.totalNumOfPixelsOnMap, ...
        'All maps should have the same number of pixels!');

    switch lower(mapType)
        case {'blockage', 'blockagedist'}
            boolsCovPs = ~isnan(curM);
        case {'coverage', 'pathlossbyveg', 'pathlosswithveg'}
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
numOfGridPts = 10^4;
allCdfXs = vertcat(cdfXs{:});
% Remove invalid grid points.
allCdfXs = allCdfXs(isfinite(allCdfXs));
gridValues = linspace(min(allCdfXs), max(allCdfXs), numOfGridPts)';

covRatios = cell(numMaps, 1);
% Find the axis range.
xMax = -inf;
switch lower(mapType)
    case {'blockage', 'coverage', 'pathlossbyveg', 'pathlosswithveg'}
        xMin = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2);
    case 'blockagedist'
        xMin = inf;
        xMinLog = 1;
        gridValues = sort([xMinLog; gridValues]);
end

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

    % Make sure the jump at zero, if present, is considered.
    if (validCdfXs(1)==0) && (validCdfXs(2)==0)
        indicesUniqueXs(1) = 2;
    end

    if length(indicesUniqueXs)==1 && length(validCdfVs)==2
        % The all 1 case.
        covRatios{idxMap} = ...
            ones(size(gridValues)).*validCdfVs(2);
        covRatios{idxMap}(1) = validCdfVs(1);
    else
        if max(validCdfXs(indicesUniqueXs)) == max(allCdfXs)
            % Current max x reaches the maximum possible value. No need to
            % expolate.
            covRatios{idxMap} ...
                = interp1(validCdfXs(indicesUniqueXs), ...
                validCdfVs(indicesUniqueXs), gridValues);
        else
            covRatios{idxMap} ...
                = interp1( ...
                [validCdfXs(indicesUniqueXs); infinitePathLossInDb], ...
                [validCdfVs(indicesUniqueXs); 1], gridValues);
        end
    end
end

if flagManualXlim
    xMin = min(xMin, simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(1));
    xMax = max(xMax, simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2));
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

    if flagShowFSPL && (idxMap==numMaps)
        curMarker = 'k-';
    else
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
    end

    hsCdf(idxMap) = plot( ...
        extendedGridVs(curBoolsValidPts), ...
        curExtendedCdfVs(curBoolsValidPts).*100, ...
        curMarker, 'LineWidth', lineWidth); %#ok<AGROW>
end
axis auto; curAxisAuto = axis; axis tight; curAxisTight = axis;
curYLimToSet = curAxisTight(3:4) ...
    + (curAxisAuto(3:4)-curAxisTight(3:4)).*0.5;
curYLimToSet(1) = max(curYLimToSet(1), 0);

switch lower(mapType)
    case {'blockage', 'coverage', 'pathlossbyveg', 'pathlosswithveg'}
        xlabel('Maximum Allowed Path Loss (dB)');
        curLoc = 'NorthWest';

        if flagManualXlim
            xlimToSet = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB;
        else
            xlimToSet = [xMin xMax];
        end
        axis([xlimToSet curYLimToSet]);

    case 'blockagedist'
        xlabel('Maximum Allowed Blockage Distance (m)');
        set(gca, 'XScale', 'log');
        curLoc = 'NorthEast';
        axis([xMinLog xMax curYLimToSet]);
end
grid on; grid minor; box on;

ylabel('Coverage Ratio Gain (%)');
if flagShowFSPL
    hLeg = legend(hsCdf, ...
        [arrayfun(@(n) ['Relay at ', num2str(n), ' m'], ...
        cdfMeta.rxHeightsInM, 'UniformOutput', false); ...
        'Reference FSPL'], ...
        'Location', curLoc);
else
    hLeg = legend(hsCdf, ...
        arrayfun(@(n) ['Relay at ', num2str(n), ' m'], ...
        cdfMeta.rxHeightsInM, 'UniformOutput', false), ...
        'Location', curLoc);
end
delete(hsCdf(1));
set(hLeg, 'AutoUpdate', 'off');
line(xlim, [0 0], 'Color', 'k');
% line([xMax xMax], ylim, 'Color', 'k');
end
% EOF