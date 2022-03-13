function [ hFigCdf, cdfMeta ] = plotEmpiricalCdfForCoverage( ...
    simState, simConfigs, mapType, flagVisible, ...
    flagManualXlim, flagShowFSPL, refFsplVs)
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
%     default be evaluated in 2D (or 3D) based on simState.mapGridXYPts and
%     simState.CellAntsXyhEffective, unless refFsplVs is present.
%   - refFsplVs
%     Optional. If refFsplVs is present and flagShowFSPL is true, this
%     function will skip the FSPL evaluation and use directly refFsplVs.
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
% Update 20191112: If flagResizeFigForPublication is set to true in the
% base workspace, the figure will be resized for publication.
%
% Yaguang Zhang, Purdue, 05/18/2021

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

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

if ~exist('flagManualXlim', 'var')
    flagManualXlim = false;
end

if ~exist('flagShowFSPL', 'var')
    flagShowFSPL = false;
end
% 2D: Distances between the Tx and the Rx will be evaluated only based on
%     (x, y).
% 3D: ... based on (x, y, altitude=elevation+height).
fsplType = '2D';

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

% Get information for plotting empirical CDFs. We need a cell, cdfVs,
% holding the empirical cumulative distribution function (cdf) values,
% evaluated at the points in the cooresponding element of the cell cdfXs.
numMaps = length(maps);
numsOfCovPixels = deal(nan(numMaps, 1));
[cdfXs, cdfVs] = deal(cell(numMaps, 1));

if flagShowFSPL
    if ~exist('refFsplVs', 'var')
        switch lower(fsplType)
            case '2d'
                [~, dists] = dsearchn(simState.CellAntsXyhEffective(:,1:2), ...
                    simState.mapGridXYPts);
            case '3d'
                % TODO.
            otherwise
                error(['Unsupported FSPL evaluation type ', fsplType, '!']);
        end
        refFsplVs = fspl(dists, simConfigs.CARRIER_WAVELENGTH_IN_M);
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

% Output the meta data for the coverage ratio vs RX height plot.
cdfMeta.numsOfCovPixels = numsOfCovPixels;
cdfMeta.coverageRatio ...
    = cdfMeta.numsOfCovPixels./cdfMeta.totalNumOfPixelsOnMap;

% Compute the ranges and points to show in the CDF plots.
xMax = -inf;
switch lower(mapType)
    case {'blockage', 'coverage', 'pathlossbyveg', 'pathlosswithveg'}
        xMin = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2);
    case 'blockagedist'
        xMin = inf;
end

[xs, ys] = deal(cell(numMaps, 1));
for idxMap = 1:numMaps
    maxPtIdxToShow = find(~isfinite(cdfXs{idxMap}), 1)-1;
    if isempty(maxPtIdxToShow)
        maxPtIdxToShow = length(cdfXs{idxMap});
    end

    xs{idxMap} = cdfXs{idxMap}(1:maxPtIdxToShow);
    ys{idxMap} = cdfVs{idxMap}(1:maxPtIdxToShow);

    xMax = max(cdfXs{idxMap}(maxPtIdxToShow), xMax);
    % Find the last point with y less than 5% and keep a record of the
    % corresponding x. We will use the minimus x for the axis xMin.
    curXMinIdx = find(ys{idxMap}<0.05, 1, 'last');
    if ~isempty(curXMinIdx)
        xMin = min(xMin, xs{idxMap}(curXMinIdx));
    end
end

if flagManualXlim
    xMin = min(xMin, simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(1));
    xMax = max(xMax, simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB(2));
end
% For plotting.
markers = {'-', '-.', '--', ':'};
lineWidth = 1;
numOfMarkers = length(markers);

hFigCdf = figure('Visible', flagVisible, ...
    'Position', [0, 0, desiredFigSizeInPixel]);
hold on; set(gca, 'fontWeight', 'bold');

% Plots.
for idxMap = 1:numMaps
    curXs = xs{idxMap};
    curYs = ys{idxMap};

    if curXs(1)>xMin
        curXs = [xMin; curXs]; %#ok<AGROW>
        curYs = [curYs(1); curYs]; %#ok<AGROW>
    end

    if curXs(end)<xMax
        if strcmpi(mapType, 'blockagedist') && (curYs(end)~=1)
            warning(['Map #', num2str(idxMap), ...
                '/', num2str(numMaps), ...
                ': the max Y is expected to be 1!', ...
                ' We have ', num2str(curYs(end)), ' instead.']);
        end
        curXs = [curXs; xMax]; %#ok<AGROW>
        curYs = [curYs; curYs(end)]; %#ok<AGROW>
    end

    if strcmpi(mapType, 'blockagedist')
        % Adjust xs for log scale.
        xMin = 1;
        idxFirstNonZeroX = find(curXs~=0, 1, 'first');
        if ~isempty(idxFirstNonZeroX)
            curXs = [xMin; curXs(idxFirstNonZeroX:end)];

            idxFirstNonZeroY = find(curYs~=0, 1, 'first');
            curYs = [curYs(idxFirstNonZeroY); curYs(idxFirstNonZeroX:end)];
        else
            curXs = [];
            curYs = [];
        end
    end

    if flagShowFSPL && (idxMap==numMaps)
        curMarker = 'k-';
    else
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
    end
    hsCdf(idxMap) = plot(curXs, curYs, ...
        curMarker, 'LineWidth', lineWidth); %#ok<AGROW>
end

axis tight;
if flagManualXlim
    xlimToSet = simConfigs.ALLOWED_PATH_LOSS_RANGE_IN_DB;
else
    xlimToSet = [xMin xMax];
end

try
    xlim(xlimToSet);
catch
    warning(['Unable to set xlim to ', num2str(xlimToSet), '!'])
end
grid on; grid minor; box on;

switch lower(mapType)
    case {'blockage', 'coverage', 'pathlossbyveg', 'pathlosswithveg'}
        xlabel('Path Loss (dB)');
        curLoc = 'NorthWest';
    case 'blockagedist'
        xlabel('Blockage Distance (m)');
        set(gca, 'XScale', 'log');
        curLoc = 'SouthEast';
end

ylabel('Empirical CDF');
if flagShowFSPL
    legend(hsCdf, ...
        [arrayfun(@(n) ['Relay at ', num2str(n), ' m'], ...
        cdfMeta.rxHeightsInM, 'UniformOutput', false); ...
        'Reference FSPL'], ...
        'Location', curLoc);
else
    legend(hsCdf, ...
        arrayfun(@(n) ['Relay at ', num2str(n), ' m'], ...
        cdfMeta.rxHeightsInM, 'UniformOutput', false), ...
        'Location', curLoc);
end
end
% EOF