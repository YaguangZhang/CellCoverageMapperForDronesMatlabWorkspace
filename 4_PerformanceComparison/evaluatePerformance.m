% EVALUATEPERORMANCE Evaluate and compare performance for different
% coverage maps.
%
% Note that the empirical CDF is essentially coverage ratio, so just the
% first plot generated in this script is good enough for performance
% comparison.
%
% Yaguang Zhang, Purdue, 07/26/2019

clear; clc; close all;

%% Configurations

% Locate the current working directory.
cd(fileparts(mfilename('fullpath')));
[~,folderNameToSaveResults,~]=fileparts(pwd);
cd('..'); addpath('lib');
curFileName = mfilename;
fileNameHintRuler = hintScriptName(curFileName);

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
setPath;

% Paths to saved maps to be processed and their label.
flagTest = false;
if flagTest
    foldersToProcessDic = containers.Map({'acre'}, ...
        {'1_CoverageMapsForAcre'});
else
    foldersToProcessDic = containers.Map({'acre', ...
        'tipp_pathloss', 'tipp_blockage'}, ...
        {'1_CoverageMapsForAcre', ...
        '2_CoverageMapsForTipp', '3_BlockageMapsForTipp'});
end

%% Before Processing the Data

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', folderNameToSaveResults);

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

%% Evaluate Performance

disp('    Evaluating coverage maps ...');

% For plotting.
markers = {'-', '-.', '--', ':'};
lineWidth = 1;
numOfMarkers = length(markers);
maxPathLossInDb = 1000;

numFoldersToProcess = foldersToProcessDic.Count;
dirLabels = foldersToProcessDic.keys;
dirsToProcess = foldersToProcessDic.values;
for idxFolder = 1:numFoldersToProcess
    curLabel = dirLabels{idxFolder};
    curParentFolderName = dirsToProcess{idxFolder};
    curDirToProcess = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', curParentFolderName);
    
    curMapFile = rdir(fullfile(curDirToProcess, '*_CoverageMaps.mat'));
    
    if isempty(curMapFile)
        error(['No coverage map file found for ', curLabel, '!']);
    elseif length(curMapFile)~=1
        error(['Multiple coverage map files found for ', curLabel, '!']);
    end
    
    coverageMapsInfo = load(curMapFile.name);
    heightsInspected = coverageMapsInfo.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;
    maps = coverageMapsInfo.coverageMapsEHata;
    
    % Generate a csv file for coverage ratio.
    numMaps = length(maps);
    [rxHeightInM, totalNumsOfPixels, numsOfCovPixels] ...
        = deal(nan(numMaps, 1));
    
    % Also get information for plotting empirical CDFs.
    [cdfXs, cdfVs] = deal(cell(numMaps, 1));
    
    for idxMap = 1:numMaps
        rxHeightInM(idxMap) = heightsInspected(idxMap);
        
        curM = maps{idxMap};
        totalNumsOfPixels(idxMap) = length(curM(:));
        
        boolsCovPs = ~isnan(curM(:));
        numsOfCovPixels(idxMap) = sum(boolsCovPs);
        
        curMNanToInf = curM; curMNanToInf(~boolsCovPs) = inf;
        [cdfVs{idxMap}, cdfXs{idxMap}] = ecdf(curMNanToInf(:));
    end
    coverageRatio = numsOfCovPixels./totalNumsOfPixels;
    
    % For exporting results.
    curFileToSaveResults = fullfile(pathToSaveResults, ...
        [curLabel, '_performance']);
    
    % Output coverage ratio information to a csv file.
    curCsvFile = [curFileToSaveResults, '_coverageRatio.csv'];
    curData.rxHeightInM = rxHeightInM;
    curData.totalNumsOfPixels = totalNumsOfPixels;
    curData.numsOfCovPixels = numsOfCovPixels;
    curData.coverageRatio = coverageRatio;
    struct2csv(curData, curCsvFile);
    
    % Output empirical CDF information to a figure.
    curPngFile = [curFileToSaveResults, '_empiricalCdf.png'];
    hsCdf = nan(1, numMaps);
    
    hCurCdf = figure; hold on;
    maxX = -inf;
    for idxMap = 1:numMaps
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
        maxPtIdxToShow = find(~isfinite(cdfXs{idxMap}), 1)-1;
        if isempty(maxPtIdxToShow)
            maxPtIdxToShow = length(cdfXs{idxMap});
        end
        xs = [cdfXs{idxMap}(1:maxPtIdxToShow); ...
            maxPathLossInDb];
        ys = [cdfVs{idxMap}(1:maxPtIdxToShow); ...
            cdfVs{idxMap}(maxPtIdxToShow)];
        hsCdf(idxMap) = plot(xs, ys, ...
            curMarker, 'LineWidth', lineWidth);
        maxX = max(cdfXs{idxMap}(maxPtIdxToShow), maxX);
    end
    axis tight; curAxis = axis; axis([curAxis(1) maxX curAxis(3:4)]);
    grid on; grid minor;
    xlabel('Path Loss (dB)'); ylabel('Empirical CDF');
    legend(hsCdf, ...
        arrayfun(@(n) ['RX at ', num2str(n), ' m'], heightsInspected, ...
        'UniformOutput', false), ...
        'Location', 'NorthWest');
    
    saveas(hCurCdf, curPngFile);
end

% Redo the evaluation with a path loss upper bound.
maxValidPathLoss = 120;
for idxFolder = 1:numFoldersToProcess
    curLabel = dirLabels{idxFolder};
    curParentFolderName = dirsToProcess{idxFolder};
    curDirToProcess = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', curParentFolderName);
    
    curMapFile = rdir(fullfile(curDirToProcess, '*_CoverageMaps.mat'));
    
    if isempty(curMapFile)
        error(['No coverage map file found for ', curLabel, '!']);
    elseif length(curMapFile)~=1
        error(['Multiple coverage map files found for ', curLabel, '!']);
    end
    
    coverageMapsInfo = load(curMapFile.name);
    heightsInspected = coverageMapsInfo.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;
    maps = coverageMapsInfo.coverageMapsEHata;
    
    % Generate a csv file for coverage ratio.
    numMaps = length(maps);
    [rxHeightInM, totalNumsOfPixels, numsOfCovPixels] ...
        = deal(nan(numMaps, 1));
    
    % Also get information for plotting empirical CDFs.
    [cdfXs, cdfVs] = deal(cell(numMaps, 1));
    
    for idxMap = 1:numMaps
        rxHeightInM(idxMap) = heightsInspected(idxMap);
        
        curM = maps{idxMap};
        curM(curM>maxValidPathLoss) = nan;
        
        totalNumsOfPixels(idxMap) = length(curM(:));
        
        boolsCovPs = ~isnan(curM(:));
        numsOfCovPixels(idxMap) = sum(boolsCovPs);
        
        curMNanToInf = curM; curMNanToInf(~boolsCovPs) = inf;
        [cdfVs{idxMap}, cdfXs{idxMap}] = ecdf(curMNanToInf(:));
    end
    coverageRatio = numsOfCovPixels./totalNumsOfPixels;
    
    % For exporting results.
    curFileToSaveResults = fullfile(pathToSaveResults, ...
        [curLabel, '_performance']);
    
    % Output coverage ratio information to a csv file.
    curCsvFile = [curFileToSaveResults, ...
        '_coverageRatio_max_', num2str(maxValidPathLoss), 'dB.csv'];
    curData.rxHeightInM = rxHeightInM;
    curData.totalNumsOfPixels = totalNumsOfPixels;
    curData.numsOfCovPixels = numsOfCovPixels;
    curData.coverageRatio = coverageRatio;
    struct2csv(curData, curCsvFile);
    
    % Output empirical CDF information to a figure.
    curPngFile = [curFileToSaveResults, ...
        '_empiricalCdf_max_', num2str(maxValidPathLoss), 'dB.png'];
    hsCdf = nan(1, numMaps);
    
    hCurCdf = figure; hold on;
    maxX = -inf;
    for idxMap = 1:numMaps
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};
        maxPtIdxToShow = find(~isfinite(cdfXs{idxMap}), 1)-1;
        if isempty(maxPtIdxToShow)
            maxPtIdxToShow = length(cdfXs{idxMap});
        end
        xs = [cdfXs{idxMap}(1:maxPtIdxToShow); ...
            maxPathLossInDb];
        ys = [cdfVs{idxMap}(1:maxPtIdxToShow); ...
            cdfVs{idxMap}(maxPtIdxToShow)];
        hsCdf(idxMap) = plot(xs, ys, ...
            curMarker, 'LineWidth', lineWidth);
        maxX = max(cdfXs{idxMap}(maxPtIdxToShow), maxX);
    end
    axis tight; curAxis = axis; axis([curAxis(1) maxX curAxis(3:4)]);
    grid on; grid minor;
    xlabel('Path Loss (dB)'); ylabel('Empirical CDF');
    legend(hsCdf, ...
        arrayfun(@(n) ['RX at ', num2str(n), ' m'], heightsInspected, ...
        'UniformOutput', false), ...
        'Location', 'NorthWest');
    
    saveas(hCurCdf, curPngFile);
end

% Coverage ratio vs max allowable path loss.
maxAllowablePathLosses = 50:150;

% We will use these path loss thresholds as the x axis.
curXs = maxAllowablePathLosses;
numCurXs = length(curXs);
for idxFolder = 1:numFoldersToProcess
    curLabel = dirLabels{idxFolder};
    curParentFolderName = dirsToProcess{idxFolder};
    curDirToProcess = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'PostProcessingResults', curParentFolderName);
    
    curMapFile = rdir(fullfile(curDirToProcess, '*_CoverageMaps.mat'));
    
    if isempty(curMapFile)
        error(['No coverage map file found for ', curLabel, '!']);
    elseif length(curMapFile)~=1
        error(['Multiple coverage map files found for ', curLabel, '!']);
    end
    
    coverageMapsInfo = load(curMapFile.name);
    heightsInspected = coverageMapsInfo.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;
    maps = coverageMapsInfo.coverageMapsEHata;
    
    % For coverage ratios.
    numMaps = length(maps);
    [rxHeightInM, totalNumsOfPixels] ...
        = deal(nan(numMaps, 1));
    coverageRatios = cell(numMaps, 1);
    
    for idxMap = 1:numMaps        
        curCovRats = zeros(numCurXs, 1);
        
        rxHeightInM(idxMap) = heightsInspected(idxMap);
        totalNumsOfPixels(idxMap) = length(curM(:));
        
        curM = maps{idxMap};
        for curIdxX = 1:numCurXs
            curMaxAllowedPL = curXs(curIdxX);
            
            curCovRats(curIdxX) ...
                = sum(curM(:)<=curMaxAllowedPL)./totalNumsOfPixels(idxMap);
        end
        coverageRatios{idxMap} = curCovRats;
    end
    
    % For exporting results.
    curFileToSaveResults = fullfile(pathToSaveResults, ...
        [curLabel, '_performance']);
    
    % Output coverage information to a figure.
    curPngFile = [curFileToSaveResults, ...
        '_coverageVsMaxAllowablePL.png'];
    hsCovRat = nan(1, numMaps);
    
    hCurCovRat = figure; hold on;
    for idxMap = 1:numMaps
        curMarker = markers{mod(idxMap-1, numOfMarkers)+1};        

        xs = curXs;
        ys = coverageRatios{idxMap};
        hsCovRat(idxMap) = plot(xs, ys, ...
            curMarker, 'LineWidth', lineWidth);
    end
    axis tight;
    grid on; grid minor;
    xlabel('Path Loss (dB)'); ylabel('Coverage Ratio');
    legend(hsCovRat, ...
        arrayfun(@(n) ['RX at ', num2str(n), ' m'], heightsInspected, ...
        'UniformOutput', false), ...
        'Location', 'NorthWest');
    
    saveas(hCurCovRat, curPngFile);
end


disp('    Done!');

%% Clear Things Up if Necessary

% Print out a ruler to indicate everything is done.
disp(fileNameHintRuler);

% EOF