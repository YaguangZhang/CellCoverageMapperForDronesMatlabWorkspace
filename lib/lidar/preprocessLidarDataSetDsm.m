function [lidarFileRelDirs, xYBoundryPolygons, lonLatBoundryPolygons] ...
    = preprocessLidarDataSetDsm(ABS_PATH_TO_LOAD_LIDAR, ...
    DEG2UTM_FCT, UTM2DEG_FCT, FLAG_FORCE_REPROCESSING_DATA) %#ok<INUSL>
%PREPROCESSLIDARDATASETDSM Preprocess the digital surface model (DSM) LiDAR
%data set located at ABS_PATH_TO_LOAD_LIDAR.
%
% This is a modified version of preprocessIndianaLidarDataSetDsm.m. The key
% difference is that the data are not limited to Indiana State, as long as
% the projection information is properly provided in their spatial
% referencing information (R) under the projected coordinate system used
% (R.ProjectedCRS).
%
% We will load the LiDAR data, obtain the sample locations, obtain the
% elevation data for them, and save the results in .mat files.
%
% Inputs:
%   - ABS_PATH_TO_LOAD_LIDAR
%     The absolute path to the .tif LiDAR data obtained from
%       https://lidar.jinha.org/
%     Note that the LiDAR data should include in their spatial referencing
%     information (R) the projected coordinate system used
%     (R.ProjectedCRS). This is required to correctly convert the unit
%     survey foot in the LiDAR data to the metric system.
%   - DEG2UTM_FCT, UTM2DEG_FCT
%     The functions to use to convert (lat, lon) to UTM (x, y) and back,
%     respectively, i.e.: (x, y) = DEG2UTM_FCT(lat, lon); (lat, lon) =
%     UTM2DEG_FCT(x, y).
%
% Outputs:
%   - lidarFileRelDirs
%     A column cell with the relative paths (relative to
%     ABS_PATH_TO_LOAD_LIDAR) for .img LiDAR data files processed. For
%     example,
%       fullfile(ABS_PATH_TO_LOAD_LIDAR, lidarFileRelDirs{1})
%     will output the absolute path to the first LiDAR data file processed.
%   - xYBoundryPolygons, lonLatBoundryPolygons
%     The polygon boundries (Matlab built-in polyshape) for the files
%     processed indicating the area they cover, in terms of UTM (x,y) and
%     GPS (lon, lat), respectively.
%
% Update 20210506: Break lidarXYZ to lidarXs, lidarYs, and lidarZs to avoid
% data type conversion.
%
% Yaguang Zhang, Purdue, 03/25/2022

ABS_DIR_TO_SAVE_RESULTS = fullfile(ABS_PATH_TO_LOAD_LIDAR, 'metaInfo.mat');
flagDatasetProcessed = exist(ABS_DIR_TO_SAVE_RESULTS, 'file');

% Set this to be false to reuse history processing results.
if ~exist('FLAG_FORCE_REPROCESSING_DATA', 'var')
    FLAG_FORCE_REPROCESSING_DATA = false;
end

% Set this to be true to generate figures for debugging. Because reusing
% history processing results will skip loading all the data needed for
% plotting, we will not generate figures if the data set of interest is
% already processed and FLAG_FORCE_REPROCESSING_DATA is false.
FLAG_GEN_DEBUG_FIGS = false; %#ok<NASGU>
%     (~flagDatasetProcessed) ...
%       || FLAG_FORCE_REPROCESSING_DATA;

% Any LiDAR z value too big or too small will be discarded (set to NaN).
maxAllowedAbsLidarZ = 10^38; %#ok<NASGU>

% When fetching the elevation data from USGS, we will get more data than
% the boundary box. This helps avoid NaN interpolation results near the
% edge. Note: 1/3 arc second approximately correspond to 10 meters.
usgsBoxPadInM = 15; %#ok<NASGU>

% Display more details for warnings.
warning on verbose;
% Throw an error instead of issuing a warning when workers die.
errIdWorkerAborted = 'MATLAB:remoteparfor:ParforWorkerAborted';
warning('error', errIdWorkerAborted);

[~, datasetName] = fileparts(ABS_PATH_TO_LOAD_LIDAR);

disp(' ')
disp(['    Preprocessing Indiana LiDAR dataset ', datasetName, ' ...'])

if flagDatasetProcessed ...
        && (~FLAG_FORCE_REPROCESSING_DATA)
    disp('        The specified dataset has been processed before.')
    disp('        Loading history results ...')
    load(ABS_DIR_TO_SAVE_RESULTS, ...
        'lidarFileRelDirs', 'xYBoundryPolygons', 'lonLatBoundryPolygons');
else
    tifFileHandles = rdir(fullfile( ...
        ABS_PATH_TO_LOAD_LIDAR, '**', '*.tif'), ...
        '', ABS_PATH_TO_LOAD_LIDAR);
    % Randomize the file list to balance the work load for each parfor
    % chunk. This can be done because when ordered by name, files close to
    % each other tend to require similar amount of resource such as RAM and
    % some of them require a lot more than others.
    rng(1);
    tifFileHandles = tifFileHandles(randperm(length(tifFileHandles)));
    lidarFileRelDirs = {tifFileHandles(:).name}';

    numLidarFiles = length(lidarFileRelDirs);

    [xYBoundryPolygons, lonLatBoundryPolygons] ...
        = deal(cell(numLidarFiles,1));

    % Generate full path strings to the Matlab .mat cache files for saving
    % and reusing the preprocessing results.
    curDirToSaveLidarResults = fullfile( ...
        ABS_PATH_TO_LOAD_LIDAR, ...
        'MatlabCache');
    if ~exist(curDirToSaveLidarResults, 'dir')
        mkdir(curDirToSaveLidarResults)
    end
    fullPathToSaveLidarResults = cell(numLidarFiles, 1);
    for idxF = 1:numLidarFiles
        [~, curLidarFileName] = fileparts(lidarFileRelDirs{idxF});
        fullPathToSaveLidarResults{idxF} = fullfile( ...
            curDirToSaveLidarResults, [curLidarFileName, '.mat']);
    end

    % Load results for files that are already processed.
    boolsFilesProcessed = cellfun(@(matPath) exist(matPath, 'file'), ...
        fullPathToSaveLidarResults)';
    indicesFilesProcessed = find(boolsFilesProcessed);
    numOfFilesProcessed = length(indicesFilesProcessed);

    fullPathToProcessedLidarResults ...
        = fullPathToSaveLidarResults(indicesFilesProcessed);
    curBoolsFilesProcessed = false(1, numOfFilesProcessed);
    [curXYBoundryPolygons, curLonLatBoundryPolygons] ...
        = deal(cell(numOfFilesProcessed,1));

    parfor idxPar = 1:numOfFilesProcessed
        tic;

        disp(['        File # ', ...
            num2str(indicesFilesProcessed(idxPar)), ...
            '/', num2str(numLidarFiles), ...
            ': Loading history results ...']);

        curFullPathToSaveLidarResults ...
            = fullPathToProcessedLidarResults{idxPar};

        % Try reusing the history results.
        if exist(curFullPathToSaveLidarResults, 'file') ...
                && (~FLAG_FORCE_REPROCESSING_DATA)
            % Try loading the .mat file.
            [xYBoundryPolygonCache, lonLatBoundryPolygonCache, ...
                warnMsgCache] = loadCachedLidarDsmTile( ...
                curFullPathToSaveLidarResults);

            if isempty(warnMsgCache)
                curXYBoundryPolygons{idxPar} = xYBoundryPolygonCache;
                curLonLatBoundryPolygons{idxPar} ...
                    = lonLatBoundryPolygonCache;
                curBoolsFilesProcessed(idxPar) = true;
            else
                warning('Failed in loading history data!');
                disp('            Aborted.');
            end
        end

        toc;
    end

    boolsFilesProcessed(indicesFilesProcessed) = curBoolsFilesProcessed;
    xYBoundryPolygons(indicesFilesProcessed) ...
        = curXYBoundryPolygons;
    lonLatBoundryPolygons(indicesFilesProcessed) ...
        = curLonLatBoundryPolygons;

    % Processing other files.
    indicesFilesToProcess = find(~boolsFilesProcessed);
    numOfFilesToProcess = length(indicesFilesToProcess);

    fullPathToSaveNewLidarResults ...
        = fullPathToSaveLidarResults(indicesFilesToProcess); %#ok<NASGU>
    [curXYBoundryPolygons, curLonLatBoundryPolygons] ...
        = deal(cell(numOfFilesToProcess,1));

    load(fullfile(fileparts(mfilename('fullpath')), ...
        'anomalyTileNameWithProjCRS.mat'), ...
        'anomalyTileNameWithProjCRS'); %#ok<NASGU>

    % The workers may die in the parfor loop, which will cause reruns of
    % the whole loop and go through files that are already processeed. We
    % will (i) check whether .mat cache files are already available to
    % avoid unnecessary reprocessing attempts, and (ii) break the jobs into
    % chunks and restart the pool if free RAM is too low.
    curCluster = gcp;
    maxNumOfWorkers = curCluster.NumWorkers;

    [numOfWorkersToUseStart, curNumOfWorkersToUse] ...
        = deal(max(floor(maxNumOfWorkers*0.9), 1));
    numOfFsPerChunk = curNumOfWorkersToUse*3;

    for idxChunk = 1:ceil(numOfFilesToProcess/numOfFsPerChunk)
        chunkStart = 1 + numOfFsPerChunk*(idxChunk-1); %#ok<NASGU>
        chunkEnd = min(numOfFsPerChunk*idxChunk, ...
            numOfFilesToProcess); %#ok<NASGU>

        % Restart the pool if free RAM is too low.
        guardMemAvailOnLinux;
        % Reset the number of workers to use if it is low.
        if curNumOfWorkersToUse<numOfWorkersToUseStart
            curNumOfWorkersToUse = numOfWorkersToUseStart;
        end

        flagCurrentChunkDone = false;
        while ~flagCurrentChunkDone
            try
                parforProcLidarTiles;
                flagCurrentChunkDone = true;
                curNumOfWorkersToUse = min(curNumOfWorkersToUse+1, ...
                    maxNumOfWorkers);
            catch err
                warning('Error in parfor!');
                disp(err);

                % Aggressively decrease the number of workers.
                curNumOfWorkersToUse = max(curNumOfWorkersToUse-5, 1);

                % Consider restarting the pool if any worker dies.
                if strcmp(err.identifier, errIdWorkerAborted)
                    curCluster = gcp('nocreate');
                    % Restart pool only when available worker is (i) less
                    % than max and (ii) less than needed.
                    if (curCluster.NumWorkers<maxNumOfWorkers) ...
                            && (curCluster.NumWorkers<curNumOfWorkersToUse)
                        delete(curCluster); gcp;
                    end
                end
            end
        end
    end

    xYBoundryPolygons(indicesFilesToProcess) ...
        = curXYBoundryPolygons;
    lonLatBoundryPolygons(indicesFilesToProcess) ...
        = curLonLatBoundryPolygons;

    % Sort the result by tile name before saving it.
    [lidarFileRelDirs, indicesNewOrder] = sort(lidarFileRelDirs);
    xYBoundryPolygons = xYBoundryPolygons(indicesNewOrder);
    lonLatBoundryPolygons = lonLatBoundryPolygons(indicesNewOrder);

    save(ABS_DIR_TO_SAVE_RESULTS, ...
        'lidarFileRelDirs', 'xYBoundryPolygons', 'lonLatBoundryPolygons');

    warning('on', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
    warning('on', 'MATLAB:polyshape:repairedBySimplify');
    warning('on', 'MATLAB:MKDIR:DirectoryExists');
end

disp('    Done!')

end
% EOF