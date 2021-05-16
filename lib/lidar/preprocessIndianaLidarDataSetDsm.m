function [lidarFileRelDirs, xYBoundryPolygons, lonLatBoundryPolygons] ...
    = preprocessIndianaLidarDataSetDsm(ABS_PATH_TO_LOAD_LIDAR, ...
    DEG2UTM_FCT, UTM2DEG_FCT)
%PREPROCESSINDIANALIDARDATASETDSM Preprocess the digital surface model
%(DSM) LiDAR data set located at ABS_PATH_TO_LOAD_LIDAR.
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
%     (R.ProjectedCRS.Name):
%       - "NAD83(2011) / Indiana East (ftUS)"
%           Data for east Indiana.
%       - "NAD83(2011) / Indiana West (ftUS)"
%           Data for west Indiana.
%     This is required to correctly convert the unit survey foot in the
%     LiDAR data to the metric system.
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
% Yaguang Zhang, Purdue, 05/06/2021

ABS_DIR_TO_SAVE_RESULTS = fullfile(ABS_PATH_TO_LOAD_LIDAR, 'metaInfo.mat');
flagDatasetProcessed = exist(ABS_DIR_TO_SAVE_RESULTS, 'file');

% Set this to be false to reuse history processing results.
FLAG_FORCE_REPROCESSING_DATA = false;
% Set this to be true to generate figures for debugging. Because reusing
% history processing results will skip loading all the data needed for
% plotting, we will not generate figures if the data set of interest is
% already processed and FLAG_FORCE_REPROCESSING_DATA is false.
FLAG_GEN_DEBUG_FIGS = false;
%     (~flagDatasetProcessed) ...
%       || FLAG_FORCE_REPROCESSING_DATA;

% Any LiDAR z value too big or too small will be discarded (set to NaN).
maxAllowedAbsLidarZ = 10^38;

% When fetching the elevation data from USGS, we will get more data than
% the boundary box. This helps avoid NaN interpolation results near the
% edge. Note: 1/3 arc second approximately correspond to 10 meters.
usgsBoxPadInM = 15;

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
            
            % Clear last warning message.
            lastwarn('');
            
            % Make sure the file can be loaded properly.
            try
                historyResult = load(curFullPathToSaveLidarResults, ...
                    'xYBoundryPolygon', 'lonLatBoundryPolygon', ...
                    'getLiDarZFromXYFct');
                xYBoundryPolygon = historyResult.xYBoundryPolygon;
                lonLatBoundryPolygon ...
                    = historyResult.lonLatBoundryPolygon;
            catch err
                disp('            There was an error!')
                dispErr(err);
                warning('The history result .mat file is invalid!');
            end
            
            % Make sure the function getLiDarZFromXYFct in the file has a
            % valid workspace. More specifically, the function handle
            % fctLonLatToLidarStatePlaneXY should be struct with field
            % 'STATE_PLANE_CODE_TIPP'.
            getLiDarZFromXYFctDetails = functions( ...
                historyResult.getLiDarZFromXYFct);
            fctLonLatToLidarStatePlaneXYDetails ...
                = functions(getLiDarZFromXYFctDetails.workspace{1} ...
                .fctLonLatToLidarStatePlaneXY);
            if ~isfield( ...
                    fctLonLatToLidarStatePlaneXYDetails.workspace{1}, ...
                    'STATE_PLANE_CODE_TIPP')
                warning('STATE_PLANE_CODE_TIPP not found!')
            end
            
            % Check whether there is any warning in loading the desired
            % data.
            [warnMsg, ~] = lastwarn;
            if isempty(warnMsg)
                curXYBoundryPolygons{idxPar} = xYBoundryPolygon;
                curLonLatBoundryPolygons{idxPar} = lonLatBoundryPolygon;
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
        = fullPathToSaveLidarResults(indicesFilesToProcess);
    [curXYBoundryPolygons, curLonLatBoundryPolygons] ...
        = deal(cell(numOfFilesToProcess,1));
    
    parfor idxPar = 1:numOfFilesToProcess
        idxF = indicesFilesToProcess(idxPar);
        
        % We will ignore the warning for function handles and polyshape.
        warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
        warning('off', 'MATLAB:polyshape:repairedBySimplify');
        
        % Load our Python module.
        py_addpath(fullfile(pwd, 'lib', 'python'));
        
        try
            tic;
            
            disp(['        File # ', ...
                num2str(idxF), '/', num2str(numLidarFiles), ...
                ': Processing raw LiDAR data ...']);
            
            curFullPathToSaveLidarResults ...
                = fullPathToSaveNewLidarResults{idxPar};
            [curLidarFileParentRelDir, curLidarFileName] ...
                = fileparts(lidarFileRelDirs{idxF}); %#ok<PFBNS>
            
            %====== START OF LIDAR DATA PROCESSING ======
            % Load LiDAR data.
            curLidarFileAbsDir = fullfile(ABS_PATH_TO_LOAD_LIDAR, ...
                curLidarFileParentRelDir, [curLidarFileName, '.tif']);
            [lidarDataImg, R] = readgeoraster(curLidarFileAbsDir);
            lidarDataImg(abs(lidarDataImg(:))>maxAllowedAbsLidarZ) ...
                = nan;
            % Reset samples with exact zero values. Regions with water seem
            % to always have a z value of 0.
            lidarDataImg(lidarDataImg(:)==0) = -inf;
            
            [lidarRasterXLabels, lidarRasterYLabels] ...
                = pixcenters(geotiffinfo(curLidarFileAbsDir));
            
            % Convert survery feet to meter.
            lidarDataImg = distdim(lidarDataImg, 'ft', 'm');
            % Convert raster (row, col) to (lat, lon). First, get the
            % projected coordinate system (CRS).
            coordinateSystem = convertStringsToChars( ...
                lower(R.ProjectedCRS.Name));
            % Remove white space.
            coordinateSystem(isspace(coordinateSystem)) = [];
            
            % Some sanity checks for the data set.
            try
                assert(contains(coordinateSystem, lower('NAD')) ...
                    && contains(coordinateSystem, '83'), ...
                    ['Expecting CRS NAD83! (', ...
                    convertStringsToChars(R.ProjectedCRS.Name), ')']);
                % Check if this is a data set for Indiana.
                assert(contains(coordinateSystem, ...
                    lower('Indiana')), ...
                    ['Expecting Indiana data! (', ...
                    convertStringsToChars(R.ProjectedCRS.Name), ')']);
            catch err
                disp(['        There was an error ', ...
                    'validating the projection name!']);
                dispErr(err);
                
                disp(['            Trying to ', ...
                    'validate the projection by its parameters...']);
                coordinateSystem = validateProjectionByParameter( ...
                    R.ProjectedCRS.ProjectionParameters);
                if strcmp(coordinateSystem, 'unknown')
                    R.ProjectedCRS.ProjectionParameters
                    error('Projection validation failed!')
                end
            end
            
            if contains(coordinateSystem, 'East','IgnoreCase', true)
                % The state plane code for the data. A better name would be
                % `STATE_PLANE_CODE_IN`, but we are reusing history code
                % and this piece of information will not be exposed to the
                % user, so it is fine to keep this field name unchanged.
                STATE_PLANE_CODE_TIPP = 'indiana east';
            elseif contains(coordinateSystem, 'West', ...
                    'IgnoreCase', true)
                STATE_PLANE_CODE_TIPP = 'indiana west';
            else
                error( ...
                    ['Unknown projected coordinate system: ', ...
                    R.ProjectedCRS.Name, '!']);
            end
            [lidarRasterXs, lidarRasterYs] ...
                = meshgrid(lidarRasterXLabels, lidarRasterYLabels);
            [lidarLons, lidarLats] ...
                = sp_proj(STATE_PLANE_CODE_TIPP, 'inverse', ...
                lidarRasterXs(:), lidarRasterYs(:), 'sf');
            
            % Store the new (x,y,z) data.
            lidarLats = lidarLats(:);
            lidarLons = lidarLons(:);
            [lidarXs, lidarYs] ...
                = DEG2UTM_FCT(lidarLats, lidarLons); %#ok<PFBNS>
            
            % Note: the data type is single for elements in lidarZs.
            lidarZs = lidarDataImg(:);
            % Note: the final matrix lidarXYZ has "single" elements, even
            % though both lidarXs and lidarYs are of the type "double".
            %   lidarXYZ = [lidarXs, lidarYs, lidarDataImg(:)];
            
            % Find the polygon boundaries.
            xYBoundryPolygonIndices = boundary(lidarXs, lidarYs);
            xYBoundryPolygon ...
                = polyshape([lidarXs(xYBoundryPolygonIndices), ...
                lidarYs(xYBoundryPolygonIndices)]);
            % Note: the (lon, lat) boundary is generated independently.
            lonLatBoundryPolygonIndices ...
                = boundary(lidarLons, lidarLats);
            lonLatBoundryPolygon ...
                = polyshape(...
                [lidarLons(lonLatBoundryPolygonIndices), ...
                lidarLats(lonLatBoundryPolygonIndices)]);
            
            assert( (xYBoundryPolygon.NumRegions == 1) ...
                && (lonLatBoundryPolygon.NumRegions == 1), ...
                'Generated boundaries should have only one region!');
            
            % Create a function to get LiDAR z from UTM coordinates.
            fctLonLatToLidarStatePlaneXY ...
                = @(lon, lat) sp_proj(STATE_PLANE_CODE_TIPP, ...
                'forward', lon, lat, 'sf');
            getLiDarZFromStatePlaneXYFct = @(spXs, spYs) ...
                interp2(lidarRasterXs, lidarRasterYs, ...
                lidarDataImg, spXs, spYs);
            getLiDarZFromXYFct ...
                = @(xs, ys) genRasterLidarZGetter( ...
                getLiDarZFromStatePlaneXYFct, ...
                fctLonLatToLidarStatePlaneXY, ...
                xs, ys, UTM2DEG_FCT);
            
            disp('            Generating elevation information ...');
            
            % Create the ranges for lat/lon in the UTM space. This helps
            % avoid NaN interpolation results near the edge.
            xRange = [min(lidarXs)-usgsBoxPadInM, ...
                max(lidarXs)+usgsBoxPadInM];
            yRange = [min(lidarYs)-usgsBoxPadInM, ...
                max(lidarYs)+usgsBoxPadInM];
            [latRMin, lonRMin] = UTM2DEG_FCT(xRange(1), yRange(1));
            [latRMax, lonRMax] = UTM2DEG_FCT(xRange(2), yRange(2));
            latRange = [latRMin, latRMax];
            lonRange = [lonRMin, lonRMax];
            
            % For avoid reading in incomplete raw terrain files that are
            % being downloaded by other worker, we will try loading the
            % data a few times.
            numTrials = 0;
            maxNumTrialsAllowed = 3; % 30;
            timeToWaitBeforeTryAgainInS = 0; % 30;
            % Default directory to save USDA data.
            usdaDataDir = 'usgsdata';
            % Backup directory to store USDA data when the data in the
            % default directory do not work (possibly because of the
            % collision of workers trying to store the same data file at
            % the same dir). We will assign one folder for each worker to
            % avoid any future collisions.
            usdaDataBackupDir = fullfile('usgsdata_forEachWorker', ...
                ['worker_', num2str(labindex)]);
            fctFetchRegion = @(latR, lonR) ...
                fetchregion(latR, lonR, ...
                'display', true, 'dataDir', usdaDataDir);
            while ~isinf(numTrials)
                try
                    % Use USGS 1/3 arc-second (~10m) resolution data for US
                    % terrain elevation.
                    region = fctFetchRegion(latRange, lonRange);
                    rawElevData = region.readelevation(...
                        latRange, lonRange, ...
                        'sampleFactor', 1, ...
                        'display', true);
                    numTrials = inf;
                catch err
                    disp('            There was an error!')
                    dispErr(err);
                    
                    if strcmpi(err.identifier, ...
                            'MATLAB:subsassigndimmismatch')
                        % A workaround for a bug in the terrain elevation
                        % libary. Some downloaded tiles have more data than
                        % what the raster needs. Typically, if the size of
                        % the raster is m x n, the data fetched could be
                        % (m+1) x (n+1). We will only use m x n of the data
                        % fetched, then.
                        fctFetchRegion = @(latR, longR) ...
                            fetchAnomalyRegion(latR, longR, ...
                            'display', true, 'dataDir', usdaDataDir);
                    else
                        % Use the backup USGS cache folder to avoid
                        % collision.
                        fctFetchRegion = @(latR, longR) ...
                            fetchregion(latR, longR, ...
                            'display', true, ...
                            'dataDir', usdaDataBackupDir);
                    end
                    
                    numTrials = numTrials+1;
                    warning(['Error fetching elevation info for ', ...
                        'file # ', num2str(idxF), '/', ...
                        num2str(numLidarFiles), ...
                        ' (trial # ', num2str(numTrials), ')!']);
                    if numTrials == maxNumTrialsAllowed
                        error( ...
                            ['Error fetching elevation info for ', ...
                            'file # ', num2str(idxF), '/', ...
                            num2str(numLidarFiles), ...
                            ' (trial # ', num2str(numTrials), ')!']);
                    end
                    pause(timeToWaitBeforeTryAgainInS);
                end
            end
            
            disp(['        Fitting elevation data ', ...
                'in the (lon, lat) system ...']);
            % Order the raw elevation data so that both lat and lon are
            % monotonically increasing.
            [rawElevDataLonsSorted, rawElevDataLatsSorted, ...
                rawElevDataElesSorted] ...
                = sortGridMatrixByXY(...
                rawElevData.longs, rawElevData.lats, ...
                rawElevData.elev);
            
            % Create a grid for the elevation data.
            [usgsElevDataLons, usgsElevDataLats] = meshgrid( ...
                rawElevDataLonsSorted, rawElevDataLatsSorted);
            
            % For very tiny LiDAR data tiles, we may not have enough
            % elevation data to carry out interp2. If that happens, we will
            % use the LiDAR data grid for the elevation, too.
            try
                % Interperlate the data with lat and lon.
                lidarEles = interp2( ...
                    usgsElevDataLons, usgsElevDataLats, ...
                    rawElevDataElesSorted, lidarLons, lidarLats);
                
                % For nan results, fetch the elevation data from USGS.
                boolsNanEles = isnan(lidarEles);
                lidarEles(boolsNanEles) ...
                    = queryElevationPointsFromUsgsInChunks( ...
                    lidarLats(boolsNanEles), lidarLons(boolsNanEles));
                
                % Create a function to get elevation from UTM coordinates
                % in the same way.
                getEleFromXYFct = @(xs, ys) ...
                    genUtmEleGetter( ...
                    usgsElevDataLats, usgsElevDataLons, ...
                    rawElevDataElesSorted, xs, ys, UTM2DEG_FCT);
            catch err
                disp(['        There was an error ', ...
                    'interperlating the elevation data!']);
                dispErr(err);
                
                % Fetch the key stored by plot_google_maps. If this does
                % not work, please run plot_google_maps with your Google
                % Maps API once and try again.
                googleApiKeyFile = load(fullfile( ...
                    fileparts(which('plot_google_map')), ...
                    'api_key.mat'));
                lidarEles = getElevations(lidarLats, lidarLons, ...
                    'key', googleApiKeyFile.apiKey);
                lidarRasterEles ...
                    = reshape(lidarEles, size(lidarRasterXs));
                
                getLiDarZFromStatePlaneXYFct = @(spXs, spYs) ...
                    interp2(lidarRasterXs, lidarRasterYs, ...
                    lidarRasterEles, spXs, spYs);
                getEleFromXYFct ...
                    = @(xs, ys) genRasterLidarZGetter( ...
                    getLiDarZFromStatePlaneXYFct, ...
                    fctLonLatToLidarStatePlaneXY, ...
                    xs, ys, UTM2DEG_FCT);
            end
            
            parsave(curFullPathToSaveLidarResults, ...
                lidarXs, lidarYs, lidarZs, ...
                xYBoundryPolygon, getLiDarZFromXYFct, ...
                lidarLats, lidarLons, lidarEles, getEleFromXYFct, ...
                lonLatBoundryPolygon,  ...
                STATE_PLANE_CODE_TIPP, DEG2UTM_FCT, UTM2DEG_FCT);
            %====== END OF LIDAR DATA PROCESSING ======
            
            curXYBoundryPolygons{idxPar} = xYBoundryPolygon;
            curLonLatBoundryPolygons{idxPar} = lonLatBoundryPolygon;
            
            if FLAG_GEN_DEBUG_FIGS
                close all; %#ok<UNRCH>
                
                lidarZsWaterRm = lidarZs;
                lidarZsWaterRm(isinf(lidarZsWaterRm)) = nan;
                lidarXYZ = [lidarXs, lidarYs, double(lidarZsWaterRm)];
                
                % Boundary on map.
                figure;
                plot(lonLatBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot_google_map('MapType', 'satellite');
                
                % (lidarLons, lidarLats, lidarZs).
                figure; hold on;
                plot(lonLatBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarLons, lidarLats, lidarZsWaterRm]);
                view(2);
                plot_google_map('MapType', 'satellite');
                
                % lidarXYZ.
                figure; hold on;
                plot(xYBoundryPolygon, 'FaceColor','red','FaceAlpha',0.1);
                plot3k(lidarXYZ);
                axis equal; view(2);
                
                % A preview of all the elevation data fetched.
                dispelev(rawElevData, 'mode', 'latlong');
                plot_google_map('MapType', 'satellite');
                
                % lidarZ - getLiDarZFromXYFct(lidarX, lidarY).
                figure; hold on;
                plot(xYBoundryPolygon, 'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), ...
                    lidarZsWaterRm - getLiDarZFromXYFct(lidarXs, lidarYs)]);
                axis equal; view(2);
                
                % lidarEles - getEleFromXYFct(lidarX, lidarY).
                figure; hold on;
                plot(xYBoundryPolygon, 'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), ...
                    lidarEles - getEleFromXYFct(lidarXs, lidarYs)]);
                axis equal; view(2);
                
                % lidarZ - lidarEles.
                figure; hold on;
                plot(xYBoundryPolygon, 'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), ...
                    lidarZsWaterRm - lidarEles]);
                axis equal; view(2);
            end
            
            toc;
            disp('        Done!');
        catch err
            disp('        There was an error!')
            dispErr(err);
            error(...
                ['Error processing LiDAR data for ', ...
                'file # ', num2str(idxF), '/', ...
                num2str(numLidarFiles), '!']')
        end
        
        warning('on', 'MATLAB:dispatcher:UnresolvedFunctionHandle');
        warning('on', 'MATLAB:polyshape:repairedBySimplify');
    end
    
    xYBoundryPolygons(indicesFilesToProcess) ...
        = curXYBoundryPolygons;
    lonLatBoundryPolygons(indicesFilesToProcess) ...
        = curLonLatBoundryPolygons;
    
    save(ABS_DIR_TO_SAVE_RESULTS, ...
        'lidarFileRelDirs', 'xYBoundryPolygons', 'lonLatBoundryPolygons');
end

disp('    Done!')

end
% EOF