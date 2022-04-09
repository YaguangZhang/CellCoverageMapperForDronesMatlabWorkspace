%PARFORPROCLIDARTILE A snippet to process LiDAR tiles via parfor.
%
% The code is put into one script for a clearer structure of
% preprocessLidarDataSetDsm.m
%
% Yaguang Zhang, Purdue, 03/27/2022

if ~exist('fileNameHintRuler', 'var')
    fileNameHintRuler = ' ------------------------- ';
end
if ~exist('datetimeFormat', 'var')
    datetimeFormat = 'yyyy/mm/dd HH:MM:ss';
end

disp(' ')
disp(fileNameHintRuler)
disp(['[', datestr(now, datetimeFormat), ...
    '] Started new round of parforProcLidarTiles.'])
disp(['    Current curNumOfWorkersToUse = ', ...
    num2str(curNumOfWorkersToUse), ' ...'])
disp(fileNameHintRuler)

% For debugging:
%   for idxPar = chunkStart:chunkEnd
parfor (idxPar = chunkStart:chunkEnd, curNumOfWorkersToUse)
    tic;

    idxF = indicesFilesToProcess(idxPar);
    curFullPathToSaveLidarResults ...
        = fullPathToSaveNewLidarResults{idxPar};

    % Try reusing the history results.
    flagReuseCacheSuccess = false;
    if exist(curFullPathToSaveLidarResults, 'file') ...
            && (~FLAG_FORCE_REPROCESSING_DATA)
        % Try loading the .mat file specified by
        % curFullPathToSaveLidarResults.
        [xYBoundryPolygon, lonLatBoundryPolygon, warnMsg] ...
            = loadCachedLidarDsmTile( ...
            curFullPathToSaveLidarResults);

        if isempty(warnMsg)
            flagReuseCacheSuccess = true;
        else
            warning('Failed in loading history data!');
            disp('            Aborted.');
        end
    end

    if ~flagReuseCacheSuccess
        try
            % We will ignore the warning for function handles, polyshape,
            % and fetchregion.
            warning('off', ...
                'MATLAB:dispatcher:UnresolvedFunctionHandle');
            warning('off', 'MATLAB:polyshape:repairedBySimplify');
            warning('off', 'MATLAB:MKDIR:DirectoryExists');

            % Load our Python module.
            py_addpath(fullfile(pwd, 'lib', 'python'));

            disp(['        File # ', ...
                num2str(idxF), '/', num2str(numLidarFiles), ...
                ': Processing raw LiDAR data ...']);

            [curLidarFileParentRelDir, curLidarFileName] ...
                = fileparts(lidarFileRelDirs{idxF}); %#ok<PFBNS>

            %====== START OF LIDAR DATA PROCESSING ======
            % Load LiDAR data.
            curLidarFileAbsDir = fullfile( ...
                ABS_PATH_TO_LOAD_LIDAR, ...
                curLidarFileParentRelDir, ...
                [curLidarFileName, '.tif']);
            [lidarDataImg, spatialRef] ...
                = readgeoraster(curLidarFileAbsDir);
            lidarDataImg(abs( ...
                lidarDataImg(:))>maxAllowedAbsLidarZ) = nan;
            % Reset samples with exact zero values. Regions with water seem
            % to always have a z value of 0.
            lidarDataImg(lidarDataImg(:)==0) = -inf;

            % Convert survery feet to meter.
            lidarDataImg = distdim(lidarDataImg, 'ft', 'm');

            % Convert raster (row, col) to (lat, lon).
            try
                % Essentailly meshgrid matrices.
                [lidarRasterXs, lidarRasterYs] ...
                    = worldGrid(spatialRef);
                % Column vectors.
                [lidarLats, lidarLons] = projinv( ...
                    spatialRef.ProjectedCRS, ...
                    lidarRasterXs(:), lidarRasterYs(:));
            catch err
                disp(['        There was an error converting ', ...
                    'raster (row, col) to (lat, lon)!']);
                disp('        Current spatialRef.ProjectedCRS:');
                disp(spatialRef.ProjectedCRS);
                disp(err);
                disp(' ');
                disp(['        We will use projection ', ...
                    'reference info from other tiles ...']);
                % Convert the coordinates with projection information from
                % other tiles.
                [~, curTileName] = fileparts(curLidarFileAbsDir);
                idxProjToUse = find(cellfun( ...
                    @(n) contains(n, curTileName), ...
                    anomalyTileNameWithProjCRS(:,1)), ...
                    1, 'first'); %#ok<PFBNS>
                spatialRef.ProjectedCRS ...
                    = anomalyTileNameWithProjCRS{idxProjToUse,2};
                [lidarLats, lidarLons] = projinv( ...
                    spatialRef.ProjectedCRS, ...
                    lidarRasterXs(:), lidarRasterYs(:));
            end

            % Store the new (x,y,z) data.
            [lidarXs, lidarYs] ...
                = DEG2UTM_FCT(lidarLats, lidarLons);

            % Note: the data type is single for elements in lidarZs.
            lidarZs = lidarDataImg(:);
            % Note: the final matrix lidarXYZ has "single" elements, even
            % though both lidarXs and lidarYs are of the type "double".
            %   lidarXYZ = [lidarXs, lidarYs, lidarDataImg(:)];

            % Find the polygon boundaries.
            xYBoundryPolygonIndices = boundary(lidarXs, lidarYs);
            xYBoundryPolygon = polyshape( ...
                [lidarXs(xYBoundryPolygonIndices), ...
                lidarYs(xYBoundryPolygonIndices)]);
            % Note: the (lon, lat) boundary is based on the UTM (x, y)
            % result to avoid running the slow command boundary.
            lonLatBoundryPolygon = polyshape(...
                [lidarLons(xYBoundryPolygonIndices), ...
                lidarLats(xYBoundryPolygonIndices)]);

            assert( (xYBoundryPolygon.NumRegions == 1) ...
                && (lonLatBoundryPolygon.NumRegions == 1), ...
                'Generated boundaries should have only 1 region!');

            % Create a function to get LiDAR z from UTM coordinates.
            fctLonLatToLidarStatePlaneXY ...
                = @(lon, lat) projfwd( ...
                spatialRef.ProjectedCRS, lat, lon);
            getLiDarZFromStatePlaneXYFct = @(spXs, spYs) ...
                interp2(lidarRasterXs, lidarRasterYs, ...
                lidarDataImg, spXs, spYs);
            getLiDarZFromXYFct ...
                = @(xs, ys) genRasterLidarZGetter( ...
                getLiDarZFromStatePlaneXYFct, ...
                fctLonLatToLidarStatePlaneXY, ...
                xs, ys, UTM2DEG_FCT);

            disp('            Generating ele information ...');

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
            usdaDataBackupDir = fullfile( ...
                'usgsdata_forEachWorker', ...
                ['worker_', num2str(labindex)]);
            fctFetchRegion = @(latR, lonR) ...
                ... fetchregion(latR, lonR, ...
                fetchAnomalyRegion(latR, lonR, ...
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
                    disp('            Error from fetchRegion!')
                    dispErr(err);

                    if ~strcmpi(err.identifier, ...
                            'MATLAB:subsassigndimmismatch')
                        % Use the backup USGS cache folder to avoid
                        % collision.
                        fctFetchRegion = @(latR, lonR) ...
                            ... fetchregion(latR, lonR, ...
                            fetchAnomalyRegion(latR, lonR, ...
                            'display', true, ...
                            'dataDir', usdaDataBackupDir);
                    end

                    numTrials = numTrials+1;
                    warning(['fetchRegion: ' ...
                        'Failed fetching elevation info ', ...
                        'for file # ', num2str(idxF), '/', ...
                        num2str(numLidarFiles), ...
                        ' (trial # ', num2str(numTrials), ')!']);
                    if numTrials == maxNumTrialsAllowed
                        error(['fetchRegion: ' ...
                            'Error fetching elevation info ', ...
                            'for file # ', num2str(idxF), '/', ...
                            num2str(numLidarFiles), ...
                            ' (trial # ', ...
                            num2str(numTrials), ')!']);
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
                    lidarLats(boolsNanEles), ...
                    lidarLons(boolsNanEles));

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
                lidarLats, lidarLons, ...
                lidarEles, getEleFromXYFct, ...
                lonLatBoundryPolygon, spatialRef, ...
                DEG2UTM_FCT, UTM2DEG_FCT);

            if FLAG_GEN_DEBUG_FIGS
                close all;

                lidarZsWaterRm = lidarZs;
                lidarZsWaterRm(isinf(lidarZsWaterRm)) = nan;
                lidarXYZ = [lidarXs, lidarYs, ...
                    double(lidarZsWaterRm)];

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
                plot(xYBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot3k(lidarXYZ);
                axis equal; view(2);

                % A preview of all the elevation data fetched.
                dispelev(rawElevData, 'mode', 'latlong');
                plot_google_map('MapType', 'satellite');

                % lidarZ - getLiDarZFromXYFct(lidarX, lidarY).
                figure; hold on;
                plot(xYBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), lidarZsWaterRm ...
                    - getLiDarZFromXYFct(lidarXs, lidarYs)]);
                axis equal; view(2);

                % lidarEles - getEleFromXYFct(lidarX, lidarY).
                figure; hold on;
                plot(xYBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), ...
                    lidarEles-getEleFromXYFct(lidarXs, lidarYs)]);
                axis equal; view(2);

                % lidarZ - lidarEles.
                figure; hold on;
                plot(xYBoundryPolygon, ...
                    'FaceColor','red','FaceAlpha',0.1);
                plot3k([lidarXYZ(:,1:2), ...
                    lidarZsWaterRm - lidarEles]);
                axis equal; view(2);
            end
            %====== END OF LIDAR DATA PROCESSING ======
        catch err
            disp('        There was an error!')
            dispErr(err);
            disp(['File which caused the error:', ...
                lidarFileRelDirs{idxF}])
            toc;
            error(...
                ['Error processing LiDAR data for ', ...
                'file # ', num2str(idxF), '/', ...
                num2str(numLidarFiles), '!']')
        end

    end

    curXYBoundryPolygons{idxPar} = xYBoundryPolygon;
    curLonLatBoundryPolygons{idxPar} = lonLatBoundryPolygon;

    toc;
    disp('        Done!');
end

% EOF