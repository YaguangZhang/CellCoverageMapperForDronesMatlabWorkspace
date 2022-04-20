%GETPOTENTIALGATEWAYLOCS Find existing high structures based on LiDAR data
%as potential locations for deploying gateways.
%
% A high structure needs to be of a height of at least
% MIN_CANDIDATE_HEIGHT_IN_M and at most MAX_CANDIDATE_HEIGHT_IN_M. This is
% decided by comparing the locally highest LiDAR point with both (i) the
% ground elevation and (ii) the minimum lowest LiDAR point. Both
% differences need to be big enough but not too big.
%
% The structure can not be a tree.
%
% We will inspect the data (roughly) in a circle with the radius
% LOCAL_CIRCLE_R_IN_M centered at the point of interest. The candidate
% stucture height test also (to some degree) make sure the structure is
% contained locally in this circle.
%
% MATLAB 2021b is required for vegetation classification.
%
% Yaguang Zhang, Purdue, 10/01/2021

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;
% Set this to true for extra plots.
%#ok<*UNRCH>
DEBUG = true;

%% Script Parameters

% The minimum and maximum height in meters for valid candidate structures.
% Reference values:
%   - the Purdue Water Tower
%       is 142 ft (43.28 m) tall and 92 ft (28.04) wide (the max radius is
%       about 22 m measured on Google Maps);
%   - cell towers
%       are typically 50 to 200 ft (15.24 to 60.96 m) high;
%   - the tallest building in Indiana
%       is the Salesforce Tower and its antenna spire is of a height of 811
%       ft (247 m);
%   - utility pole
%       can reach heights of 120 ft (36.58 m) or more to satisfy clearance
%       requirements.
MIN_CANDIDATE_HEIGHT_IN_M = 20;
MAX_CANDIDATE_HEIGHT_IN_M = 250;
LOCAL_CIRCLE_R_IN_M = 30;
MIN_NEAR_EDGE_R_IN_M = LOCAL_CIRCLE_R_IN_M*0.9;

% To avoid extremely narrow high peaks, which are likely to be caused by
% LiDAR data error, we will make sure the closest points are all high
% enough (>=MIN_HEIGHT_RATIO_OF_CLOSEST_HIGH_PT*structHInMOfCandLoc).
NUM_OF_CLOSEST_PT = 8;          % Excluding the candidate location.
MIN_HEIGHT_RATIO_OF_CLOSEST_HIGH_PT = 0.5;
MIN_NUM_OF_CLOSEST_HIGH_PT = 2; % Excluding the candidate location.

% For any candidate location, we will make sure there is no other candidate
% locations too nearby.
MIN_CANDIDATE_DIST_IN_M = 10;

% Change PRESET to run the code for different areas.
%   - 'ACRE_EXACT'
%     Purdue Agronomy Center for Research and Education (ACRE).
%   - 'ACRE_EXTENDED_2KM' / 'ACRE_EXTENDED_4KM' / 'ACRE_EXTENDED_8KM'
%     Extended area for ACRE_EXACT with 2 km radius (about 1.2 miles) / 4
%     km radius (about 2.5 miles) / 8 km radius (about 5 miles).
%   - 'Purdue'
%     A square area covering the West Lafayette Purdue campus.
PRESET = 'Purdue';

% The LiDAR data set to use. Currently we only suppor the 2019 Indiana
% state-wide digital surface model (DSM) data from:
%       https://lidar.jinha.org/
% Set this to "IN_DSM_2019"/"IN_DSM_2019_DEMO" for the complete/a demo data
% set. Please refer to the Preprocessing LiDAR Data section for a complete
% list of supported LiDAR data sets.
LIDAR_DATA_SET_TO_USE = 'IN_DSM_2019';

% The grid resolution to cover the area of interest. It should be
% comparable with the LiDAR data resolution. Note that the 2019 DSM data
% set has a resolution of 5 feet (1.524m).
GRID_RESOLUTION_IN_M = 1.5;

% The zone label to use in the UTM (x, y) system.
UTM_ZONE = '16 T';

% For figures.
scatterPtSize = 12;
markerThickLineWidth = 2;

%% Reuse History Grid If Possible

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', ['CandLocs_', PRESET, ...
    '_HRange', num2str(MIN_CANDIDATE_HEIGHT_IN_M), ...
    'm_to_', num2str(MAX_CANDIDATE_HEIGHT_IN_M), ...
    'm_LocalCircleR_', num2str(LOCAL_CIRCLE_R_IN_M), ...
    'm_MinCandDist_', num2str(MIN_CANDIDATE_DIST_IN_M), '_m']);
if ~exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% The absolute path for debug results.
pathToSaveDebugResults = fullfile(pathToSaveResults, 'DEBUG');
if ~exist(pathToSaveDebugResults, 'dir')~=7
    mkdir(pathToSaveDebugResults);
end

dirToSaveSimConfigs = fullfile(pathToSaveResults, 'simConfigs.mat');
if exist(dirToSaveSimConfigs, 'file')
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Loading history results ...'])
    histSimConfigs = load(dirToSaveSimConfigs);
    simConfigs = histSimConfigs.simConfigs;

    disp(['    [', datestr(now, datetimeFormat), '] Done!'])
end

%% Preprocessing LiDAR Data
% Note: the first time of this may take a long time, depending on the size
% of the LiDAR data set, but (1) it supports recovery from interruptions,
% and (2) once we have gone through all the data once, loading the
% information would be very fast.

% For GPS and UTM conversions.
[simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone] ...
    = genUtmConvertersForFixedZone(UTM_ZONE);

% Set the dir to find the LiDAR data set.
switch LIDAR_DATA_SET_TO_USE
    case 'IN_DSM_2019'
        dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar_2019', 'IN', 'DSM');
    otherwise
        error(['Unkown LiDAR data set ', LIDAR_DATA_SET_TO_USE, '!'])
end

% Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need to
% remove preprocessIndianaLidarDataSet from path after things are done.
addpath(fullfile(pwd, 'lib', 'lidar'));
[lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessLidarDataSetDsm(dirToLidarFiles, ...
    simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone);
rmpath(fullfile(pwd, 'lib', 'lidar'));
lidarFileAbsDirs = cellfun(@(d) ...
    [dirToLidarFiles, strrep(d, '\', filesep)], ...
    lidarFileRelDirs, 'UniformOutput', false);

%% Configurations

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Configuring the simulation for PRESET ', PRESET, ' ...'])

if strcmp(PRESET(1:4), 'ACRE')
    % Read in the ACRE boundary.
    [utmXYBoundaryOfAcre, ~] = extractBoundaryFromKmzFile( ...
        fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'Lidar', 'ACRE', 'AcreExactBoundaryRaw', 'ACRE.kmz'));
    % Convert the boundary to (lat, lon) for plotting.
    [gpsLatsBoundaryOfAcre, gpsLonsBoundaryOfAcre] ...
        = simConfigs.utm2deg_speZone( ...
        utmXYBoundaryOfAcre(:,1), ...
        utmXYBoundaryOfAcre(:,2));
    gpsLatLonBoundaryOfAcre ...
        = [gpsLatsBoundaryOfAcre, gpsLonsBoundaryOfAcre];
end

switch PRESET
    case 'ACRE_EXACT'
        UTM_X_Y_BOUNDARY_OF_INTEREST = utmXYBoundaryOfAcre;
    case {'ACRE_EXTENDED_2KM', 'ACRE_EXTENDED_4KM', 'ACRE_EXTENDED_8KM'}
        idxStartOfKM = strfind(PRESET, 'KM');
        padDistInKm = str2double(PRESET(15:(idxStartOfKM-1)));
        % Extend the boundary.
        UTM_X_Y_BOUNDARY_OF_INTEREST = polybuffer( ...
            polyshape(utmXYBoundaryOfAcre), padDistInKm*1000);
        UTM_X_Y_BOUNDARY_OF_INTEREST = ...
            UTM_X_Y_BOUNDARY_OF_INTEREST.Vertices;
    case 'Purdue'
        UTM_X_Y_BOUNDARY_OF_INTEREST ...
            = constructUtmRectanglePolyMat(...
            [40.467341, -87.015762; ...
            40.501484, -86.979905]);
    otherwise
        error(['Unsupported preset "', PRESET, '"!'])
end

% Make sure the UTM boundary is closed.
if ~all( UTM_X_Y_BOUNDARY_OF_INTEREST(1, :) ...
        == UTM_X_Y_BOUNDARY_OF_INTEREST(end, :) )
    UTM_X_Y_BOUNDARY_OF_INTEREST(end+1, :) ...
        = UTM_X_Y_BOUNDARY_OF_INTEREST(1, :);
end

% Boundary in GPS (lat, lon) system.
if ~exist('GPS_LAT_LON_BOUNDARY_OF_INTEREST', 'var')
    [gpsLatsBoundaryOfInterest, gpsLonsBoundaryOfInterest] ...
        = simConfigs.utm2deg_speZone( ...
        UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
    GPS_LAT_LON_BOUNDARY_OF_INTEREST ...
        = [gpsLatsBoundaryOfInterest, gpsLonsBoundaryOfInterest];
end

% Centroids for the LiDAR files in UTM.
lidarFileXYCentroids ...
    = extractCentroidsFrom2DPolyCell(lidarFileXYCoveragePolyshapes);

% The .mat copies for the LiDAR data. For the 2019 dataset, they are stored
% in a cache folder.
lidarMatFileAbsDirs = lidarFileAbsDirs;
for idxMatF = 1:length(lidarMatFileAbsDirs)
    [lidarMatFPath, lidarMatFName, ~] ...
        = fileparts(lidarMatFileAbsDirs{idxMatF});
    lidarMatFileAbsDirs{idxMatF} = fullfile(lidarMatFPath, '..', ...
        'MatlabCache', [lidarMatFName, '.mat']);
end

disp(['    [', datestr(now, datetimeFormat), '] Done!'])

%% Create a Grid for the Extended Area of Interest
% We will extend the area of interest by LOCAL_CIRCLE_R_IN_M to make sure
% locations on the edge are treated the same as those far away from the
% edge.

if ~exist(dirToSaveSimConfigs, 'file')
    disp(' ')
    disp(['    [', datestr(now, datetimeFormat), ...
        '] Creating a grid for the extended area of interest ...'])

    disp(['        [', datestr(now, datetimeFormat), ...
        '] Generating grid points ...'])
    utmXYBoundaryOfInterestExt = polybuffer( ...
        polyshape(UTM_X_Y_BOUNDARY_OF_INTEREST), LOCAL_CIRCLE_R_IN_M);
    utmXYBoundaryOfInterestExt = utmXYBoundaryOfInterestExt.Vertices;

    mapMinX = min(utmXYBoundaryOfInterestExt(:,1));
    mapMaxX = max(utmXYBoundaryOfInterestExt(:,1));
    mapMinY = min(utmXYBoundaryOfInterestExt(:,2));
    mapMaxY = max(utmXYBoundaryOfInterestExt(:,2));

    mapWidthInM = mapMaxX-mapMinX;
    mapHeightInM = mapMaxY-mapMinY;

    mapXLabels = constructAxisGrid( ...
        mean([mapMaxX, mapMinX]), ...
        floor((mapMaxX-mapMinX)./GRID_RESOLUTION_IN_M), ...
        GRID_RESOLUTION_IN_M);
    mapYLabels = constructAxisGrid( ...
        mean([mapMaxY, mapMinY]), ...
        floor((mapMaxY-mapMinY)./GRID_RESOLUTION_IN_M), ...
        GRID_RESOLUTION_IN_M);
    [mapXs, mapYs] = meshgrid(mapXLabels,mapYLabels);

    % Discard map grid points out of the area of interest.
    boolsMapGridPtsToKeep = inpoly2([mapXs(:), mapYs(:)], ...
        utmXYBoundaryOfInterestExt);

    simConfigs.mapGridXYPts = [mapXs(boolsMapGridPtsToKeep), ...
        mapYs(boolsMapGridPtsToKeep)];

    % Convert UTM (x, y) to (lat, lon).
    [mapGridLats, mapGridLons] = simConfigs.utm2deg_speZone( ...
        simConfigs.mapGridXYPts(:,1), simConfigs.mapGridXYPts(:,2));
    simConfigs.mapGridLatLonPts = [mapGridLats, mapGridLons];

    disp(['        [', datestr(now, datetimeFormat), ...
        '] Fetching ground elevation and LiDAR z ...'])

    % Fetch ground elevation and LiDAR z data for the grid points.
    [simConfigs.mapGridGroundEles, simConfigs.mapGridLidarZs, ...
        curEleForNanPts] ...
        = generateProfileSamps( ...
        simConfigs.mapGridXYPts, simConfigs.utm2deg_speZone, ...
        lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
        lidarMatFileAbsDirs, 'both');
    boolsNanMapGridGroundEles = isnan(simConfigs.mapGridGroundEles);
    simConfigs.mapGridGroundEles(boolsNanMapGridGroundEles) ...
        = curEleForNanPts(boolsNanMapGridGroundEles);

    disp(['        [', datestr(now, datetimeFormat), ...
        '] Saving results ...'])
    save(dirToSaveSimConfigs, 'simConfigs', '-v7.3');

    % Generate overview plots for the grid.
    hOverviewForGrid = figure; hold on;
    scatter3(simConfigs.mapGridLatLonPts(:, 2), ...
        simConfigs.mapGridLatLonPts(:, 1), ...
        simConfigs.mapGridGroundEles, ...
        scatterPtSize, simConfigs.mapGridGroundEles, 'filled');
    xlabel('Longitude'); ylabel('Latitude');
    hCb = colorbar; hCb.Label.String = 'Ground Elevation (m)';
    plot_google_map('MapType', 'hybrid');
    saveas(hOverviewForGrid, fullfile(pathToSaveResults, ...
        'overviewForGrid_eles.jpg'));
    close(hOverviewForGrid);

    hOverviewForGrid = figure; hold on;
    scatter3(simConfigs.mapGridLatLonPts(:, 2), ...
        simConfigs.mapGridLatLonPts(:, 1), ...
        simConfigs.mapGridLidarZs, ...
        scatterPtSize, simConfigs.mapGridLidarZs, 'filled');
    xlabel('Longitude'); ylabel('Latitude');
    hCb = colorbar; hCb.Label.String = 'LiDAR z (m)';
    plot_google_map('MapType', 'hybrid');
    saveas(hOverviewForGrid, fullfile(pathToSaveResults, ...
        'overviewForGrid_lidarZs.jpg'));
    close(hOverviewForGrid);

    hOverviewForGrid = figure; hold on;
    structHs = simConfigs.mapGridLidarZs - simConfigs.mapGridGroundEles;
    scatter3(simConfigs.mapGridLatLonPts(:, 2), ...
        simConfigs.mapGridLatLonPts(:, 1), ...
        structHs, ...
        scatterPtSize, structHs, 'filled');
    xlabel('Longitude'); ylabel('Latitude');
    hCb = colorbar; hCb.Label.String = 'Struct Height (m)';
    plot_google_map('MapType', 'hybrid');
    saveas(hOverviewForGrid, fullfile(pathToSaveResults, ...
        'overviewForGrid_structHs.jpg'));
    close(hOverviewForGrid);

    disp(['    [', datestr(now, datetimeFormat), '] Done!'])
end

%% Terrain Classification

% ptCloud = pointCloud([simConfigs.mapGridXYPts, simConfigs.mapGridLidarZs]);
%
% [groundPtsIdx,nonGroundPtCloud,groundPtCloud] = segmentGroundSMRF(ptCloud);
% % Visualize ground and non-ground points in green and magenta, respectively
% figure
% pcshowpair(nonGroundPtCloud,groundPtCloud)
% title('Classified Ground and Non-Ground Points')
%
% neighbors = 10;
% [normals,curvatures,neighInds] = helperExtractFeatures(nonGroundPtCloud, ...
%     neighbors);
%
% % Specify the normal threshold and curvature threshold
% normalThresh = 0.85;
% curveThresh = 0.02;
% % Classify the points into building and vegetation
% labels = helperClassify(normals,curvatures,neighInds, ...
%     normalThresh,curveThresh);
%
% maskData = grdTruthLabels>=2 & grdTruthLabels<=6;
%
% % Compress low, medium, and high vegetation to a single value
% grdTruthLabels(grdTruthLabels>=3 & grdTruthLabels<=5) = 4;
% % Update grdTruthLabels for metrics calculation
% grdTruthLabels(grdTruthLabels == 2) = 1;
% grdTruthLabels(grdTruthLabels == 4) = 2;
% grdTruthLabels(grdTruthLabels == 6) = 3;
%
% estimatedLabels = zeros(ptCloud.Count,1);
% estimatedLabels(groundPtsIdx) = 1;
% estimatedLabels(labels == 1) = 2;
% estimatedLabels(labels == 2) = 3;
%
% grdTruthLabels = grdTruthLabels(maskData);
% estimatedLabels = estimatedLabels(maskData);
%
% ptCloud = select(ptCloud,maskData);
% hFig = figure('Position',[0 0 900 400]);
% axMap1 = subplot(1,2,1,'Color','black','Parent',hFig);
% axMap1.Position = [0 0.2 0.5 0.55];
% pcshow(ptCloud.Location,grdTruthLabels,'Parent',axMap1)
% axis off
% title(axMap1,'Aerial Lidar Data with Ground Truth Labels')
% axMap2 = subplot(1,2,2,'Color','black','Parent',hFig);
% axMap2.Position = [0.5,0.2,0.5,0.55];
% pcshow(ptCloud.Location,estimatedLabels,'Parent',axMap2)
% axis off
% title(axMap2,'Aerial Lidar Data with Classified Labels')
%
% confusionMatrix = segmentationConfusionMatrix(estimatedLabels,double(grdTruthLabels));
% ssm = evaluateSemanticSegmentation({confusionMatrix}, ...
%     {'Ground' 'Vegetation' 'Building'},'Verbose',0);
% disp(ssm.DataSetMetrics)
% disp(ssm.ClassMetrics)

%% Search for Candidate Locs

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Searching for candidate locations ...'])

mapGridStructHInM = ...
    simConfigs.mapGridLidarZs - simConfigs.mapGridGroundEles;
numOfGridPts = length(mapGridStructHInM);

% Candidate points need to be (i) in the original area of interest and (ii)
% of heights (Lidar z - ground ele) in the desired range.
boolsGridPtIsCand = inpoly2(simConfigs.mapGridXYPts, ...
    UTM_X_Y_BOUNDARY_OF_INTEREST) ...
    & (mapGridStructHInM >= MIN_CANDIDATE_HEIGHT_IN_M) ...
    & (mapGridStructHInM <= MAX_CANDIDATE_HEIGHT_IN_M);

% For printing out progress info.
candCnt = 0;
progCnt = 0;
totalCands = sum(boolsGridPtIsCand);
numOfCandsProgress = ceil(totalCands/10);
for idxPt = find(boolsGridPtIsCand)'
    candCnt = candCnt+1;
    progCnt = progCnt+1;

    if progCnt>=numOfCandsProgress
        disp(['        [', datestr(now, datetimeFormat), ...
            '] ', num2str(candCnt/totalCands*100), ...
            '% (', num2str(candCnt), '/', num2str(totalCands), ') ...'])
        progCnt = 0;
    end

    % Now we need to make sure the candidate is a local maxima.
    curGridPtXY = simConfigs.mapGridXYPts(idxPt, :);
    curGridPtLatLon = simConfigs.mapGridLatLonPts(idxPt, :);
    curGridPtLidarZ = simConfigs.mapGridLidarZs(idxPt);

    % Fetch the points in the local circle. Note that rangesearch_fast is
    % way faster than rangesearch. One test we did costed rangesearch_fast
    % 0.023713 s, while rangesearch used 3.491118 s. The outputs are only a
    % little different near the edge. The results from rangesearch contains
    % more points there.
    %     indicesNearbyPts = rangesearch( ...
    %         simConfigs.mapGridXYPts, curGridPtXY, LOCAL_CIRCLE_R_IN_M);
    [indicesNearbyPts, distsNearbyPts] = rangesearch_fast( ...
        curGridPtXY, LOCAL_CIRCLE_R_IN_M, ...
        simConfigs.mapGridXYPts);

    latLonsNearbyPts = simConfigs.mapGridLatLonPts(indicesNearbyPts, :);
    lidarZsNearbyPts = simConfigs.mapGridLidarZs(indicesNearbyPts);
    structHsNearbyPts = mapGridStructHInM(indicesNearbyPts);

    % Make sure the current grid point of interest is the highest location
    % in the inspected circle.
    if curGridPtLidarZ<max(lidarZsNearbyPts)
        boolsGridPtIsCand(idxPt) = false;
        continue;
    end

    % Make sure the LiDAR z difference between the lowest point in the
    % circle and the grid point of interest is also within the expected
    % candidate height range.
    maxLidarZDiffInM = curGridPtLidarZ - min(lidarZsNearbyPts);
    if (maxLidarZDiffInM>MAX_CANDIDATE_HEIGHT_IN_M) ...
            || (maxLidarZDiffInM<MIN_CANDIDATE_HEIGHT_IN_M)
        boolsGridPtIsCand(idxPt) = false;
        continue;
    end

    % Make sure points near the edge of the circle are low enough compared
    % with the grid point of interest. Note that there is a bug in the
    % original rangesearch function (it is fixed in rangesearch_fast) where
    % the output distance is not correct.
    boolsNearEdgePts = distsNearbyPts>=MIN_NEAR_EDGE_R_IN_M;
    if min(curGridPtLidarZ-lidarZsNearbyPts(boolsNearEdgePts)) ...
            < MIN_CANDIDATE_HEIGHT_IN_M
        boolsGridPtIsCand(idxPt) = false;
        if DEBUG
            debugInfo = 'Rejected: High Edge Pts';
        else
            continue;
        end
    end

    % Make sure the nearest points are high enough to avoid very narrow
    % single peaks.
    [distsNearbyPtsSorted, indicesDistsNearbyPtsSorted] ...
        = sort(distsNearbyPts);
    % Note that we skip the first point because that should be the grid
    % point of interest itself with a distance of zero.
    indicesClosestPts = ...
        indicesDistsNearbyPtsSorted(2:(NUM_OF_CLOSEST_PT+1));
    if ~( ...
            sum( structHsNearbyPts(indicesClosestPts) ...
            >=( MIN_HEIGHT_RATIO_OF_CLOSEST_HIGH_PT ...
            *mapGridStructHInM(idxPt) ) ...
            ) ...
            >= MIN_NUM_OF_CLOSEST_HIGH_PT ...
            )
        boolsGridPtIsCand(idxPt) = false;
        if DEBUG
            if exist('debugInfo', 'var')
                debugInfo = [debugInfo, ' & Low Closest Pts']; %#ok<AGROW>
            else
                debugInfo = 'Rejected: Low Closest Pts';
            end
        else
            continue;
        end
    end

    hOverviewMap = figure; hold on;
    hClosestPts = plot( ...
        latLonsNearbyPts(indicesClosestPts, 2), ...
        latLonsNearbyPts(indicesClosestPts, 1), 'rs', ...
        'LineWidth', markerThickLineWidth);
    hEdgePts = plot( ...
        latLonsNearbyPts(boolsNearEdgePts, 2), ...
        latLonsNearbyPts(boolsNearEdgePts, 1), 'rd');
    hCurPt = plot(curGridPtLatLon(2), curGridPtLatLon(1), 'wo', ...
        'LineWidth', markerThickLineWidth);
    hSca3 = scatter3(latLonsNearbyPts(:,2), latLonsNearbyPts(:,1), ...
        lidarZsNearbyPts, ...
        scatterPtSize, lidarZsNearbyPts, 'filled');
    xlabel('Longitude'); ylabel('Latitude');
    hCb = colorbar; hCb.Label.String = 'LiDAR z (m)';
    legend([hCurPt, hClosestPts, hEdgePts], ...
        'Point of Interest', 'Closest Points', 'Edge Points', ...
        'AutoUpdate','off');
    hLegends = findobj(gcf, 'Type', 'Legend');
    for idxLegend = 1:length(hLegends)
        set(hLegends(idxLegend).BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.8]));
    end
    if exist('debugInfo', 'var')
        title(debugInfo);
        clearvars debugInfo;
    end
    plot_google_map('MapType', 'hybrid');
    % plot_google_map('MapType', 'hybrid', 'MapScale', true);

    if DEBUG
        saveas(hOverviewMap, fullfile(pathToSaveDebugResults, ...
            ['candGridPt_', num2str(idxPt), '.jpg']));
    end
    if boolsGridPtIsCand(idxPt)
        saveas(hOverviewMap, fullfile(pathToSaveResults, ...
            ['candGridPt_', num2str(idxPt), '.jpg']));
    end

    view(3);
    if DEBUG
        saveas(hOverviewMap, fullfile(pathToSaveDebugResults, ...
            ['candGridPt_', num2str(idxPt), '_3D.jpg']));
    end
    if boolsGridPtIsCand(idxPt)
        saveas(hOverviewMap, fullfile(pathToSaveResults, ...
            ['candGridPt_', num2str(idxPt), '_3D.jpg']));
    end

    view(2);
    if exist('gpsLatLonBoundaryOfAcre', 'var')
        plot(gpsLatLonBoundaryOfAcre(:,2), ...
            gpsLatLonBoundaryOfAcre(:,1), 'w:', ...
            'LineWidth', markerThickLineWidth);
    end
    plot(GPS_LAT_LON_BOUNDARY_OF_INTEREST(:,2), ...
        GPS_LAT_LON_BOUNDARY_OF_INTEREST(:,1), 'w-', ...
        'LineWidth', markerThickLineWidth);
    axis auto; plot_google_map('MapType', 'hybrid');
    if DEBUG
        saveas(hOverviewMap, fullfile(pathToSaveDebugResults, ...
            ['candGridPt_', num2str(idxPt), '_ZoomedOut.jpg']));
    end
    if boolsGridPtIsCand(idxPt)
        saveas(hOverviewMap, fullfile(pathToSaveResults, ...
            ['candGridPt_', num2str(idxPt), '_ZoomedOut.jpg']));
    end

    close(hOverviewMap);
end

disp(['    [', datestr(now, datetimeFormat), '] Done!'])

%% Overview Plot for the Candidate Locations

% EOF