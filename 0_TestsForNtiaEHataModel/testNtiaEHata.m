% TESTNTIAEHATA Compute the path loss via NTIA eHata for several Tx & Rx
% pairs.
%
% We compare the results with those from NTIA to make sure the NTIA library
% is set up properly.
%
% Yaguang Zhang, Purdue, 06/10/2019

clear; clc; close all;

%% Configurations

% Locate the current working directory.
cd(fileparts(mfilename('fullpath')));
cd('..'); addpath('lib');
curFileName = mfilename;
fileNameHintRuler = hintScriptName(curFileName);

% 'Matlab' or 'CPlusPlus'.
LIBRARY_TO_USE = 'CPlusPlus';

% Load the NTIA eHata library first, if necessary, to avoid the "unable to
% find ehata" error.
if strcmpi(LIBRARY_TO_USE, 'cplusplus')
    loadEHataCppLib;
end

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
setPath;

% Set this to be true to generate figures.
FLAG_GEN_FIGS = true;

% The absolute path to the Lidar .las file.
ABS_PATH_TO_TIPP_LIDAR_TIF = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', 'Tipp', 'output_be.tif');

% The absolute path to the antenna infomation file.
ABS_PATH_TO_TIPP_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'RandomizedCarrierSitesv2.csv');

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '0_TestsForNtiaEHataModel');

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_EVALUATION_MAT = fullfile(pathToSaveResults, ...
    'ntiaEHataEvaluationLog.mat');
ABS_PATH_TO_SAVE_EVALUATION_LOG = fullfile(pathToSaveResults, ...
    'ntiaEHataEvaluationLog.csv');

% Paths for data reuse.
pathToFetchGeoResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '2_CoverageMapsForTipp');
[~, TIPP_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_TIPP_LIDAR_TIF);
ABS_PATH_TO_SAVE_LIDAR = fullfile(pathToFetchGeoResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_LiDAR.mat']);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(pathToFetchGeoResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_Elevations.mat']);

% The UTM zone expected for the Tippecanoe data.
UTM_ZONE = '16 T';

% The State plane code for the Tippecanoe data.
STATE_PLANE_CODE_TIPP = 'indiana west';

% The raster LiDAR data we have has a spacial resolution of 5 survey feet
% (~1.524 m). We will down sample this to speed up the computation.
LIDAR_DOWN_SAMPLE_FACTOR_PER_SIDE = 80; % New resolution: ~121.92 m.

% The area of interest in terms of (minLat, maxLat, minLon, maxLon) for
% generating the coverage maps. This is a little smaller than the whole
% county because we do not have all the LiDAR data to cover that large of
% an area.
LAT_LON_RANGES_OF_INTEREST ...
    = [40.242314, 40.536134, -87.083795, -86.700450];

% Set this to be true to try to fetch elevation data from Google Maps.
% Currently, we do not have enough quota for this.
FLAG_FTECH_ELE_FROM_GOOGLE = false;

% Parameters for the extented Hata model.
TERRAIN_RES_IN_M = 50; % Terrain profile resolution.
CARRIER_FREQUENCY_IN_MHZ = 1900;
DEFAULT_TX_ANT_HEIGHT_IN_M = 200;

% Rx heights for different algorithms to inspect.
RX_ANT_HEIGHTS_IN_M_FOR_EHATA = [1.5; 10];

switch lower(LIBRARY_TO_USE)
    case 'cplusplus'
        % The quantile percent not exceeded of the signal.
        %   Limits: 0 < reliability < 1
        NTIA_EHATA_RELIABILITY = 0.95;
        % The NLCD environment code.
        NTIA_EHATA_ENVIRO_CODE = 82; % Cultivated Crops.
    case 'matlab'
        REGION = 'Suburban';
    otherwise
        error(['Unsupported library: ', model, '!'])
end

% For coordinating the color used in plot3k.
EXPECTED_PL_RANGE_IN_DB = [35, 150];

% We will test a few manually chosen locations.
numOfTestSets = 2;
indicesCellAntsToTest = [1, 2];
rxLatLonsToTest =[40.436052, -86.962681; ...
    40.470587, -86.986578; ...
    40.445117, -87.032908; ...
    40.533119, -86.918199; ...
    40.377548, -87.053936; ...
    40.256138, -87.064580; ...
    40.260070, -86.713350];

%% Before Processing the Data

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

%% Load the Lidar Data

disp(' ')
disp('    Loading LiDAR data ...')

if exist(ABS_PATH_TO_SAVE_LIDAR, 'file')==2
    disp('        Loading previous results ...')
    load(ABS_PATH_TO_SAVE_LIDAR);
else
    disp(' ')
    disp('        Processing the raw LiDAR data set ...')
    
    [lidarDataImg, ~] ...
        = geotiffread(ABS_PATH_TO_TIPP_LIDAR_TIF);
    [lidarRasterXLabels, lidarRasterYLabels] ...
        = pixcenters(geotiffinfo(ABS_PATH_TO_TIPP_LIDAR_TIF));
    
    % Down sample the LiDAR data.
    lidarDataImg ...
        = lidarDataImg(1:LIDAR_DOWN_SAMPLE_FACTOR_PER_SIDE:end, ...
        1:LIDAR_DOWN_SAMPLE_FACTOR_PER_SIDE:end);
    lidarRasterXLabels ...
        = lidarRasterXLabels(1:LIDAR_DOWN_SAMPLE_FACTOR_PER_SIDE:end);
    lidarRasterYLabels ...
        = lidarRasterYLabels(1:LIDAR_DOWN_SAMPLE_FACTOR_PER_SIDE:end);
    
    % Prepare the results for building a grid.
    [ lidarRasterXLabels, lidarRasterYLabels, lidarDataImg ] ...
        = sortGridMatrixByXY(lidarRasterXLabels', lidarRasterYLabels', ...
        lidarDataImg);
    
    % Convert raster (row, col) to (lat, lon).
    [lidarRasterXs, lidarRasterYs] ...
        = meshgrid(lidarRasterXLabels, lidarRasterYLabels);
    [lidarLons, lidarLats] ...
        = sp_proj(STATE_PLANE_CODE_TIPP, 'inverse', ...
        lidarRasterXs(:), lidarRasterYs(:), 'sf');
    
    % Store the new (x,y,z) data.
    lidarLats = lidarLats(:);
    lidarLons = lidarLons(:);
    [lidarXs, lidarYs, lidarRasterZones] ...
        = deg2utm(lidarLats, lidarLons);
    assert(all(strcmp( ...
        mat2cell(lidarRasterZones, ...
        ones(1, length(lidarRasterZones)), 4), ...
        UTM_ZONE)), ...
        'LiDAR raster locations are not all in the expected UTM zone!');
    lidarXYZ = [lidarXs, lidarYs, lidarDataImg(:)];
    
    % Create a function to get LiDAR z from UTM coordinates.
    fctLonLatToLidarStatePlaneXY ...
        = @(lon, lat) sp_proj(STATE_PLANE_CODE_TIPP, ...
        'forward', lon, lat, 'sf');
    getLiDarZFromStatePlaneXYFct = @(spXs, spYs) ...
        interp2(lidarRasterXs, lidarRasterYs, ...
        lidarDataImg, spXs, spYs);
    getLiDarZFromXYFct ...
        = @(xs, ys) genRasterLidarZGetter( ...
        getLiDarZFromStatePlaneXYFct, fctLonLatToLidarStatePlaneXY, ...
        xs, ys, UTM_ZONE);
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_LIDAR, ...
        'lidarXYZ', 'lidarLats', 'lidarLons', 'getLiDarZFromXYFct');
end

[lidarNumSamps, ~] = size(lidarXYZ);

disp('    Done!')

%% Load the (lat, lon, ele) Information

disp(' ')
disp('    Generating elevation information ...')

if exist(ABS_PATH_TO_SAVE_ELEVATIONS, 'file')==2
    disp('        Loading previous results ...')
    load(ABS_PATH_TO_SAVE_ELEVATIONS);
else
    disp('        Fetching elevation data ...')
    
    latRange = [min(lidarLats) max(lidarLats)];
    lonRange = [min(lidarLons) max(lidarLons)];
    % Use USGS 1/3 arc-second (~10m) resolution data for US terrain
    % elevation.
    region = fetchregion(latRange, lonRange, 'display', true);
    rawElevData = region.readelevation(latRange, lonRange, ...
        'sampleFactor', 1, ...
        'display', true);
    
    % Save a preview of the elevation data fetched.
    dispelev(rawElevData, 'mode', 'latlong');
    hElePreview = gcf;
    axis tight;
    
    pathToSaveFig = fullfile(pathToSaveResults, 'ElevationsRaw.jpg');
    saveas(hElePreview, pathToSaveFig);
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, 'rawElevData', ...
        'lidarLats', 'lidarLons');
end

% we will only try to generate the elevation information when it is not
% available.
if (~exist('lidarEles', 'var') || all(isnan(lidarEles)))
    if FLAG_FTECH_ELE_FROM_GOOGLE
        lidarEles = nan;
        % Use Google as the elevation source for better resolution.
        % However, this approach is limited by the quota Google provides.
        disp('        Fetching elevation information from Google Maps ...')
        % Fetch the elevation data from Google Maps.
        GOOGLE_MAPS_API = 'AIzaSyDlkaE_QJxvRJpTutWG0N-LCvoT0e7FPHE';
        
        countTrials = 0;
        while isnan(lidarEles)
            countTrials = countTrials+1;
            disp(['            Trial # ', num2str(countTrials), ' ...']);
            lidarEles = getElevationsFromGoogle(lidarLats, lidarLons, ...
                GOOGLE_MAPS_API);
        end
    else
        disp('        Fitting elevation data in the (lon, lat) system ...')
        % Order the raw elevation data so that both lat and lon are
        % monotonically increasing.
        [rawElevDataLonsSorted, rawElevDataLatsSorted, ...
            rawElevDataElesSorted] ...
            = sortGridMatrixByXY(rawElevData.longs, rawElevData.lats, ...
            rawElevData.elev);
        
        % Create a grid for the elevation data.
        [tippElevDataLons, tippElevDataLats] = meshgrid( ...
            rawElevDataLonsSorted, rawElevDataLatsSorted);
        
        % Interperlate the data with lat and lon.
        lidarEles = interp2(tippElevDataLons, tippElevDataLats, ...
            rawElevDataElesSorted, lidarLons, lidarLats);
        
        % Create a function to get elevation from UTM coordinates in the
        % same way.
        getEleFromXYFct = @(xs, ys) ...
            genUtmEleGetter(tippElevDataLats, tippElevDataLons, ...
            rawElevDataElesSorted, xs, ys, UTM_ZONE);
    end
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, ...
        'lidarEles', 'getEleFromXYFct', ...
        'FLAG_FTECH_ELE_FROM_GOOGLE', '-append');
end

disp('    Done!')

%% Load Antenna Info

disp(' ')
disp('    Loading cellular antenna information ...')

cellAntsLatLon = csvread(ABS_PATH_TO_TIPP_CELL_ANTENNAS_CSV, 1, 1);

% Keep only the cell towers we need.
boolSCellAntsToKeep ...
    = cellAntsLatLon(:,1)>=LAT_LON_RANGES_OF_INTEREST(1) ...
    & cellAntsLatLon(:,1)<=LAT_LON_RANGES_OF_INTEREST(2) ...
    & cellAntsLatLon(:,2)>=LAT_LON_RANGES_OF_INTEREST(3) ...
    & cellAntsLatLon(:,2)<=LAT_LON_RANGES_OF_INTEREST(4);
cellAntsLatLon = cellAntsLatLon(boolSCellAntsToKeep, :);

[numOfCellAnts, ~] = size(cellAntsLatLon);

% Add UTM x and y.
[cellAntsXs, cellAntsYs, cellAntsZones] ...
    = deg2utm(cellAntsLatLon(:,1), cellAntsLatLon(:,2));
assert(all(strcmp( ...
    mat2cell(cellAntsZones, ones(1, numOfCellAnts), 4), ...
    UTM_ZONE)), ...
    'Cell antenna locations are not all in the expected UTM zone!');

% We will use the default cell tower antenna height to compute their
% altitude.
cellAntsAlts = getEleFromXYFct(cellAntsXs, cellAntsYs) ...
    + DEFAULT_TX_ANT_HEIGHT_IN_M;

cellAntsLatLonXYAlt = [cellAntsLatLon ...
    cellAntsXs cellAntsYs cellAntsAlts];

disp('    Done!')

%% Evaluate Path Loss

cellAntsLatLonToTest = cellAntsLatLon(indicesCellAntsToTest, :);
[numRxLocsToTest, ~] = size(rxLatLonsToTest);
pathLosses = cell(numOfTestSets);
% Prepare a matrix for outputting the results as a csv file.
csvTitle = {'txLat', 'txLon', 'txAntennaHeightM', ...
    'rxLat', 'rxLon', 'rxAntennaHeightM', ...
    'txLocToRxLocDistM', 'pathLossDb'};
csvResultsAsMat = nan(numOfTestSets.*numRxLocsToTest, length(csvTitle));
csvRowCnt = 1;
for idxTestSet = 1:numOfTestSets
    idxCellAntenna = indicesCellAntsToTest(idxTestSet);
    curTxLatLon = cellAntsLatLonToTest(idxTestSet, :);
    curRxAntH = RX_ANT_HEIGHTS_IN_M_FOR_EHATA(idxTestSet);
    
    pathLosses{idxTestSet} = nan(numRxLocsToTest, 1);
    for idxRx = 1:numRxLocsToTest
        curRxLatLon = rxLatLonsToTest(idxRx, :);
        
        % Parameters needed by the NTIA eHata library.
        curBaseAntXY = cellAntsLatLonXYAlt(idxCellAntenna, 3:4);
        % Get the altitude (ground elevation + antenna height) of the base
        % station antenna.
        curBaseAntAltInM = cellAntsLatLonXYAlt(idxCellAntenna, 5);
        % Get the ground elevation for the base station antenna.
        curBaseAntEleInM ...
            = getEleFromXYFct(curBaseAntXY(1), curBaseAntXY(2));
        
        curBaseAntHeightInM = curBaseAntAltInM - curBaseAntEleInM;
        
        [areaOfInterestXs, areaOfInterestYs, curRxZone] ...
            = deg2utm(curRxLatLon(1), curRxLatLon(2));
        assert(strcmp(curRxZone, UTM_ZONE), ...
            'RX location is not in the expected UTM zone!');
        
        [pathLosses{idxTestSet}(idxRx), ~, ~] ...
            = genMedianBasicPropLossMaps( ...
            CARRIER_FREQUENCY_IN_MHZ, ...
            curBaseAntXY, curBaseAntHeightInM, ...
            areaOfInterestXs, areaOfInterestYs, curRxAntH, ...
            getEleFromXYFct, getLiDarZFromXYFct, ...
            NTIA_EHATA_ENVIRO_CODE, ...
            TERRAIN_RES_IN_M, LIBRARY_TO_USE, ...
            NTIA_EHATA_RELIABILITY, 0);
        
        csvResultsAsMat(csvRowCnt, :) ...
            = [curTxLatLon(1), curTxLatLon(2), curBaseAntHeightInM, ...
            curRxLatLon(1), curRxLatLon(2), curRxAntH, ...
            norm(curBaseAntXY-[areaOfInterestXs areaOfInterestYs]), ...
            pathLosses{idxTestSet}(idxRx)];
        
        csvRowCnt = csvRowCnt+1;
    end
end

hCsvFile = fopen(ABS_PATH_TO_SAVE_EVALUATION_LOG,'w');
fprintf(hCsvFile,'%s,', csvTitle{1:end-1});
fprintf(hCsvFile,'%s\n', csvTitle{end});
fclose(hCsvFile);
dlmwrite(ABS_PATH_TO_SAVE_EVALUATION_LOG, csvResultsAsMat, ...
    'delimiter', ',', 'precision', '%.6f', '-append');

save(ABS_PATH_TO_SAVE_EVALUATION_MAT, 'CARRIER_FREQUENCY_IN_MHZ', ...
    'NTIA_EHATA_RELIABILITY', 'NTIA_EHATA_ENVIRO_CODE', ...
    'TERRAIN_RES_IN_M', 'csvTitle', 'csvResultsAsMat');

%% Plot Results on a Map

if FLAG_GEN_FIGS
    for idxTestSet = 1:numOfTestSets
        hCurPLMap = figure; hold on;
        
        hAreaOfInt = plot3([LAT_LON_RANGES_OF_INTEREST(3); ...
            LAT_LON_RANGES_OF_INTEREST(3); ...
        LAT_LON_RANGES_OF_INTEREST(4); LAT_LON_RANGES_OF_INTEREST(4); ...
        LAT_LON_RANGES_OF_INTEREST(3)], ...
        [LAT_LON_RANGES_OF_INTEREST(1); LAT_LON_RANGES_OF_INTEREST(2); ...
        LAT_LON_RANGES_OF_INTEREST(2); LAT_LON_RANGES_OF_INTEREST(1); ...
        LAT_LON_RANGES_OF_INTEREST(1)], ...
        ones(5,1).*max(lidarXYZ(:,3)), ...
        '--', 'LineWidth', 1.5, 'Color', ones(1,3).*0.7);
    
        hTxs = plot(cellAntsLatLonToTest(idxTestSet,2), ...
            cellAntsLatLonToTest(idxTestSet,1), 'rx', ...
            'MarkerSize', 12, 'LineWidth', 3);
        boolsNonInfPL = ~isinf(pathLosses{idxTestSet});
        plot3k([rxLatLonsToTest(boolsNonInfPL,2), ...
            rxLatLonsToTest(boolsNonInfPL,1), ...
            pathLosses{idxTestSet}(boolsNonInfPL)], ...
            'Marker', {'.', 30}, ...
            'Labels', {'', 'Longitude (degrees)', 'Latitude (degrees)', ...
                '', 'Path Loss (dB)'});
        grid on; view(2); legend(hTxs, 'TX');
        plotGoogleMapAfterPlot3k(gcf, 'satellite');
        
        pathToSaveFig = fullfile(pathToSaveResults, ...
            ['testSet', num2str(idxTestSet), 'Results.png']);
        saveas(hCurPLMap, pathToSaveFig);
    end
end

% EOF