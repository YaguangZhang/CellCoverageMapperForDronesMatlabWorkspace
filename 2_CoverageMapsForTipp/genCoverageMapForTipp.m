% GENCOVERAGEMAPFORTIPP Calculate path loss and generate coverage maps for
% Tippecanoe County, Indiana.
%
% If new simulation is needed, please delete previous result files to avoid
% progress recovery.
%
% Yaguang Zhang, Purdue, 06/19/2019

clear; clc; close all;

%% Configurations

dataTimeStrStart = datestr(datetime('now'));
timerValueStart = tic;

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
    'PostProcessingResults', '2_CoverageMapsForTipp');

% Configure other paths accordingly.
[~, TIPP_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_TIPP_LIDAR_TIF);
ABS_PATH_TO_SAVE_LIDAR = fullfile(pathToSaveResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_LiDAR.mat']);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(pathToSaveResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_Elevations.mat']);
ABS_PATH_TO_SAVE_COVERAGE_MAPS = fullfile(pathToSaveResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_CoverageMaps.mat']);
% We will save a cell, with each element being a boolean vectore indicating
% for a height of interest whether any of the cells have been calculated or
% not.
ABS_PATH_TO_SAVE_COMP_PROGRESS = fullfile(pathToSaveResults, ...
    [TIPP_LIDAR_LAS_FILENAME, '_ComputationProgress.mat']);

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
GRID_UNIT_LENGTH_IN_M = 250;

% Set this to be true to try to fetch elevation data from Google Maps.
% Currently, we do not have enough quota for this.
FLAG_FTECH_ELE_FROM_GOOGLE = false;

% Parameters for the extented Hata model.
TERRAIN_RES_IN_M = 25; % Terrain profile resolution.
CARRIER_FREQUENCY_IN_MHZ = 1900;
DEFAULT_TX_ANT_HEIGHT_IN_M = 50;

% Rx heights for different algorithms to inspect.
RX_ANT_HEIGHTS_IN_M_FOR_EHATA = [1.5; 10; 30; 50; 70; 90; 100];

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
EXPECTED_PL_RANGE_IN_DB = [35, 200];

%% Before Processing the Data

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

flagProcessInterrupted = false;

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
    
    % Convert survery feet to meter.
    lidarDataImg = distdim(lidarDataImg, 'ft', 'm');
    
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
        
        if FLAG_GEN_FIGS
            % Plot the interperlated elevation for the LiDAR data.
            hInterpElesForLiDAR = figure;
            plot3k([lidarLons lidarLats lidarEles], 'ColorBar', false);
            xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
            view(2); axis equal; axis tight;
            
            pathToSaveFig = fullfile(pathToSaveResults, ...
                'ElevationsInterpolated.jpg');
            saveas(hInterpElesForLiDAR, pathToSaveFig);
            
            % Test getEleFromXYFct with the LiDAR data locations.
            lidarElesFromUtm ...
                = getEleFromXYFct(lidarXYZ(:,1), lidarXYZ(:,2));
            
            % Plot the interperlated elevation for the LiDAR data.
            hInterpElesFromUtmForLiDAR = figure;
            plot3k([lidarLons lidarLats lidarElesFromUtm], 'ColorBar', false);
            xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
            view(2); axis equal; axis tight;
            
            pathToSaveFig = fullfile(pathToSaveResults, ...
                'ElevationsInterpolatedFromUtmCoors.jpg');
            saveas(hInterpElesFromUtmForLiDAR, pathToSaveFig);
            
            % Plot the difference between the two interperlated elevation
            % sets.
            hInterpElesDiffForLiDAR = figure;
            plot3k([lidarLons lidarLats lidarEles-lidarElesFromUtm], ...
                'Labels', {'', 'Longitude (degrees)', 'Latitude (degrees)', ...
                '', 'lidarEles-lidarElesFromUtm'});
            view(2); axis equal; axis tight;
            
            pathToSaveFig = fullfile(pathToSaveResults, ...
                'ElevationsInterpolatedDiff.jpg');
            saveas(hInterpElesDiffForLiDAR, pathToSaveFig);
        end
    end
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_ELEVATIONS, ...
        'lidarEles', 'getEleFromXYFct', ...
        'FLAG_FTECH_ELE_FROM_GOOGLE', '-append');
end

disp('    Done!')

%% Generate a Few Plots for the LiDAR Data

if FLAG_GEN_FIGS
    disp(' ')
    disp('    Generating plots for LiDAR Data ...')
    
    hZ = figure;
    plot3k(lidarXYZ, 'Labels', ...
        {'', 'x (m)', 'y (m)', '', 'z (m)'});
    grid on; view(2); axis equal; axis tight;
    
    pathToSaveFig = fullfile(pathToSaveResults, 'lidarZ.png');
    saveas(hZ, pathToSaveFig);
    
    downSampFactor = 10;
    hZDownSamped = figure;
    plot3k(lidarXYZ(1:downSampFactor:end, :), 'Labels', ...
        {['LiDAR z (Down sampled by ', num2str(downSampFactor) ')'], ...
        'x (m)', 'y (m)', '', 'z (m)'});
    grid on; view(2); axis equal; axis tight;
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        ['lidarZDownSamped_', num2str(downSampFactor), '.png']);
    saveas(hZDownSamped, pathToSaveFig);
    
    hDiffBetweenLidarZAndEle = figure;
    plot3k([lidarXYZ(:,1:2) lidarXYZ(:,3)-lidarEles], 'Labels', ...
        {'LiDAR z minus Elevation', ...
        'x (m)', 'y (m)', '', 'z-ele (m)'});
    grid on; view(2); axis equal; axis tight
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        'diffBetweenLidarZAndEle.png');
    saveas(hDiffBetweenLidarZAndEle, pathToSaveFig);
    
    disp('    Done!')
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

%% Area of Interest

disp(' ')
disp('    Locating area of interest ...')
% Determine the area of interest in the UTM system.
[rangeOfInterestXs, rangeOfInterestYs, rangeOfInterestZones] ...
    = deg2utm(LAT_LON_RANGES_OF_INTEREST(1:2), ...
    LAT_LON_RANGES_OF_INTEREST(3:4));
assert(all(strcmp( ...
    mat2cell(rangeOfInterestZones, ...
    ones(1, 2), 4), UTM_ZONE)), ...
    ['Locations in the area of interest are not all', ...
    ' in the expected UTM zone!']);

if FLAG_GEN_FIGS
    % Show the antennas and the area of interest on the Lidar data.
    hCellAntsOnZ = figure;
    plot3k(lidarXYZ, 'ColorBar', false);
    hold on;
    hAreaOfInt = plot3([rangeOfInterestXs(1); rangeOfInterestXs(1); ...
        rangeOfInterestXs(2); rangeOfInterestXs(2); rangeOfInterestXs(1)], ...
        [rangeOfInterestYs(1); rangeOfInterestYs(2); ...
        rangeOfInterestYs(2); rangeOfInterestYs(1); rangeOfInterestYs(1)], ...
        ones(5,1).*max(lidarXYZ(:,3)), ...
        '--', 'LineWidth', 1.5, 'Color', ones(1,3).*0.7);
    hCellAnt = plot3(cellAntsXs, cellAntsYs, ...
        cellAntsLatLonXYAlt(:,5), '*r');
    grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    view(2); axis equal; axis tight;
    legend([hAreaOfInt, hCellAnt], 'Area of Interest', 'Cellular Antennas');
    
    pathToSaveFig = fullfile(pathToSaveResults, ...
        'areaOfInterestAndCellAntsOnLidarZ.png');
    saveas(hCellAntsOnZ, pathToSaveFig);
end

% Create a grid for the area of interest.
areaOfInterestXs ...
    = rangeOfInterestXs(1):GRID_UNIT_LENGTH_IN_M...
    :rangeOfInterestXs(2);
areaOfInterestYs ...
    = rangeOfInterestYs(1):GRID_UNIT_LENGTH_IN_M...
    :rangeOfInterestYs(2);

disp('    Done!')

%% Generate Path Loss Maps

genPathLossMaps;

%% Combine Path Loss Maps

RX_ANT_HEIGHTS_TO_INSPECT_IN_M = RX_ANT_HEIGHTS_IN_M_FOR_EHATA;
combinePathLossMaps;

%% Clear Things Up if Necessary
% This should be put at the very end of a script using the cpp ehata lib.
if libisloaded('ehata')
    unloadlibrary('ehata');
end

% Print out a ruler to indicate everything is done.
disp(fileNameHintRuler);

toc(timerValueStart);

%% Log Execution Time

logExecTimeForCovMapGen;

% EOF