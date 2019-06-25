% GENCOVERAGEMAPFORACRE Calculate path loss and generate coverage maps for
% Acre.
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
ABS_PATH_TO_ACRE_LIDAR_LAS = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', 'ACRE', 'points.las');

% The absolute path to the antenna infomation file.
ABS_PATH_TO_ACRE_CELL_ANTENNAS_CSV = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'ACRE', 'CellAntennasLatLonAlt.csv');

% The absolute path to save results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '1_CoverageMapsForAcre');

% Configure other paths accordingly.
[~, ACRE_LIDAR_LAS_FILENAME, ~] = fileparts(ABS_PATH_TO_ACRE_LIDAR_LAS);
ABS_PATH_TO_SAVE_LIDAR = fullfile(pathToSaveResults, ...
    [ACRE_LIDAR_LAS_FILENAME, '_LiDAR.mat']);
ABS_PATH_TO_SAVE_ELEVATIONS = fullfile(pathToSaveResults, ...
    [ACRE_LIDAR_LAS_FILENAME, '_Elevations.mat']);

% The UTM zone expected for the ACRE data.
UTM_ZONE = '16 T';

% The area of interest in terms of (minLat, maxLat, minLon, maxLon) for
% generating the coverage maps.
LAT_LON_RANGES_OF_INTEREST ...
    = [40.467341, 40.501484, -87.015762, -86.979905];
GRID_UNIT_LENGTH_IN_M = 25;

% Set this to be true to try to fetch elevation data from Google Maps.
% Currently, we do not have enough quota for this.
FLAG_FTECH_ELE_FROM_GOOGLE = false;

% Parameters for the extented Hata model.
TERRAIN_RES_IN_M = 10; % Terrain profile resolution.
CARRIER_FREQUENCY_IN_MHZ = 1900;

% Rx heights for different algorithms to inspect.
RX_ANT_HEIGHTS_IN_M_FOR_EHATA = [1.5; 10; 50; 100];

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
    
    lidarData = lasdata(ABS_PATH_TO_ACRE_LIDAR_LAS);
    lidarXYZ = [lidarData.x, lidarData.y, lidarData.z];
    
    % The stored intensity is unit16 and we need double for plotting.
    lidarIntensity = double(lidarData.get_intensity);
    
    % Create a function to get LiDAR z from UTM coordinates.
    getLiDarZFromXYFct ...
        = scatteredInterpolant(lidarData.x, lidarData.y, lidarData.z);
    
    % Save the results.
    save(ABS_PATH_TO_SAVE_LIDAR, ...
        'lidarXYZ', 'lidarIntensity', ...
        'getLiDarZFromXYFct');
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
    disp('        Converting UTM to (lat, lon) ...')
    % UTM to (lat, lon).
    [lidarLats, lidarLons] = utm2deg(lidarXYZ(:,1), lidarXYZ(:,2), ...
        repmat(UTM_ZONE, lidarNumSamps, 1));
    
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
        [acreElevDataLons, acreElevDataLats] = meshgrid( ...
            rawElevDataLonsSorted, rawElevDataLatsSorted);
        
        % Interperlate the data with lat and lon.
        lidarEles = interp2(acreElevDataLons, acreElevDataLats, ...
            rawElevDataElesSorted, lidarLons, lidarLats);
        
        % Create a function to get elevation from UTM coordinates in the
        % same way.
        getEleFromXYFct = @(xs, ys) ...
            genUtmEleGetter(acreElevDataLats, acreElevDataLons, ...
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
    
    hIntensity = figure;
    plot3k([lidarXYZ(:,1), lidarXYZ(:,2), lidarIntensity], 'Labels', ...
        {'', 'x (m)', 'y (m)', '', 'LiDAR Intenstiy'});
    grid on; view(2); axis equal; axis tight;
    
    pathToSaveFig = fullfile(pathToSaveResults, 'lidarIntensity.png');
    saveas(hIntensity, pathToSaveFig);
    
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

cellAntsLatLonAlt = csvread(ABS_PATH_TO_ACRE_CELL_ANTENNAS_CSV);
[numOfCellAnts, ~] = size(cellAntsLatLonAlt);

% Add UTM x and y.
[cellAntsXs,cellAntsYs,cellAntsZones] ...
    = deg2utm(cellAntsLatLonAlt(:,1), cellAntsLatLonAlt(:,2));
assert(all(strcmp( ...
    mat2cell(cellAntsZones, ones(1, numOfCellAnts), 4), ...
    UTM_ZONE)), ...
    'Cell antenna locations are not all in the expected UTM zone!');
cellAntsLatLonXYAlt = [cellAntsLatLonAlt(:,1:2) ...
    cellAntsXs cellAntsYs cellAntsLatLonAlt(:,3)];

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
    hCellAnt = plot3(cellAntsXs, cellAntsYs, cellAntsLatLonAlt(:,3), '*r');
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

disp(' ')
disp('    Generating coverage maps ...')

% Generate a propogation map for each height to inspect.
numOfHs = length(RX_ANT_HEIGHTS_IN_M_FOR_EHATA);
[coverageMapsEHata, coverageMapsEHataXLabels, coverageMapsEHataYLabels, ...
    towerPathLossMapsEHata, towerPathLossMapsEHataXLabels, ...
    towerPathLossMapsEHataYLabels] ...
    = deal(cell(numOfHs, 1));

disp(' ')
disp('        Computing path loss via the extended Hata model ...')
for idxH = 1:numOfHs
    disp(' ')
    disp('            Closing figures to save RAM ...')
    close all;
    disp('            Done!')
    disp(' ')
    
    disp(['            Height #', num2str(idxH), '/', num2str(numOfHs)]);
    curRxAntH = RX_ANT_HEIGHTS_IN_M_FOR_EHATA(idxH);
    % Generate a propogation map for each cell tower. We will first
    % generate the path loss maps, one for each cell tower location.
    [pathLossMaps, pathLossMapXLabels, pathLossMapYLabels] ...
        = deal(cell(numOfCellAnts, 1));
    for idxCellAntenna = 1:numOfCellAnts
        disp(['                Cell antenna #', ...
            num2str(idxCellAntenna), '/', num2str(numOfCellAnts)]);
        curBaseAntXY = cellAntsLatLonXYAlt(idxCellAntenna, 3:4);
        % Get the altitude (ground elevation + antenna height) of the base
        % station antenna.
        curBaseAntAltInM = cellAntsLatLonXYAlt(idxCellAntenna, 5);
        % Get the ground elevation for the base station antenna.
        curBaseAntEleInM ...
            = getEleFromXYFct(curBaseAntXY(1), curBaseAntXY(2));
        
        curBaseAntHeightInM = curBaseAntAltInM - curBaseAntEleInM;
        
        switch lower(LIBRARY_TO_USE)
            case 'cplusplus'
                [pathLossMaps{idxCellAntenna}, ...
                    pathLossMapXLabels{idxCellAntenna}, ...
                    pathLossMapYLabels{idxCellAntenna}] ...
                    = genMedianBasicPropLossMaps( ...
                    CARRIER_FREQUENCY_IN_MHZ, ...
                    curBaseAntXY, curBaseAntHeightInM, ...
                    areaOfInterestXs, areaOfInterestYs, curRxAntH, ...
                    getEleFromXYFct, getLiDarZFromXYFct, ...
                    NTIA_EHATA_ENVIRO_CODE, ...
                    TERRAIN_RES_IN_M, LIBRARY_TO_USE, ...
                    NTIA_EHATA_RELIABILITY);
                curFigFileName = [ ...
                    'CellTowerPathLoss_RxHeight_', num2str(curRxAntH), ...
                    '_TxCell_', num2str(idxCellAntenna), ...
                    '_eHataLib_', LIBRARY_TO_USE, ...
                    '_TerrainCode_', num2str(NTIA_EHATA_ENVIRO_CODE), ...
                    '_Reliability_', num2str(NTIA_EHATA_RELIABILITY)];
            case 'matlab'
                [pathLossMaps{idxCellAntenna}, ...
                    pathLossMapXLabels{idxCellAntenna}, ...
                    pathLossMapYLabels{idxCellAntenna}] ...
                    = genMedianBasicPropLossMaps( ...
                    CARRIER_FREQUENCY_IN_MHZ, ...
                    curBaseAntXY, curBaseAntHeightInM, ...
                    areaOfInterestXs, areaOfInterestYs, curRxAntH, ...
                    getEleFromXYFct, getLiDarZFromXYFct, REGION, ...
                    TERRAIN_RES_IN_M, LIBRARY_TO_USE);
                curFigFileName = [ ...
                    'CellTowerPathLoss_RxHeight_', num2str(curRxAntH), ...
                    '_TxCell_', num2str(idxCellAntenna), ...
                    '_eHataLib_', LIBRARY_TO_USE, ...
                    '_Terrain_', REGION];
        end
        
        if FLAG_GEN_FIGS
            % Generate a figure on Google map to show the path loss map.
            [curXs, curYs] ...
                = meshgrid(pathLossMapXLabels{idxCellAntenna}, ...
                pathLossMapYLabels{idxCellAntenna});
            curBoolsToShow ...
                = pathLossMaps{idxCellAntenna}(:) ...
                >=EXPECTED_PL_RANGE_IN_DB(1) ...
                & pathLossMaps{idxCellAntenna}(:) ...
                <=EXPECTED_PL_RANGE_IN_DB(2);
            [curLatsToShow, curLonsToShow] ...
                = utm2deg(curXs(curBoolsToShow), curYs(curBoolsToShow), ...
                repmat(UTM_ZONE, sum(curBoolsToShow), 1));
            
            hCurPLMap = figure;
            plot3k([curLonsToShow curLatsToShow ...
                pathLossMaps{idxCellAntenna}(curBoolsToShow)], ...
                'Labels', {'', ...
                'Longitude (degrees)', 'Latitude (degrees)', ...
                '', 'Path Loss (dB)'}, ...
                'ColorRange', EXPECTED_PL_RANGE_IN_DB);
            grid on; view(2); axis tight;
            plotGoogleMapAfterPlot3k(gcf, 'satellite');
            
            pathToSaveFig = fullfile(pathToSaveResults, ...
                [curFigFileName, '.png']);
            saveas(hCurPLMap, pathToSaveFig);
            
            hCurTransPLMap = figure;
            curAlpha = 0.7;
            curDotSize = 3;
            scatter(curLonsToShow, curLatsToShow, ...
                curDotSize.*ones(length(curLatsToShow),1), ...
                pathLossMaps{idxCellAntenna}(curBoolsToShow), ...
                'Marker', '.', ...
                'MarkerEdgeAlpha', curAlpha);
            curCB = colorbar; curCB.Label.String = 'Path Loss (dB)';
            xlabel('Longitude (degrees)');
            ylabel('Latitude (degrees)')
            grid on; view(2); axis tight;
            plot_google_map('MapType', 'satellite');
            
            pathToSaveFig = fullfile(pathToSaveResults, ...
                [curFigFileName, '_Transparent.png']);
            saveas(hCurTransPLMap, pathToSaveFig);
        end
    end
    
    % Save the results.
    towerPathLossMapsEHata{idxH} = pathLossMaps;
    towerPathLossMapsEHataXLabels{idxH} = pathLossMapXLabels;
    towerPathLossMapsEHataYLabels{idxH} = pathLossMapYLabels;
end

disp('    Done!')

%% Combine Path Loss Maps

disp(' ')
disp('        Combinning path loss maps ...')

for idxH = 1:numOfHs
    curRxAntH = RX_ANT_HEIGHTS_IN_M_FOR_EHATA(idxH);
    
    [coverageMapsEHata{idxH}, coverageMapsEHataXLabels{idxH}, ...
        coverageMapsEHataYLabels{idxH}] = combineTowerPathLossMaps( ...
        towerPathLossMapsEHata{idxH}, ...
        towerPathLossMapsEHataXLabels{idxH}, ...
        towerPathLossMapsEHataYLabels{idxH}, ...
        EXPECTED_PL_RANGE_IN_DB);
    
    if FLAG_GEN_FIGS
        % Generate a figure on Google map to show the path loss map.
        [curXs, curYs] ...
            = meshgrid(coverageMapsEHataXLabels{idxH}, ...
            coverageMapsEHataYLabels{idxH});
        [curLats, curLons] ...
            = utm2deg(curXs(:), curYs(:), ...
            repmat(UTM_ZONE, length(curXs(:)), 1));
        
        hCurPLMap = figure;
        plot3k([curLons curLats ...
            coverageMapsEHata{idxH}(:)], ...
            'Labels', {'', ...
            'Longitude (degrees)', 'Latitude (degrees)', ...
            '', 'Path Loss (dB)'});
        grid on; view(2); plotGoogleMapAfterPlot3k(gcf, 'satellite');
        
        pathToSaveFig = fullfile(pathToSaveResults, [ ...
            'Coverage_RxHeight_', num2str(curRxAntH), ...
            '_eHataLib_', LIBRARY_TO_USE, '.png']);
        saveas(hCurPLMap, pathToSaveFig);
    end
end

disp('    Done!')

%% Clear Things Up if Necessary
% This should be put at the very end of a script using the cpp ehata lib.
if libisloaded('ehata')
    unloadlibrary('ehata');
end

% Print out a ruler to indicate everything is done.
disp(fileNameHintRuler);

% EOF