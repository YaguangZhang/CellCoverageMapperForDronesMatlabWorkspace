%RUNMULTIPLELORACOVSIMSFORACREPTS2021 Run multiple analyzeCellularCoverage
%for point-by-point simulation results for the 2021 ACRE LoRaWAN dataset.
%
% This script can also be used to update the figures for completed
% simulations.
%
% Yaguang Zhang, Purdue, 06/20/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Initialization

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Initializing ...'])

% Predefined simulation group label.
%   - 'AcreLoRaWan2021'
%     The ACRE LoRaWAN results collected in 2021: one tower gateway +
%     multiple vechiles. Based on ThingsBoard records.
%   - 'AcreLoRaWan2021AgIT'
%     Same as 'acreLoRaWan2021', except that the data source is Purdue Ag
%     IT, which powers ThingsBoard and is in fact more complete (with less
%     missing points).
SIM_GROUP_PRESET = 'AcreLoRaWan2021AgIT';

pathToPostProcessingResultsFolder ...
    = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    ['PostProcessingResults_', SIM_GROUP_PRESET]);
if ~exist(pathToPostProcessingResultsFolder, 'dir')
    mkdir(pathToPostProcessingResultsFolder);
end
pathToSaveSimManDiary = fullfile( ...
    pathToPostProcessingResultsFolder, 'simManDiary.txt');
diary(pathToSaveSimManDiary);

pathToSaveCustomSimMetas = fullfile( ...
    pathToPostProcessingResultsFolder, 'CustomSimMetas.mat');

% Preset of interest.
PRESET = 'CustomSim';

%------------------------
% Required for CustomSim
%------------------------
% Carrier frequencies.
%	- 1900 MHz
%     For cellular 4G LTE
%   - C-Band: 3700 MHz (band n77) and 4700 MHz (band n79)
%     For cellular 5G sub 6G
%   - 7000 MHz (7 GHz) and 13000 MHz (13GHz)
%     For broadband wireless backhaul.
%   - mmWave 28000 MHz (28 GHz)
%     For cellular 5G millimeter wave
CustomSimMetaDefault.CARRIER_FREQUENCY_IN_MHZ = 915;
% Tower location.
CustomSimMetaDefault.ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
    ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'PurdueAcreLoraWanTowers.csv');
% Area of interest is to be extended Tipp.
pathToStateInfoFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries');
extTippBoundary = load(fullfile(pathToStateInfoFolder, ...
    'Tipp_Extended', 'boundary.mat'));
CustomSimMetaDefault.UTM_X_Y_BOUNDARY_OF_INTEREST ...
    = extTippBoundary.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;

% Rx location information TBD.
CustomSimMetaDefault.mapGridLatLonPts = nan;
CustomSimMetaDefault.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = nan;

%------------------------
% Optional for CustomSim
%------------------------
% Enable simulations for vegetation.
CustomSimMetaDefault.FLAG_EVAL_LOSS_THROUGH_VEG = true;
% Dir to save results.
CustomSimMetaDefault.pathToPostProcessingResultsFolder ...
    = pathToPostProcessingResultsFolder;

% For feedback.
CustomSimMetaDefault.SIM_GROUP_PRESET = SIM_GROUP_PRESET;

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Create CustomSimMetas

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Setting up simulations ...'])

switch lower(SIM_GROUP_PRESET)
    case lower('AcreLoRaWan2021')
        ABS_PATH_TO_RX_LOCS = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'AcreRxInfo', '2021', 'oyster_data_for_sim.csv');
    case lower('AcreLoRaWan2021AgIT')
        ABS_PATH_TO_RX_LOCS = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'AcreRxInfo', '2021', 'oyster_data_ag_it.csv');
    otherwise
        error(['Unsupported SIM_GROUP_PRESET: ', SIM_GROUP_PRESET, '!']);
end

rxGidLonLatHInMs = readmatrix(ABS_PATH_TO_RX_LOCS);
% Make sure the records loaded from the RX loc file are sorted by heights.
rxGidLonLatHInMs = sortrows(rxGidLonLatHInMs, 4);

rxHInMsUnique = unique(rxGidLonLatHInMs(:,4));
numOfSims = length(rxHInMsUnique);
CustomSimMetas = cell(numOfSims, 1);

for idxSim = 1:length(CustomSimMetas)
    curRxHInM = rxHInMsUnique(idxSim);

    CustomSimMetas{idxSim} = CustomSimMetaDefault;
    CustomSimMetas{idxSim}.mapGridLatLonPts = rxGidLonLatHInMs( ...
        rxGidLonLatHInMs(:,4)==curRxHInM, [3,2]);
    CustomSimMetas{idxSim}.RX_ANT_HEIGHTS_TO_INSPECT_IN_M = curRxHInM;
end

save(pathToSaveCustomSimMetas, ...
    'SIM_GROUP_PRESET', 'PRESET', 'ABS_PATH_TO_RX_LOCS', ...
    'CustomSimMetas', 'rxGidLonLatHInMs', 'rxHInMsUnique', 'numOfSims');

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Run Sims

for idxSim = 1:length(CustomSimMetas)
    CustomSimMeta = CustomSimMetas{idxSim};

    disp(' ')
    disp(['[', datestr(now, datetimeFormat), ...
        '] Running sim for ', CustomSimMeta.SIM_GROUP_PRESET, ...
        ' with PRESET ', PRESET, ' (', ...
        num2str(idxSim), ') ...'])

    try
        diary off;
        analyzeCellularCoverage;
        diary(pathToSaveSimManDiary);
    catch err
        diary(pathToSaveSimManDiary);
        disp(getReport(err))
        rethrow(err);
    end

    disp(['[', datestr(now, datetimeFormat), ...
        '] Done!'])
end

%% Aggregate Results

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Aggregating simulation results ...'])

% The simulation would clear some of the parameters we need.
load(pathToSaveCustomSimMetas);

% Retrieve simulation results.
[rxEHataPLsInDb, rxAccBlkDistsInM, ...
    rxAccBlkDistsInM_Terrain, rxAccBlkDistsInM_Veg] ...
    = deal(cell(numOfSims, 1));
for idxSim = 1:length(CustomSimMetas)
    CustomSimMeta = CustomSimMetas{idxSim};
    curRxHInM = rxHInMsUnique(idxSim);

    % Locate the simulation result folder based on the name pattern.
    curSimResultsFolderName = ['Simulation_', PRESET, ...
        '_Carrier_', num2str(CustomSimMeta.CARRIER_FREQUENCY_IN_MHZ), ...
        'MHz_LiDAR_*', '_idxSim_', num2str(idxSim)];

    if ~exist('pathToPostProcessingResultsFolder', 'var')
        pathToPostProcessingResultsFolder ...
            = CustomSimMeta.pathToPostProcessingResultsFolder;
    else
        assert(strcmp(pathToPostProcessingResultsFolder, ...
            CustomSimMeta.pathToPostProcessingResultsFolder), ...
            'Results are expected to be saved in the same folder!');
    end

    curSimResultFolder = rdir( ...
        fullfile(pathToPostProcessingResultsFolder, ...
        curSimResultsFolderName), 'isdir');
    assert(length(curSimResultFolder)==1, ...
        'Exactly one sim result folder is expected!');

    dirToSimState = fullfile( ...
        curSimResultFolder.name, 'simState.mat');
    curSimState = load(dirToSimState, 'simState');

    assert(length(curSimState.simState.coverageMaps) == 1, ...
        'Expecting one and only one map!');

    rxEHataPLsInDb{idxSim} = curSimState.simState.coverageMaps{1};
    rxAccBlkDistsInM{idxSim} = curSimState.simState.blockageDistMaps{1};
    rxAccBlkDistsInM_Terrain{idxSim} ...
        = curSimState.simState.blockageByTerrainDistMaps{1};
    rxAccBlkDistsInM_Veg{idxSim} ...
        = curSimState.simState.blockageByVegDistMaps{1};

    if length(rxEHataPLsInDb{idxSim}) ...
            ~= sum(rxGidLonLatHInMs(:,4)==curRxHInM)
        disp(['Num of pts in current sim results: ', ...
            num2str(length(rxEHataPLsInDb{idxSim}))]);
        disp(['Num of pts expected based on RX height:', ...
            num2str(sum(rxGidLonLatHInMs(:,4)==curRxHInM))])
        error('Sim results have unexpected number of points!');
    end
end

% Arrange results into a table for easy exportation.
gid = rxGidLonLatHInMs(:,1);
lon = rxGidLonLatHInMs(:,2);
lat = rxGidLonLatHInMs(:,3);
install_height_m = rxGidLonLatHInMs(:,4);
ehata_path_loss_db = vertcat(rxEHataPLsInDb{:});
accumulate_blockage_dist_m = vertcat(rxAccBlkDistsInM{:});
accumulate_blockage_dist_m_terrain = vertcat(rxAccBlkDistsInM_Terrain{:});
accumulate_blockage_dist_m_vegetation = vertcat(rxAccBlkDistsInM_Veg{:});

% LiDAR Data Info
dirToLidarFiles = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar_2019', 'IN', 'DSM');

% Preprocess .img/.tif LiDAR data. To make Matlab R2019b work, we need to
% remove preprocessIndianaLidarDataSet from path after things are done.
addpath(fullfile(pwd, 'lib', 'lidar'));
[lidarFileRelDirs, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessIndianaLidarDataSetDsm(dirToLidarFiles, ...
    deg2utm_speZone, utm2deg_speZone);
rmpath(fullfile(pwd, 'lib', 'lidar'));
lidarFileAbsDirs = cellfun(@(d) ...
    [dirToLidarFiles, strrep(d, '\', filesep)], ...
    lidarFileRelDirs, 'UniformOutput', false);

% Extra information on the LiDAR data set.
%   - Overall boundries for the area covered by the LiDAR data set in UTM.
% lidarFilesXYCoveragePolyshape ...
%     = mergePolygonsForAreaOfInterest(lidarFileXYCoveragePolyshapes, 1);
%   - Centroids for the LiDAR files in UTM.
lidarFileXYCentroids ...
    = extractCentroidsFrom2DPolyCell(lidarFileXYCoveragePolyshapes);
%   - The .mat copies for the LiDAR data. For the 2019 dataset, they are
%   stored in a cache folder.
lidarMatFileAbsDirs = lidarFileAbsDirs;
for idxMatF = 1:length(lidarMatFileAbsDirs)
    [lidarMatFPath, lidarMatFName, ~] ...
        = fileparts(lidarMatFileAbsDirs{idxMatF});
    lidarMatFileAbsDirs{idxMatF} = fullfile(lidarMatFPath, '..', ...
        'MatlabCache', [lidarMatFName, '.mat']);
end

% Extra information.
[utm_x, utm_y] = deg2utm_speZone(lat, lon);
utm_zone = repmat(simConfigs.UTM_ZONE, length(utm_x), 1);

dirToLoadRxLocEle = fullfile(pathToPostProcessingResultsFolder, ...
    'rxLocEleCache.mat');
if exist(dirToLoadRxLocEle, 'file')
    load(dirToLoadRxLocEle);
else
    ele_usgs_m = generateProfileSamps( ...
        [utm_x, utm_y], utm2deg_speZone, ...
        lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
        lidarMatFileAbsDirs, 'elevation');
    save(dirToLoadRxLocEle, 'ele_usgs_m');
end

% We have only one gateway in 2021.
towerUtmXYH = simState.CellAntsXyhEffective(1,:);

dist_2d_squared = (utm_x - towerUtmXYH(1)).^2 ...
    + (utm_y - towerUtmXYH(2)).^2;
dist_2d_m = sqrt(dist_2d_squared);
dist_3d_m = sqrt(dist_2d_squared + (install_height_m - towerUtmXYH(3)).^2);

simResultsTable = table(gid, lon, lat, install_height_m, ...
    ehata_path_loss_db, accumulate_blockage_dist_m, ...
    accumulate_blockage_dist_m_terrain, ...
    accumulate_blockage_dist_m_vegetation, ...
    utm_x, utm_y, utm_zone, dist_2d_m, ele_usgs_m, dist_3d_m);

writetable(simResultsTable, ...
    fullfile(pathToPostProcessingResultsFolder, ...
    ['simResults_', SIM_GROUP_PRESET, '.csv']));

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Debugging Plots

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Aggregating simulation results ...'])

% Path loss.
curData = ehata_path_loss_db;
colorRangeToSet = [60, 150];

curDataRange = [min(curData), max(curData)];
disp(['    Current data value range: ', ...
    num2str(curDataRange(1)), ', ', num2str(curDataRange(2))]);

assert(colorRangeToSet(1)<=curDataRange(1), ...
    colorRangeToSet(2)>=curDataRange(2), ...
    'colorRangeToSet does not cover curDataRange!');

figure;
plot3k([lon, lat, curData], 'ColorRange', colorRangeToSet, ...
    'Labels', {'eHata Path Loss Overview', ...
    'Longitude (degree)', 'Latitude (degree)', '', 'Path Loss (dB)'});
plot_google_map('MapType', 'hybrid');
zlim([0, colorRangeToSet(2)]); view(2);

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_eHata');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);

% LoS blockage.
curData = accumulate_blockage_dist_m;
colorRangeToSet = [0, 8000];

curDataRange = [min(curData), max(curData)];
disp(['    Current data value range: ', ...
    num2str(curDataRange(1)), ', ', num2str(curDataRange(2))]);

assert(colorRangeToSet(1)<=curDataRange(1), ...
    colorRangeToSet(2)>=curDataRange(2), ...
    'colorRangeToSet does not cover curDataRange!');

figure;
plot3k([lon, lat, curData], 'ColorRange', colorRangeToSet, ...
    'Labels', {'Accumulate LoS Blockage Distance Overview', ...
    'Longitude (degree)', 'Latitude (degree)', '', ...
    'Blockage Distance (m)'}); view(2);
plot_google_map('MapType', 'hybrid');
zlim([0, colorRangeToSet(2)]);

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_LoSBlockageDist');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);

curPointAlpha = 0.1;
figure;
hSca = scatter(dist_2d_m, accumulate_blockage_dist_m, ...
    'MarkerFaceColor','b','MarkerEdgeColor','b');
alpha(hSca, curPointAlpha);
set(gca, 'XScale', 'log', 'YScale', 'log');
title('Accumulate LoS Blockage Distance Over 2D Distance');
xlabel('2D distance (m)');
ylabel('Accumulate LoS blockage distance (m)');
grid on;

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_LoSBlockageDistOver2DDist');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);

figure;
hSca = scatter(dist_3d_m, accumulate_blockage_dist_m, ...
    'MarkerFaceColor','b','MarkerEdgeColor','b');
alpha(hSca, curPointAlpha);
set(gca, 'XScale', 'log', 'YScale', 'log');
title('Accumulate LoS Blockage Distance Over 3D Distance');
xlabel('3D distance (m)');
ylabel('Accumulate LoS blockage distance (m)');
grid on;

curPathToSaveFig = fullfile(pathToPostProcessingResultsFolder, ...
    'Overview_LoSBlockageDistOver3DDist');
saveas(gcf, [curPathToSaveFig, '.fig']);
saveas(gcf, [curPathToSaveFig, '.jpg']);
disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

diary off;

% EOF