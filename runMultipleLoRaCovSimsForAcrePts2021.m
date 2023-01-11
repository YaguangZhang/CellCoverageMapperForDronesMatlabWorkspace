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
%     Same as 'AcreLoRaWan2021', except that the data source is Purdue Ag
%     IT, which powers ThingsBoard and is in fact more complete (with less
%     missing points).
%   - 'AcreLoRaWan2021AgITCleaned'
%     Same as 'AcreLoRaWan2021AgIT', but (i) with some invalid records
%     removed (e.g., a GPS sensor on a vehicle for another project is
%     ignored), and (ii) we would also want to remove power lines (via a
%     hand-drawn polygon) near ACRE in the simulations.
SIM_GROUP_PRESET = 'AcreLoRaWan2021AgITCleaned';

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
% Carrier frequency.
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
    case lower('AcreLoRaWan2021AgITCleaned')
        ABS_PATH_TO_RX_LOCS = fullfile( ...
            ABS_PATH_TO_SHARED_FOLDER, ...
            'AcreRxInfo', '2021', 'oyster_data_ag_it_cleaned.csv');
        % Also include a polygon for ignoring LiDAR z values. For any LiDAR
        % z profile points in (inside or on the edge of) this polygon, the
        % LiDAR z values will be replaced with ground elevation (from the
        % cached digital terrain models).
        CustomSimMetaDefault.polyLatLonsToUseGroundEle = [ ...
            40.470880933942098 -87.013911394462895
            40.470786314240300 -87.013592862673406
            40.473407922992202 -87.011985035545493
            40.475997122720003 -87.010383275689705
            40.478733904435600 -87.008678372206901
            40.479709980873302 -87.008077712260999
            40.481768232660500 -87.006791450558694
            40.484620145404499 -87.005019807081894
            40.485639974755898 -87.004400945319404
            40.487578069277902 -87.003181423610997
            40.491057269974000 -87.001009340170199
            40.494407099131699 -86.999316571231802
            40.497812125967400 -86.997587398660201
            40.501590672962500 -86.995688342467602
            40.503567520255103 -86.994666007105096
            40.505959496781401 -86.993434350852397
            40.509689149175799 -86.991550462840195
            40.509734125165302 -86.991914499171102
            40.506669922752103 -86.993449519032794
            40.505273282733100 -86.994163940332200
            40.504143032525597 -86.994725163008894
            40.502754413340199 -86.995435033853994
            40.501478795965397 -86.996078164705196
            40.499829453636998 -86.996939717354905
            40.499181239467902 -86.997252181872199
            40.497594126378097 -86.998065196344399
            40.496316103955699 -86.998711360831706
            40.494472848237201 -86.999654821655795
            40.492666454645203 -87.000580080663397
            40.491046888012100 -87.001429498768701
            40.489496497009902 -87.002379026865000
            40.488322144270803 -87.003110133162807
            40.486434831004999 -87.004299318510306
            40.482916165183902 -87.006471401951003
            40.479724979573000 -87.008446299045900
            40.476021352727102 -87.010754896110797
            40.473847540143197 -87.012089695990596
            40.470880933942098 -87.013911394462895];
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
    [~, lidarMatFName, ~] ...
        = fileparts(lidarMatFileAbsDirs{idxMatF});
    lidarMatFileAbsDirs{idxMatF} = fullfile(dirToLidarFiles, ...
        'MatlabCache', [lidarMatFName, '.mat']);

    assert(exist(lidarMatFileAbsDirs{idxMatF}, 'file'));
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