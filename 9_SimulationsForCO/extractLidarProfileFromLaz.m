%EXTRACTLIDARPROFILEFROMLAZ Snippet to extract a LiDAR profile for a long
%link (~14.5 km) in Colorado.
%
% Developed with Matlab R2022a.
%
% Yaguang Zhang, Purdue, 10/19/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

%% Configurations

% Path to the LiDAR laz file.
pathToLaz = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', 'CO', 'YZ_CA_CO_LiDAR', 'OneTileExample', 'co_laz.laz');

% Link of interest:
latLonStart = [39.991794, -105.274711];
latLonEnd = [40.121142, -105.250972];

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', 'CO_LongLink_LidarProfileForProfA');

if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

%% Load LiDAR Data

lasReader = lasFileReader(pathToLaz);
ptCloud = readPointCloud(lasReader);

crs = readCRS(lasReader);

% Convert UTM (x, y) to (lat, lon).
[lidarLats, lidarLons] = projinv(crs, ...
    ptCloud.Location(:, 1), ptCloud.Location(:, 2));

% Extract LiDAR z, too, for convenience.
lidarZsInM = ptCloud.Location(:, 3);

%% Construct Link of Interest in UTM

[utmXStart, utmYStart] = projfwd(crs, latLonStart(1), latLonStart(2));
[utmXEnd, utmYEnd] = projfwd(crs, latLonEnd(1), latLonEnd(2));

%% Overview Figure of the LiDAR Data

hOverviewMap = figure;
maxLidarZInM = max(lidarZsInM);
hLink = plot3([latLonStart(2), latLonEnd(2)], ...
    [latLonStart(1), latLonEnd(1)], ...
    [maxLidarZInM, maxLidarZInM], '-r');
view(2);
plot_google_map('MapType', 'hybrid');
plot3k([lidarLons, lidarLats, lidarZsInM], ...
    'Labels', {'Anemo LiDAR Data Overview', ...
    'Longitude (degree)', 'Latitude (degree)', 'LiDAR z (m)', ...
    'LiDAR z (m)'});
view(2);
legend(hLink, 'Link of Interest');

saveas(hOverviewMap, fullfile(pathToSaveResults, 'LidarOverview_Anemo'));

%% Extract Profile

%% Overview Figure of the Profile

% EOF