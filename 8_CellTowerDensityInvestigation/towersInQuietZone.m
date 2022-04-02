% TOWERSINQUIETZONE Plot an overview map for the cell towers in NRQZ.
%
% The The National Radio Quiet Zone (NRQZ) is bounded by NAD-83 meridians
% of longitude at 78d 29m 59.0s W and 80d 29m 59.2s W and latitudes of 37d
% 30m 0.4s N and 39d 15m 0.4s N
%
% Yaguang Zhang, Purdue, 04/01/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); addpath('.');
cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Settings

% NRQZ Bound
nrqzLatLimDms=[37 30 0.4; 39 15 0.4];
nrqzLonLimDms=[-80 29 59.2; -78 29 59.0];
nrqzLatLim = dms2degrees(nrqzLatLimDms);
nrqzLonLim = dms2degrees(nrqzLonLimDms);
nrqzLonLatBound = [nrqzLonLim(1), nrqzLatLim(1);
    nrqzLonLim(2), nrqzLatLim(1);
    nrqzLonLim(2), nrqzLatLim(2);
    nrqzLonLim(1), nrqzLatLim(2);
    nrqzLonLim(1), nrqzLatLim(1)];

% Path to save the plots.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', '8_CellTowerDensityInvestigation', ...
    'TowersInNRQZ');
% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

% For plotting.
areaOfInterestColor = [0.9290 0.6940 0.1250];
lightBlue = [0.3010 0.7450 0.9330];
darkBlue = [0 0.4470 0.7410];
colorEffectiveTowers = 'b';
markerEffectiveTowers = '.';
markerSizeEffectiveTowers = 10;
colorIneffectiveTowers = 'r';
markerIneffectiveTowers = 'x';
lineWidthIneffectiveTowers = 1.5;

%% Fetch Cell Towers

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Fetching tower locations ...'])

% Default to the NTIA+HIFLD cell tower locations.
ABS_PATH_TO_CELL_ANTENNAS_CSV = fullfile( ...
    ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'NtiaLayoutPlusHifldCellTs', ...
    'NtiaLayoutMergedWithHifldCellTs_Threshold_1000m_LatLonH.csv');

% Note: we use "height" to indicate the vertical distance from the ground
% to the antenna; "elevation" to indicate the ground elevation; and
% "altitude" to indicate elevation+height.
cellAntsLatLonH = csvread( ...
    ABS_PATH_TO_CELL_ANTENNAS_CSV, 1, 1); %#ok<CSVRD>

%% Towers in NRQZ

% Check which towers are in/on the NRQZ boundary.
boolsTsInNrqz = InPolygon(cellAntsLatLonH(:,2), cellAntsLatLonH(:,1), ...
    nrqzLonLatBound(:,1), nrqzLonLatBound(:,2));
effeCellAntsLons = cellAntsLatLonH(boolsTsInNrqz, 2);
effeCellAntsLats = cellAntsLatLonH(boolsTsInNrqz, 1);

inEffeCellAntsLons = cellAntsLatLonH(~boolsTsInNrqz, 2);
inEffeCellAntsLats = cellAntsLatLonH(~boolsTsInNrqz, 1);

%% Map Generation

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Generating plots ...'])

hFigOverview = figure; hold on;
hAreaOfInterest = plot(polyshape(nrqzLonLatBound), ...
    'FaceColor', areaOfInterestColor);
plot_google_map('MapType', 'hybrid');
axis manual;

hIneffeCells = plot(inEffeCellAntsLons, ...
    inEffeCellAntsLats, markerEffectiveTowers, ...
    'MarkerSize', markerSizeEffectiveTowers, ...
    'Color', colorIneffectiveTowers);
hEffeCells = plot(effeCellAntsLons, effeCellAntsLats, ...
    markerEffectiveTowers, ...
    'MarkerSize', markerSizeEffectiveTowers, ...
    'Color', colorEffectiveTowers);

xlabel('Longtitude (Degree)'); xlabel('Latitude (Degree)');
saveas(hFigOverview, fullfile(pathToSaveResults, 'tsInNrqz.fig'));
saveas(hFigOverview, fullfile(pathToSaveResults, 'tsInNrqz.jpg'));

zoom(0.5);
saveas(hFigOverview, fullfile(pathToSaveResults, 'tsInNrqz_Zoom_0_5.jpg'));

zoom(0.5);
saveas(hFigOverview, fullfile(pathToSaveResults, 'tsInNrqz_Zoom_0_25.jpg'));

disp(' ')
disp(['    [', datestr(now, datetimeFormat), ...
    '] Done!'])

% EOF