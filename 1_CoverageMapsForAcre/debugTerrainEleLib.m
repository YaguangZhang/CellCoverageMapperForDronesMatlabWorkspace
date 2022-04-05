% DEBUGTERRAINELELIB For fixing the bug (shifted tiles) in Terrain
% Elevation Library under Linux.
%
% Yaguang Zhang, Purdue, 03/25/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', 'DebuggingTerrainEleLib');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

% Save the machine type and Matlab version in output figure file names.
curVer = strrep(version, '.', '-');
curEnv = computer;

% Plot only the part of the data that will be visible.
axisToSet = [-87.0145, -87.0084, 40.4971, 40.5007];

%% Test Case for fetchregion
% Note that on Frankie (Linux) the center field part does not match with
% the satellite image. However, they match perfectly on Artsy (Windows).

lowResGridLatRange = [40.467071825508839; 40.501647105371035];
lowResGridLonRange = [-87.015786097515942; -86.976458345550711];
regionRef = fetchregion(lowResGridLatRange, lowResGridLonRange);
rawRefElevData = regionRef.readelevation( ...
    lowResGridLatRange, lowResGridLonRange, 'sampleFactor', 1);

dispelev(rawRefElevData, 'mode', 'latlong'); plot_google_map;
saveas(gcf, fullfile(pathToSaveResults, ...
    [curEnv, '_', curVer, '_FetchRegion_Overview.jpg']));

[lons, lats] = meshgrid(rawRefElevData.longs, rawRefElevData.lats);
eles = rawRefElevData.elev;
indicesToShowLon = find(lons(1,:)>=axisToSet(1), 1): ...
    find(lons(1,:)<=axisToSet(2), 1, 'last');
indicesToShowLat = find(lats(:,1)<=40.5007, 1): ...
    find(lats(:,1)>=40.4971, 1, 'last');

figure;
surf(lons(indicesToShowLat,indicesToShowLon), ...
    lats(indicesToShowLat,indicesToShowLon), ...
    eles(indicesToShowLat,indicesToShowLon), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot_google_map('MapType', 'hybrid');
view(2); zlim([0, max(eles(:))]);
axis(axisToSet);
saveas(gcf, fullfile(pathToSaveResults, ...
    [curEnv, '_', curVer, '_FetchRegion_Surf.jpg']));

%% Test Case for fetchAnomalyRegion
% The mismatch is fixed for Frankie (Linux).

lowResGridLatRange = [40.467071825508839; 40.501647105371035];
lowResGridLonRange = [-87.015786097515942; -86.976458345550711];
regionRef = fetchAnomalyRegion(lowResGridLatRange, lowResGridLonRange);
rawRefElevData = regionRef.readelevation( ...
    lowResGridLatRange, lowResGridLonRange, 'sampleFactor', 1);

dispelev(rawRefElevData, 'mode', 'latlong'); plot_google_map;
saveas(gcf, fullfile(pathToSaveResults, ...
    [curEnv, '_', curVer, '_FetchAnomalyRegion_Overview.jpg']));

[lons, lats] = meshgrid(rawRefElevData.longs, rawRefElevData.lats);
eles = rawRefElevData.elev;
indicesToShowLon = find(lons(1,:)>=axisToSet(1), 1): ...
    find(lons(1,:)<=axisToSet(2), 1, 'last');
indicesToShowLat = find(lats(:,1)<=40.5007, 1): ...
    find(lats(:,1)>=40.4971, 1, 'last');

figure;
surf(lons(indicesToShowLat,indicesToShowLon), ...
    lats(indicesToShowLat,indicesToShowLon), ...
    eles(indicesToShowLat,indicesToShowLon), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot_google_map('MapType', 'hybrid');
view(2); zlim([0, max(eles(:))]);
axis(axisToSet);
saveas(gcf, fullfile(pathToSaveResults, ...
    [curEnv, '_', curVer, '_FetchAnomalyRegion_Surf.jpg']));

%% Generate and Export a Reference Figure

[eles, curR] = readgeoraster( ...
   fullfile(pwd, 'usgsdata', 'n41w088', ...
    'usgs_ned_13_n41w088_gridfloat.flt'), ...
    'CoordinateSystemType', 'geographic');
[lats, lons] = geographicGrid(curR);

indicesToShowLon = find(lons(1,:)>=axisToSet(1), 1): ...
    find(lons(1,:)<=axisToSet(2), 1, 'last');
indicesToShowLat = find(lats(:,1)<=40.5007, 1): ...
    find(lats(:,1)>=40.4971, 1, 'last');

figure;
surf(lons(indicesToShowLat,indicesToShowLon), ...
    lats(indicesToShowLat,indicesToShowLon), ...
    eles(indicesToShowLat,indicesToShowLon), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot_google_map('MapType', 'hybrid');
view(2); zlim([0, max(eles(:))]);
axis(axisToSet);

% Save the machine type and Matlab version in the file name.
curVer = strrep(version, '.', '-');
curEnv = computer;
saveas(gcf, fullfile(pathToSaveResults, ...
    [curEnv, '_', curVer, '_Ref.jpg']));

% EOF