% TESTFRESNELZON Visualize Fresnel Zone radius for different carrier
% frequencies and TX-to-RX distances.
%
% Yaguang Zhang, Purdue, 04/27/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, '..', ...
    'PostProcessingResults', 'Test_FresnelZone');
if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

%% Evaluate the Max Fresnel Zone Radii

carrierFreqsInMHz = 900:5:28000;
tx2RxDistsInM = 10:10:30000;

% We will store the radii values into a matrix, with rows for different
% freuqencies and columns for different distances.
numOfFreqs = length(carrierFreqsInMHz);
numOfDists = length(tx2RxDistsInM);
max1stFresnelZoneRadiiMatrix = nan(numOfFreqs, numOfDists);

halfTx2RxDistsInM = tx2RxDistsInM./2;
for idxFreq = 1:numOfFreqs
    curFreqInMHz = carrierFreqsInMHz(idxFreq);
    curWaveLengthInM = physconst('LightSpeed')/curFreqInMHz/1e6;
    curFctCalcFirstFresRadii = @(d1s, d2s) sqrt( ...
        (curWaveLengthInM .* d1s .* d2s)./(d1s + d2s));

    max1stFresnelZoneRadiiMatrix(idxFreq,:) = curFctCalcFirstFresRadii( ...
        halfTx2RxDistsInM, halfTx2RxDistsInM);
end

%% Plot

[XInKm, YInGHz] = meshgrid(tx2RxDistsInM./1000, carrierFreqsInMHz./1000);
hRadiiFig = figure('Position', [0, 0, 750, 750]);
surf(XInKm, YInGHz, max1stFresnelZoneRadiiMatrix, 'EdgeColor', 'interp');
view(2); colormap hot; colorbar;
xlabel('TX-to-RX Distance (km)'); ylabel('Carrier Frequency (GHz)');
title('Max 1st Fresnel Zone Radius (m)');
axis tight; axis equal;
saveas(hRadiiFig, fullfile(pathToSaveResults, ...
    'Max1stFresnelZoneRadii.jpg'));

hRadiiFig = figure('Position', [0, 0, 750, 750]);
surf(XInKm, YInGHz, 0.6.*max1stFresnelZoneRadiiMatrix, 'EdgeColor', 'interp');
view(2); colormap hot; colorbar;
xlabel('TX-to-RX Distance (km)'); ylabel('Carrier Frequency (GHz)');
title('60% Max 1st Fresnel Zone Radius (m)')
axis tight; axis equal;
saveas(hRadiiFig, fullfile(pathToSaveResults, ...
    'Max1stFresnelZoneRadii_60Percent.jpg'));

% EOF