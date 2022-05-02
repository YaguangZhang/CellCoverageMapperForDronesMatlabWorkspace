% TESTCOMPBLKPLANDDIST A snippet for optimizing the speed of computing
% blockage path loss and distance.
%
% One could pause in the function computeBlockageAndCoveragePLs during a
% cellular coverage simulation to prepare the inputs.
%
% Yaguang Zhang, Purdue, 05/02/2022

profile on;

for idxTest = 1:1000
    approximateBlockageAndCoveragePLs( ...
        lidarProfile, terrainProfile, ...
        startXYH, endXYH, ...
        simConfigs);
end

for idxTest = 1:1000
    computeBlockageAndCoveragePLs( ...
        lidarProfile, terrainProfile, ...
        startXYH, endXYH, ...
        simConfigs);
end

profile viewer;

% EOF