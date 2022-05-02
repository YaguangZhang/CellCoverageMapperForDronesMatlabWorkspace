% TESTCOMPBLKPLANDDIST A snippet for optimizing the speed of computing
% blockage path loss and distance.
%
% One could pause in the function computeBlockageAndCoveragePLs during a
% cellular coverage simulation to prepare the inputs.
%
% Yaguang Zhang, Purdue, 05/02/2022

cd(fileparts(mfilename('fullpath')));

addpath('.'); addpath('FctsCompBlkPLAndDist');
load('testCompBlkPLAndDist.mat');

cd('..'); addpath('lib');

profile on;
for idxTest = 1:1000
    [blockagePL, blockageDistInM] ...
        = computeBlockagePLAndDist_Pt2Line(txXYAlt, rxXYAlt, ...
        lidarProfile, simConfigs);
end
profsave;

profile on;
for idxTest = 1:1000
    [blockagePL1, blockageDistInM1] ...
        = approximateBlockagePLAndDist(txXYAlt, rxXYAlt, ...
        lidarProfile, simConfigs);
end
profsave;

profile on;
for idxTest = 1:1000
    [blockagePL2, blockageDistInM2] ...
        = computeBlockagePLAndDist(txXYAlt, rxXYAlt, ...
        lidarProfile, simConfigs);
end
profsave;

disp(['Error of blockageDistInM1 (approximateBlockagePLAndDist): ', ...
    num2str(blockageDistInM1-blockageDistInM)]);
disp(['Error of blockageDistInM2 (computeBlockagePLAndDist): ', ...
    num2str(blockageDistInM2-blockageDistInM)]);

% EOF