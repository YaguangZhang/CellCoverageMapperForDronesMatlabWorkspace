% GENCOVERAGEMAPS Generate cellular coverage maps.
%
% This is the main script listing all the processing steps. Please comment
% & uncomment commands as it is needed, depending on which results require
% updates.
%
% Yaguang Zhang, Purdue, 06/10/2019

clear; clc; close all;

% Add libs to current path and set ABS_PATH_TO_SHARED_FOLDER according to
% the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(pwd)));
setPath;

%% 1_Calibration: Calibrate the Gnu Radio RX
addpath(fullfile(pwd, '1_CoverageMapsForAcre'));
genCoverageMapForAcre;

% EOF