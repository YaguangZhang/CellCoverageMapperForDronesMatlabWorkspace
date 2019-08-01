% GENCOVERAGEMAPS Generate cellular coverage maps.
%
% This is the main script listing all the processing steps. Please comment
% & uncomment commands as it is needed, depending on which results require
% updates.
%
% Yaguang Zhang, Purdue, 06/10/2019

clear; clc; close all;
cd(fileparts(mfilename('fullpath')));

%% 1_CoverageMapsForAcre: Generage Coverage Maps for ACRE
addpath(fullfile(pwd, '1_CoverageMapsForAcre'));
genCoverageMapForAcre;

%% 2_CoverageMapsForTipp: Generage Coverage Maps for Tippecanoe
addpath(fullfile(pwd, '2_CoverageMapsForTipp'));
genCoverageMapForTipp;

%% 3_BlockageMapsForTipp: Generage Blockage Maps for Tippecanoe
addpath(fullfile(pwd, '3_BlockageMapsForTipp'));
genBlockageMapForTipp;

%% 4_PerformanceComparison: Compare the Maps
addpath(fullfile(pwd, '4_PerformanceComparison'));
evaluatePerformance;

% EOF