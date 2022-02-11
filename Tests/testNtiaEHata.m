%TESTNTIAEHATA Test the eHata library implemented by NTIA.
%
% Ref: https://github.com/NTIA/ehata
%
% Note: On Windoes machines, files under lib\ext\eHataNtia will be used; on
% Linux machines, files under lib\ext\eHataNtiaLinux will be used.
%
% Yaguang Zhang, Purdue, 02/11/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Constants for the tests.
simConfigs.CARRIER_FREQUENCY_IN_MHZ = 1900;
simConfigs.LOS_FIRST_FRES_CLEAR_RATIO = 0.6;
simConfigs.NTIA_EHATA_ENVIRO_CODE = 82; % Cultivated Crops.
simConfigs.NTIA_EHATA_RELIABILITY = 0.95;

simConfigs.CARRIER_WAVELENGTH_IN_M...
    = physconst('LightSpeed')/simConfigs.CARRIER_FREQUENCY_IN_MHZ/1e6;

% For reproducibility.
rng(1);

% Expected results.
refs = [168.243071206611, 177.628295224296, 161.393628486599];
maxAllowedError = 10^(-10);

%% Tests Partially Based on ShrinkedIN Simulations
% Three different tower&use loc pairs + random terrain profiles.

% Test 1
mu = 1000;
sigma = 150;

txXYH = [533512.500684368, 4211061.64612324, 105.2];
rxXYH = [483574.728695754, 4253672.80960810, 1.5];
numOfPtsInElePro = 100;

elePro = generateProfileForEHata(txXYH, rxXYH, ...
    normrnd(mu, sigma, [numOfPtsInElePro, 1]));
txEle = elePro(end);
elePro(3) = elePro(end);
elePro(end) = txEle;
pl1 = computeCoveragePL(txXYH, rxXYH, ...
    elePro, simConfigs);

% Test 2
mu = 2500;
sigma = 400;

txXYH = [472110.424576777, 4353972.06066561, 87];
rxXYH = [471237.117000556, 4391091.77793361, 15];
numOfPtsInElePro = 200;

elePro = generateProfileForEHata(txXYH, rxXYH, ...
    normrnd(mu, sigma, [numOfPtsInElePro, 1]));
pl2 = computeCoveragePL(txXYH, rxXYH, ...
    elePro, simConfigs); 

% Test 3
mu = 1000;
sigma = 150;

txXYH = [518720.910665658, 4377861.19303991, 50];
rxXYH = [518603.718510522, 4374894.19685083, 50];
numOfPtsInElePro = 50;

elePro = generateProfileForEHata(txXYH, rxXYH, ...
    normrnd(mu, sigma, [numOfPtsInElePro, 1]));
pl3 = computeCoveragePL(txXYH, rxXYH, ...
    elePro, simConfigs);

results = [pl1, pl2, pl3];
assert(all(refs-results<=maxAllowedError), 'Tests failed!')
disp('    Tests passed!')

% EOF