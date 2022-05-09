%TESTPYTHONELEVATIONLIB Test functions unber lib which uses
%lib/python/ConcurrentWebreader.py for fetching ground elevation values
%from USGS.
%
% Yaguang Zhang, Purdue, 05/09/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, 'Tests');

diary(fullfile(pathToSaveResults, 'TestPythonElevationLib.log'));

%% Test Case: IN Boundary

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Testing queryElevationPointsFromUsgsInChunks with IN boundary ...'])

disp(['    [', datestr(now, datetimeFormat), ...
    '] Loading IN boundary ...'])
inBoundaryLatLons = loadInBoundary;

disp(['    [', datestr(now, datetimeFormat), ...
    '] Fetching 1 USGS record ...'])

tic;
eles1 = queryElevationPointsFromUsgsInChunks( ...
    inBoundaryLatLons(1,1), inBoundaryLatLons(1,2));
toc;

disp(['    [', datestr(now, datetimeFormat), ...
    '] Fetching 50 USGS record ...'])

tic;
eles2 = queryElevationPointsFromUsgsInChunks( ...
    inBoundaryLatLons(2:51,1), inBoundaryLatLons(2:51,2));
toc;

disp(['    [', datestr(now, datetimeFormat), ...
    '] Fetching 1000 USGS record ...'])

tic;
eles3 = queryElevationPointsFromUsgsInChunks( ...
    inBoundaryLatLons(1000:1999,1), inBoundaryLatLons(1000:1999,2));
toc;

disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

diary off;

% EOF