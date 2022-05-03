%TESTINPOLY Test the inpolygon function and 3rd implementations for it.
%
% Note that we are using proprietory INDOT GPS data for the tests, so the
% data will not be pushed to GitHub. However, the resultant .log file will
% be pushed there for future reference.
%
% Three alternatives to the built-in inpolygon functions are considered:
%   (1) InPolygon MEX implemenation
%       https://www.mathworks.com/matlabcentral/fileexchange/20754-fast-inpolygon-detection-mex
%   (2) INPOLY: A fast points-in-polygon test (inpoly2)
%       https://www.mathworks.com/matlabcentral/fileexchange/10391-inpoly-a-fast-points-in-polygon-test
%   (3) A custom wrapper named "InPolygon" for inpoly2
%       ./InPolyLibs/InPolygonReplacer/InPolygon.m
%
% Yaguang Zhang, Purdue, 04/28/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(pwd, 'Tests');

% Remove project lib folder and load only the libraries to test.
rmpath(genpath(fullfile(pathToSaveResults, '..', 'lib')));
addpath(genpath(fullfile(pathToSaveResults, 'InPolyLibs', 'inpoly')));

diary(fullfile(pathToSaveResults, 'TestInPolyResults.log'));

%% Create Test Pts

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Loading test dataset ...'])

tic;
load(fullfile(pathToSaveResults, 'testInPoly.mat'));
toc;

% Attach the boundary vertices to the test dataset.
gpsLonLatCoors = [gpsLonLatCoors; inBoundaryLatLons(:,2:-1:1)];

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Speed Tests

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Speed tests ...'])

disp(['    [', datestr(now, datetimeFormat), ...
    '] inpolygon ...'])
tic;
boolsGpsLonLatCoorsOutOfIn = ~inpolygon(gpsLonLatCoors(:,1), ...
    gpsLonLatCoors(:,2), ...
    inBoundaryLatLons(:,2), inBoundaryLatLons(:,1));
toc;

disp(['    [', datestr(now, datetimeFormat), ...
    '] InPolygon ...'])
absPathToCurLib = fullfile(pathToSaveResults, ...
    'InPolyLibs', 'InPolygon-MEX');
addpath(genpath(absPathToCurLib));
tic;
boolsGpsLonLatCoorsOutOfIn1 = ~InPolygon(gpsLonLatCoors(:,1), ...
    gpsLonLatCoors(:,2), ...
    inBoundaryLatLons(:,2), inBoundaryLatLons(:,1));
toc;
rmpath(genpath(absPathToCurLib));

disp(['    [', datestr(now, datetimeFormat), ...
    '] inpoly2 x100 ...'])
tic;
for idx = 1:100
    boolsGpsLonLatCoorsOutOfIn2 = ~inpoly2(gpsLonLatCoors, ...
        inBoundaryLatLons(:,2:-1:1));
end
toc;

disp(['    [', datestr(now, datetimeFormat), ...
    '] InPolygon replacer (inpoly2) x100 ...'])
absPathToCurLib = fullfile(pathToSaveResults, ...
    'InPolyLibs', 'InPolygonReplacer');
addpath(genpath(absPathToCurLib));
tic;
for idx = 1:100
    boolsGpsLonLatCoorsOutOfIn3 = ~InPolygon(gpsLonLatCoors(:,1), ...
        gpsLonLatCoors(:,2), ...
        inBoundaryLatLons(:,2), inBoundaryLatLons(:,1));
end
toc;
rmpath(genpath(absPathToCurLib));

disp(' ')
disp(['[', datestr(now, datetimeFormat), ...
    '] Done!'])

%% Verify Anomaly Results with the Built-In inpolygon

% Note: we have visualized the inpoly2 output boolsGpsLonLatCoorsOutOfIn2
% and it is actually correct; here we just want to double check the output
% of inpolygon agrees with inpoly2 when there is a disagreement between
% InPolygon and inpoly2.

disp('=========')
disp(' Summary ')

whos gpsLonLatCoors inBoundaryLatLons;

disp(['Num of misclassified pts - InPolygon: ', num2str(sum( ...
    boolsGpsLonLatCoorsOutOfIn1~=boolsGpsLonLatCoorsOutOfIn))]);
disp(['Num of misclassified pts - inpoly2: ', num2str(sum( ...
    boolsGpsLonLatCoorsOutOfIn2~=boolsGpsLonLatCoorsOutOfIn))]);
disp(['Num of misclassified pts - InPolygon replacer: ', num2str(sum( ...
    boolsGpsLonLatCoorsOutOfIn3~=boolsGpsLonLatCoorsOutOfIn))]);

disp('=========')

% EOF
