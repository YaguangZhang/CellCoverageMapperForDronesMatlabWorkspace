% SETPATH Add lib folders into Matlab path.
%
% Yaguang Zhang, Purdue, 06/10/2019

cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(genpath(fullfile(pwd, 'lib')));

% A workaround for making Matlab R2019b work. If
% preprocessIndianaLidarDataSet.m is in the path, it has trouble using
% loaded functions as they are, even though R2019a works anyway.
rmpath(fullfile(pwd, 'lib', 'lidar'));

% The absolute path to the shared folder holding the data and code. Please
% make sure it is correct for the machine to run this script.
%  - On (quite powerful) Windows Artsy:
absHomePathWinArtsy = ['D:\One Drive - Purdue\OneDrive - purdue.edu', ...
    '\OATS\CellCoverageMapper'];
%  - Local copy on the computer cluster at Purdue:
absHomePathLinuxCoverage = ['/home/coverage', ...
    '/CellCoverageMapper'];

% The absolute path to Python 3. Please make sure it is correct for the
% machine to run this script.
%  - On (quite powerful) Windows Artsy:
absPythonPathWinArtsy ...
    = ['C:\Users\Yaguang Zhang\AppData\Local\Programs', ...
    '\Python\Python37\python.exe'];
%  - Local copy on the computer cluster at Purdue:
absPythonPathLinuxCoverage = ['/usr/bin/python3.7'];

unknownComputerErrorMsg = ...
    ['Compute not recognized... \n', ...
    '    Please update setPath.m for your machine. '];
unknownComputerErrorId = 'setPath:computerNotKnown';

[~, curHostname] = system('hostname');
switch strtrim(curHostname)
    case 'Artsy'
        % ZYG's lab desktop.
        ABS_PATH_TO_SHARED_FOLDER = absHomePathWinArtsy;
        ABS_PATH_TO_PYTHON = absPythonPathWinArtsy;
    case 'coverage-compute-big'
        % The computer cluster at Purdue.
        ABS_PATH_TO_SHARED_FOLDER = absHomePathLinuxCoverage;
        ABS_PATH_TO_PYTHON = absPythonPathLinuxCoverage;
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% Make sure only the right NTIA library is added to the path.
if ispc
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtiaLinux'));
elseif isunix
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtia'));
end

% Make sure Python and its lib folder is added to path. Make sure Python is
% available. Make sure Python is available.
curPythonVersion = pyversion;
if isempty(curPythonVersion) || (~strcmp(curPythonVersion(1:3), '3.7'))
    pyversion(ABS_PATH_TO_PYTHON);
end
% Check the version again.
curPythonVersion = pyversion;
if ~strcmp(curPythonVersion(1:3), '3.7')
    error(['Loaded Python is not version 3.7.', ...
        ' Please restart Matlab and try again!']);
end
try
    py_addpath(fullfile(pwd, 'lib', 'python'));
catch err
    warning(['Error identifier: ', err.identifier]);
    warning(['Error message: ',err.message]);
    errorMsg = 'Unable to set Python path! ';
    if isunix
        errorMsg = [errorMsg, ...
            'Please make sure both python3.7 and ', ...
            'python3.7-dev are installed!']; 
    end
    error(errorMsg);
end
% EOF