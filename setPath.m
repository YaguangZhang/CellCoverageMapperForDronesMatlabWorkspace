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
%  - Local copy on the computer cluster at Purdue:
absHomePathLinuxCoverageOnFrankie = ['/home/coverage/nvme/', ...
    '/CellCoverageMapper'];

% The absolute path to Python 3. Please make sure it is correct for the
% machine to run this script.
%  - On (quite powerful) Windows Artsy:
absPythonPathWinArtsy ...
    = 'C:\Python39\python.exe';
%  - Local copy on the computer cluster at Purdue:
absPythonPathLinuxCoverage = '/usr/bin/python3.7';
%  - Local copy on the computer cluster at Purdue:
if verLessThan('matlab','9.12')
    absPythonPathLinuxCoverageOnFrankie = '/usr/bin/python3.7';
else
    absPythonPathLinuxCoverageOnFrankie = '/usr/bin/python3.8';
end

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
    case 'ygzhang'
        % The virtual machine coverage on Purdue GPU cluster Frankie.
        ABS_PATH_TO_SHARED_FOLDER = absHomePathLinuxCoverageOnFrankie;
        ABS_PATH_TO_PYTHON = absPythonPathLinuxCoverageOnFrankie;
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% Make sure only the right NTIA library is added to the path.
if ispc
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtiaLinux'));
elseif isunix
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtia'));
end

% We need Python for concurrent HTTP requests to get elevation data from
% USGS faster. Make sure Python and its lib folder is added to path.

% Make sure Python is available and it is of the correct version. Note:
%   pyenv('Version', ABS_PATH_TO_PYTHON, 'ExecutionMode', 'OutOfProcess');
% may help fix the "Segmentation Failure" error.
curPyEnv = pyenv;
if ~strcmp(curPyEnv.Executable, ABS_PATH_TO_PYTHON)
    curPyEnv = pyenv('Version', ABS_PATH_TO_PYTHON);
    curPyVersionStr = char(curPyEnv.Version);
    assert(strcmp(curPyVersionStr(1), '3'), ...
        ['Loaded Python is not version 3.x!', ...
        ' Please restart Matlab and try again!']);
end
% Make sure our Python module is available.
try
    py_addpath(fullfile(pwd, 'lib', 'python'));
catch err
    warning(['Error identifier: ', err.identifier]);
    warning(['Error message: ',err.message]);
    errorMsg = 'Unable to set Python path! ';
    if isunix
        errorMsg = [errorMsg, ...
            'Please make sure both python3.x and ', ...
            'python3.x-dev are installed!'];
    end
    error(errorMsg);
end

% EOF