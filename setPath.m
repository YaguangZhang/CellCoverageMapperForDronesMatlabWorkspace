% SETPATH Add lib folders into Matlab path.
%
% Yaguang Zhang, Purdue, 06/10/2019

cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(genpath(fullfile(pwd, 'lib')));

% The absolute path to the shared folder holding the data and code. Please
% make sure it is correct for the machine to run this script.
%  - On (quite powerful) Windows Artsy:
absPathWinArtsy = ['D:\One Drive - Purdue\OneDrive - purdue.edu', ...
    '\OATS\CellCoverageMapper'];
%  - Local copy on Windows Dell:
absPathWinDell = ['C:\Users\Zyglabs\OneDrive - purdue.edu', ...
    '\OATS\CellCoverageMapper'];
%  - Local copy on the computer cluster at Purdue:
absPathWinDell = ['/home/coverage', ...
    '/CellCoverageMapper'];

unknownComputerErrorMsg = ...
    ['Compute not recognized... \n', ...
    '    Please update setPath.m for your machine. '];
unknownComputerErrorId = 'setPath:computerNotKnown';

[~, curHostname] = system('hostname');
switch strtrim(curHostname)
    case 'Artsy'
        % ZYG's lab desktop.
        ABS_PATH_TO_SHARED_FOLDER = absPathWinArtsy;
    case 'ZygLabs-Dell'
        % ZYG's Dell laptop.
        ABS_PATH_TO_SHARED_FOLDER = absPathWinDell;
    case 'coverage-compute'
        % The computer cluster at Purdue.
        ABS_PATH_TO_SHARED_FOLDER = absPathWinDell;
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% Make sure only the right NTIA library is added to the path.
if ispc
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtiaLinux'));
elseif isunix
    rmpath(fullfile(pwd, 'lib', 'ext', 'eHataNtia'));
end

% EOF