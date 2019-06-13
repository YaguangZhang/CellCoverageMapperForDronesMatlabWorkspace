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
unknownComputerErrorMsg = ...
    ['Compute not recognized... \n', ...
    '    Please update setPath.m for your machine. '];
unknownComputerErrorId = 'setPath:computerNotKnown';
switch getenv('computername')
    case 'ARTSY'
        % ZYG's lab desktop.
        ABS_PATH_TO_SHARED_FOLDER = absPathWinArtsy;
    case 'ZYGLABS-DELL'
        % ZYG's Dell laptop.
        ABS_PATH_TO_SHARED_FOLDER = absPathWinDell;    
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end
% EOF