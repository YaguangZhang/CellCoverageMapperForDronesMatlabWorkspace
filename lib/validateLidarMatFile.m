function isValidMat = validateLidarMatFile(absPathToLidarFile, simConfigs)
%VALIDATELIDARMATFILE Validate whether the input dir is pointing a valid
%result .mat file for the LiDAR file. If not, replace it with the correct
%one.
%
% Inputs:
%   - absPathToLidarFile
%     The absolute full path to the .mat file to validate.
%   - simConfigs
%     The configuration struct for the simulation. We need fields
%     deg2utm_speZone and utm2deg_speZone for converting between GPS (lat,
%     lon) and UTM (x, y).
%
% Output:
%   - isValidMat
%     A flag indicating whether the file is valid or not.
%
% Yaguang Zhang, Purdue, 10/29/2019

% Check the file.
isValidMat = exist(absPathToLidarFile, 'file');
if isValidMat
    try
        isValidMat = matcat(absPathToLidarFile);
    catch
        isValidMat = false;
    end
end

if ~exist('simConfigs', 'var')
    simConfigs = evalin('base', 'simConfigs');
end

if ~isValidMat
    % Locate the Matlab workspace.
    [absPathToLib, ~] = fileparts(mfilename('fullpath'));
    
    % Create the cache folder.
    absPathToCacheFolder = fullfile(absPathToLib, ...
        ['cache_', char(java.util.UUID.randomUUID.toString)]);
    
    % Load the function for preprocessing the LiDAR data.
    if ~exist(absPathToCacheFolder, 'dir')
        mkdir(absPathToCacheFolder);
    end
    
    % Create the LiDAR dataset file structure.
    [absPathToLidarFolder, lidarFilename, ~] ...
        = fileparts(absPathToLidarFile);
    [~, lidarFileParentFolderName, ~] ...
        = fileparts(absPathToLidarFolder);
    absPathToCacheFolderLidarData = fullfile(absPathToCacheFolder, ...
        lidarFileParentFolderName);
    if ~exist(absPathToCacheFolderLidarData, 'dir')
        mkdir(absPathToCacheFolderLidarData);
    end
    
    % Copy the raw .tif LiDAR data to the cache folder.
    
    copyfile(fullfile(absPathToLidarFolder, [lidarFilename, '.tif']), ...
        absPathToCacheFolderLidarData);
    copyfile(fullfile(absPathToLidarFolder, [lidarFilename, '.img']), ...
        absPathToCacheFolderLidarData);
    
    % Fix the LiDAR file by reprocessing the cooresponding LiDAR raw data.
    flagPreProcessLidarFctAvailable ...
        = which('preprocessIndianaLidarDataSet');
    if isempty(flagPreProcessLidarFctAvailable)
        addpath(fullfile(absPathToLib, 'lidar'));
    end
    preprocessIndianaLidarDataSet(absPathToCacheFolder, ...
        simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone);
    
    % Copy the resulting .mat file back to replace the corrupted one.
    copyfile(fullfile(absPathToCacheFolderLidarData, ...
        [lidarFilename, '.mat']), ...
        absPathToLidarFile);
    
    % Delete the cached folder.
    rmdir(absPathToCacheFolder, 's');
    
    % Remove the library if necessary.
    if isempty(flagPreProcessLidarFctAvailable)
        rmpath(fullfile(absPathToLib, 'lidar'));
    end
end
end
% EOF