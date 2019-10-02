function [ utmXyBoundary ] ...
    = generateUtmXyBoundaryOfLidarDataset(dirToLidarFiles, simConfigs)
%GENERATEUTMXYBOUNDARYOFLIDARDATASET Generate the UTM (x, y) boundary
%polygon matrix for the input LiDAR dataset.
%
% Inputs:
%   - dirToLidarFiles
%     The absolute path to the LiDAR data set.
%   - simConfigs
%     A struct with the simulation configuration information. We need
%     fields:
%       - deg2utm_speZone, utm2deg_speZone
%        Functions to convert GPS degrees (lat, lon) from/to UTM (x, y). We
%        will polulate them later.
%
% Output:
%   - utmXyBoundary
%     The output polygon matrix, with each row being one vertex in the form
%     of UTM (x, y).
%
% Yaguang Zhang, Purdue, 10/02/2019

dirToLib = fileparts(mfilename('fullpath'));
dirToLidarPreprocessorFolder = fullfile(dirToLib, 'lidar');

% We need a function that may not be on the Matlab path. After using it, we
% will remove it from the Matlab path if necessary to leave the Matlab path
% unchanged.
flagLidarPreprocessorAvailable ...
    = exist('preprocessIndianaLidarDataSet', 'file')==2;
if ~flagLidarPreprocessorAvailable
    addpath(dirToLidarPreprocessorFolder);
end

% Preprocess the input LiDAR dataset.
[~, lidarFileXYCoveragePolyshapes, ~] ...
    = preprocessIndianaLidarDataSet(dirToLidarFiles, ...
    simConfigs.deg2utm_speZone, simConfigs.utm2deg_speZone);

% Remove the path to the LiDAR data preprocessor if neccesary.
if ~flagLidarPreprocessorAvailable
    rmpath(dirToLidarPreprocessorFolder);
end

% Overall boundries for the area covered by the LiDAR data set in UTM.
lidarFilesXYCoveragePolyshape ...
    = mergePolygonsForAreaOfInterest(lidarFileXYCoveragePolyshapes, 1);

% The area of interest is default to the area covered by the input LiDAR
% data.

[lidarFilesXYCoveragePolyshapeXs, lidarFilesXYCoveragePolyshapeYs] ...
    = boundary(lidarFilesXYCoveragePolyshape);
utmXyBoundary = [lidarFilesXYCoveragePolyshapeXs, ...
    lidarFilesXYCoveragePolyshapeYs];

assert(length(regions( ...
    polyshape(utmXyBoundary)))==1, ...
    'The area of interest should have only one region!');

end
% EOF