function [ hostname ] = parseCacheFilename(cacheFilename)
%PARSECACHEFILENAME Parse the name of a simulation cache file for more
%information.
%
% Example cache file name:
%   CovAnalysisCache_Task_Tipp_LidarSet_IN_DSM_2019_NumPix_100_TerProRes_1.5_LidProRes_1.5_HostName_ygzhang.mat
%
% Yaguang Zhang, Purdue, 01/24/2022

curLabel = 'HostName';
tokens = regexp(cacheFilename, [curLabel, '_([^_.]+)'], 'tokens');
assert(length(tokens)==1 ...
    && length(tokens{1})==1, ...
    'More than one match found!');
hostname = tokens{1}{1};

end
% EOF