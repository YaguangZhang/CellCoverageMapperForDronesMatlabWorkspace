function [LIDAR_DATA_SET_TO_USE] ...
    = verifyLidarDataSetToUse(lidarDataToUse, ABS_PATH_TO_SHARED_FOLDER)
%VERIFYLIDARDATASETTOUSE Verify the label for the LiDAR data set to use.
%
% Inputs:
%   - lidarDataToUse
%     The user defined label for the LiDAR data set to use.
%   - ABS_PATH_TO_SHARED_FOLDER
%     The work folder for the simulation.
%
% Output:
%   - LIDAR_DATA_SET_TO_USE
%     The verified label for the LiDAR data set to use.
%
% Yaguang Zhang, Purdue, 09/25/2019

switch lower(lidarDataToUse)
    case ''
        dirToLidarFilesMatIn = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar', 'IN', 'metaInfo.mat');
        if exist(dirToLidarFilesMatIn, 'file')
            LIDAR_DATA_SET_TO_USE = 'IN';
        else
            dirToLidarFilesMatTippExtended ...
                = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
                'Lidar', 'Tipp_Extended', 'metaInfo.mat');
            if exist(dirToLidarFilesMatTippExtended, 'file')
                LIDAR_DATA_SET_TO_USE = 'Tipp_Extended';
            else
                error(['None of the LiDAR data sets ', ...
                    'have been processed before. ', ...
                    'Please explicitly specify LIDAR_DATA_SET_TO_USE!']);
            end
        end
    case 'tipp_extended'
        LIDAR_DATA_SET_TO_USE = 'Tipp_Extended';
    case 'in'
        LIDAR_DATA_SET_TO_USE = 'IN';
    otherwise
        error(['Unsupported LIDAR_DATA_SET_TO_USE ', ...
            LIDAR_DATA_SET_TO_USE, '!']);
end

end
% EOF