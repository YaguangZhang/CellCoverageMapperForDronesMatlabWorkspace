function [terrainProfile, lidarProfile] ...
    = fetchTerrainAndLidarProfiles(absPathToCacheMatFile, ...
    startXY, endXY, simConfigs, ...
    lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
    lidarMatFileAbsDirs)
%FETCHTERRAINANDLIDARPROFILES Generate the terrain and LiDAR profiles
%between the start and end locations.
%
% We will use the local LiDAR/elevation data set (lidarMatFileAbsDirs)
% first, then try filling locations where the local dataset does not cover
% with online elevation data from USGS via the National Map - Elevation
% Point Query Service.
%
% Inputs:
%   - absPathToCacheMatFile
%     The absolute directory to the cache .mat file. We will load the data
%     there in if that .mat file exists; otherwise, we will create the file
%     and store the results into it.
%   - startXY, endXY
%     The UTM (x,y) locations for the start and end locations of the
%     profile.
%   - simConfigs
%     The simulation configs. We need here:
%       - simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M
%         The guaranteed spacial resolution for terrain profiles; a larger
%         value will decrease the simulation time but small obstacles may
%         get ingored.
%       - simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M
%         Similarly, the guaranteed spacial resolution for LiDAR profiles.
%       - simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE
%         The guaranteed minimum number of LiDAR z (or possibly elevation)
%         elements in one terrain profile; this will ensure non-empty
%         terrain profiles.
%       - simConfigs.utm2deg_speZone
%         The function to convert UTM (x, y) to GPS degrees (lat, lon).
%   - lidarFileXYCentroids
%     A matrix of the centroids (as rows) for the available LiDAR data
%     tiles in the output .mat files from preprocessIndianaLidarDataSet.m.
%   - lidarFileXYCoveragePolyshapes
%     A list of polyshapes for the available LiDAR data tiles in the output
%     .mat files from preprocessIndianaLidarDataSet.m.
%   - lidarMatFileAbsDirs
%     A list of the absolute directories to the output .mat files from
%     preprocessIndianaLidarDataSet.m.
%
% Outpus:
%   - terrainProfile, lidarProfile
%     The resulting terrain and LiDAR profiles.
%
% Yaguang Zhang, Purdue, 09/18/2019

if exist(absPathToCacheMatFile, 'file')
    cachedResults = load(absPathToCacheMatFile, ...
        'terrainProfile', 'lidarProfile');
    terrainProfile = cachedResults.terrainProfile;
    lidarProfile = cachedResults.lidarProfile;
else
    % Sample locations along the cellular tower and the drone for creating
    % terrain and LiDAR profiles.
    if simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M ...
            == simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M
        curTerrainProfileSampLocs = generateTerrainProfileSampLocs(...
            startXY, endXY, ...
            simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M,...
            simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE);
        [terrainProfile, lidarProfile, curEleForNanPts] ...
            = generateProfileSamps( ...
            curTerrainProfileSampLocs, simConfigs.utm2deg_speZone, ...
            lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
            lidarMatFileAbsDirs, 'both');
        
        boolsToReplaceWithEle = isnan(terrainProfile);
        terrainProfile(boolsToReplaceWithEle) ...
            = curEleForNanPts(boolsToReplaceWithEle);
        lidarProfile(boolsToReplaceWithEle) ...
            = curEleForNanPts(boolsToReplaceWithEle);
    else
        
        curTerrainProfileSampLocs = generateTerrainProfileSampLocs(...
            startXY, endXY, ...
            simConfigs.MAX_ALLOWED_TERRAIN_PROFILE_RESOLUATION_IN_M,...
            simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE);
        [terrainProfile, ~, curEleForNanTerrainPts] ...
            = generateProfileSamps( ...
            curTerrainProfileSampLocs, simConfigs.utm2deg_speZone, ...
            lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
            lidarMatFileAbsDirs, 'elevation');
        
        boolsToReplaceWithEle = isnan(terrainProfile);
        terrainProfile(boolsToReplaceWithEle) ...
            = curEleForNanTerrainPts(boolsToReplaceWithEle);
        
        curLidarProfileSampLocs = generateTerrainProfileSampLocs( ...
            startXY, endXY, ...
            simConfigs.MAX_ALLOWED_LIDAR_PROFILE_RESOLUATION_IN_M, ...
            simConfigs.MIN_NUM_OF_TERRAIN_SAMPLES_PER_PROFILE);
        [~, lidarProfile, curEleForNanLidarPts] ...
            = generateProfileSamps( ...
            curLidarProfileSampLocs, simConfigs.utm2deg_speZone, ...
            lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
            lidarMatFileAbsDirs, 'lidar');
        
        boolsToReplaceWithEle = isnan(lidarProfile);
        lidarProfile(boolsToReplaceWithEle) ...
            = curEleForNanLidarPts(boolsToReplaceWithEle);
    end
    
    save(absPathToCacheMatFile, ...
        'terrainProfile', 'lidarProfile');
end
end
% EOF