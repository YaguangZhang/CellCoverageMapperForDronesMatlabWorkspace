function [ profileTerrain, profileLidar, eleForNanLocs ] ...
    = generateProfileSamps( ...
    profileSampLocs, degToUtmFct, ...
    lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
    lidarMatFileAbsDirs, terrainDataType)
%GENERATEPROFILESAMPS Generate the terrain profile sample points (elevation
%or LiDAR).
%
% Inputs:
%   - profileSampLocs
%     A matrix of the UTM (x,y) coordinates (as rows) for the locations in
%     the terrian profile to fetch data.
%   - degToUtmFct
%     The funtion to convert UTM (x,y) to GPS (lat, lon).
%   - lidarFileXYCentroids
%     A matrix of the centroids (as rows) for the available LiDAR data
%     tiles in the output .mat files from preprocessIndianaLidarDataSet.m.
%   - lidarFileXYCoveragePolyshapes
%     A list of polyshapes for the available LiDAR data tiles in the output
%     .mat files from preprocessIndianaLidarDataSet.m.
%   - lidarMatFileAbsDirs
%     A list of the absolute directories to the output .mat files from
%     preprocessIndianaLidarDataSet.m.
%   - terrainDataType
%     A string controling what data to fetch: 'elevation', 'lidar', or
%     'both'.
%
% Outputs:
%   - profileTerrain, profileLidar
%     A column vector of the specified data (according to the available
%     LiDAR/elevation data tiles in the output .mat files from
%     preprocessIndianaLidarDataSet.m.) for the input profile sample
%     locations. If any profile sample point is outside of the area covered
%     by the available data files, an NaN element will be stored in the
%     output profile at the corresponding loction.
%   - eleForNanLocs
%     For NaN elements in the output profiles, we will fetch elevation data
%     from USGS as a fallback and store them in eleForNanLocs.
%     eleForNanLocs is of the same size as the output profiles, and has
%     valid elevation data for elements cooresponding to NaNs in the
%     profiles.
%
% Yaguang Zhang, Purdue, 09/17/2019

[indicesClosestTile, ~] = dsearchn(lidarFileXYCentroids, profileSampLocs);

[numOfProfileSampLocs,~] = size(profileSampLocs);
[profileTerrain, profileLidar, eleForNanLocs] ...
    = deal(nan(numOfProfileSampLocs,1));

% boolsProfileSampsInClosestTile = arrayfun( ...
%     @(idx) ...
%      isinterior(...
%     lidarFileXYCoveragePolyshapes{indicesClosestTile(idx)}, ...
%      profileSampLocs(idx,1), profileSampLocs(idx,2)), ...
%     1:numOfProfileSampLocs);

for curIdxTile = unique(indicesClosestTile)'
    load(lidarMatFileAbsDirs{curIdxTile});
    
    % Find the profile locations which has the closest data tile as the
    % current one ...
    boolsProfileSampsWithCurTileAsClosest = indicesClosestTile==curIdxTile;
    % ... and are indeed in the polyshape of the current tile
    boolsProfileSampsInCurTile = boolsProfileSampsWithCurTileAsClosest;
    boolsProfileSampsInCurTile( ...
        boolsProfileSampsWithCurTileAsClosest) = ...
        isinterior(...
        lidarFileXYCoveragePolyshapes{curIdxTile}, ...
        profileSampLocs(boolsProfileSampsWithCurTileAsClosest,1), ...
        profileSampLocs(boolsProfileSampsWithCurTileAsClosest,2));
    
    switch lower(terrainDataType)
        case 'elevation'
            profileTerrain(boolsProfileSampsInCurTile) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
        case 'lidar'
            profileLidar(boolsProfileSampsInCurTile) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
        case 'both'
            profileTerrain(boolsProfileSampsInCurTile) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
            profileLidar(boolsProfileSampsInCurTile) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
        otherwise
            error('Unsupported output data type!')
    end
end

boolsNeedFallbackEle = isnan(profileTerrain)&isnan(profileLidar);
[profileNanSampLats, profileNanSampLons] = degToUtmFct( ...
    profileSampLocs(boolsNeedFallbackEle,1), ...
    profileSampLocs(boolsNeedFallbackEle,2));
if any(boolsNeedFallbackEle)
    eleForNanLocs(boolsNeedFallbackEle) ...
        = queryElevationPointsFromUsgs(profileNanSampLats, ...
        profileNanSampLons);
end

end
% EOF