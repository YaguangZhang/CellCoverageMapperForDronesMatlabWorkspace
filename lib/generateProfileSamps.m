function [ profileTerrain, profileLidar, eleForNanLocs ] ...
    = generateProfileSamps( ...
    profileSampLocs, utmToDegFct, ...
    lidarFileXYCentroids, lidarFileXYCoveragePolyshapes, ...
    lidarMatFileAbsDirs, terrainDataType, lidarRasterResolutionInM )
%GENERATEPROFILESAMPS Generate the terrain profile sample points (elevation
%or LiDAR).
%
% Inputs:
%   - profileSampLocs
%     A matrix of the UTM (x,y) coordinates (as rows) for the locations in
%     the terrian profile to fetch data.
%   - utmToDegFct
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
%   - lidarRasterResolutionInM
%     Optional. The LiDAR data raster resolution in meters. We will use
%     nearest neighbor interpolation for points close enough to a data tile
%     (distance to the neareast raster point <= sqrt(2)/2 x raster
%     resolution).
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

if ~exist('lidarRasterResolutionInM', 'var')
    % The Indiana rater LiDAR data has a resolution of 5 feet (~1.524 m).
    lidarRasterResolutionInM = 1.525;
end

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
    % Only load the fields needed to save time.
    fieldsToLoad = {'getEleFromXYFct', 'getLiDarZFromXYFct', ...
        'lidarEles', 'lidarXYZ', 'xYBoundryPolygon'};
    
    try
        load(lidarMatFileAbsDirs{curIdxTile}, fieldsToLoad{:});
    catch err
        disp(' ');
        warning('Error loading .mat LiDAR file!');
        dispErr(err);
        disp(' ');
        disp('Trying to fix the file...');
        validateLidarMatFile(lidarMatFileAbsDirs{curIdxTile}, simConfigs);
        load(lidarMatFileAbsDirs{curIdxTile}, fieldsToLoad{:});
        disp('succeeded!');
        disp(' ');
    end
    
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
    
    % Find the profile locations which has the closest data tile as the
    % current one and are close enough to the current data tile by
    % extending the current tile boundary.
    if length(xYBoundryPolygon.Vertices(:,1))==4
        % Most of the time, the boundary for the tile should have exactly 4
        % vertices. In this case, we will treat it as a rectangle and
        % simply shift its vertices.
        xYBoundryPolygonExt = extendRectPolyShape(...
            xYBoundryPolygon, lidarRasterResolutionInM);
    else
        % Otherwise, we will extend the boundary more precisely.
        xYBoundryPolygonExt = polybuffer( ...
            xYBoundryPolygon, ...
            sqrt(2)/2*lidarRasterResolutionInM);
        assert(xYBoundryPolygonExt.NumRegions==1, ...
            'The extended data tile should have only one region!');
    end
    
    % We only need to consider points not in the current tile.
    boolsCandidateProfileSampsNearCurTile = ...
        boolsProfileSampsWithCurTileAsClosest ...
        &(~boolsProfileSampsInCurTile);
    boolsProfileSampsNearCurTile = ...
        boolsCandidateProfileSampsNearCurTile;
    if any(boolsProfileSampsNearCurTile)
        try
            boolsProfileSampsNearCurTile( ...
                boolsCandidateProfileSampsNearCurTile) = ...
                isinterior(...
                xYBoundryPolygonExt, ...
                profileSampLocs(boolsCandidateProfileSampsNearCurTile,1), ...
                profileSampLocs(boolsCandidateProfileSampsNearCurTile,2));
        catch err
            parsave('./deleteme.mat', ...
                xYBoundryPolygonExt, profileSampLocs, ...
                boolsCandidateProfileSampsNearCurTile);
            throw(err);
        end
    end
    
    % Find the nearest neighbors.
    indicesNearestPtInCurTile = dsearchn(lidarXYZ(:,1:2), ...
        profileSampLocs(boolsProfileSampsNearCurTile,1:2));
    switch lower(terrainDataType)
        case 'elevation'
            profileTerrain(boolsProfileSampsInCurTile) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
            profileTerrain(boolsProfileSampsNearCurTile) ...
                = lidarEles(indicesNearestPtInCurTile);
        case 'lidar'
            profileLidar(boolsProfileSampsInCurTile) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
            profileLidar(boolsProfileSampsNearCurTile) ...
                = lidarXYZ(indicesNearestPtInCurTile, 3);
        case 'both'
            profileTerrain(boolsProfileSampsInCurTile) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
            profileTerrain(boolsProfileSampsNearCurTile) ...
                = lidarEles(indicesNearestPtInCurTile);
            
            profileLidar(boolsProfileSampsInCurTile) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsProfileSampsInCurTile,1), ...
                profileSampLocs(boolsProfileSampsInCurTile,2));
            profileLidar(boolsProfileSampsNearCurTile) ...
                = lidarXYZ(indicesNearestPtInCurTile, 3);
        otherwise
            error('Unsupported output data type!')
    end
end

boolsNeedFallbackEle = isnan(profileTerrain)&isnan(profileLidar);
[profileNanSampLats, profileNanSampLons] = utmToDegFct( ...
    profileSampLocs(boolsNeedFallbackEle,1), ...
    profileSampLocs(boolsNeedFallbackEle,2));
if any(boolsNeedFallbackEle)
    eleForNanLocs(boolsNeedFallbackEle) ...
        = queryElevationPointsFromUsgsInChunks(profileNanSampLats, ...
        profileNanSampLons);
end

end
% EOF