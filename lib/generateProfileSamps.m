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

% [indicesClosestTile, ~] ...
%     = dsearchn(lidarFileXYCentroids, profileSampLocs);
indicesClosestTile = knnsearch(lidarFileXYCentroids, profileSampLocs);

[numOfProfileSampLocs,~] = size(profileSampLocs);
[profileTerrain, profileLidar, eleForNanLocs] ...
    = deal(nan(numOfProfileSampLocs,1));

% boolsProfileSampsInClosestTile = arrayfun( ...
%     @(idx) ...
%      isinterior(...
%     lidarFileXYCoveragePolyshapes{indicesClosestTile(idx)}, ...
%      profileSampLocs(idx,1), profileSampLocs(idx,2)), ...
%     1:numOfProfileSampLocs);

deltaDist = lidarRasterResolutionInM*2;
for curIdxTile = unique(indicesClosestTile)'
    % Only load the needed fields to save time.
    fieldsToLoad = {'getEleFromXYFct', 'getLiDarZFromXYFct', ...
        'lidarEles', 'lidarXs', 'lidarYs', 'lidarZs'};

    try
        load(lidarMatFileAbsDirs{curIdxTile}, fieldsToLoad{:});
    catch err
        disp(' ');
        warning('Error loading .mat LiDAR file!');
        dispErr(err);
        disp(' ');
        disp('Trying to fix the file...');
        validateLidarMatFile(lidarMatFileAbsDirs{curIdxTile});
        load(lidarMatFileAbsDirs{curIdxTile}, fieldsToLoad{:});
        disp('succeeded!');
        disp(' ');
    end

    % Find the profile locations which has the closest data tile as the
    % current one ...
    boolsProfileSampsWithCurTileAsClosest = indicesClosestTile==curIdxTile;
    % ... and are indeed in the polyshape of the current tile.
    xYBoundryPolygon = lidarFileXYCoveragePolyshapes{curIdxTile};
    boolsProfileSampsInCurTile = boolsProfileSampsWithCurTileAsClosest;
    boolsProfileSampsInCurTile( ...
        boolsProfileSampsWithCurTileAsClosest) = InPolygon( ...
        profileSampLocs(boolsProfileSampsWithCurTileAsClosest,1), ...
        profileSampLocs(boolsProfileSampsWithCurTileAsClosest,2), ...
        xYBoundryPolygon.Vertices(:,1), xYBoundryPolygon.Vertices(:,2));

    % Find the profile locations which has the closest data tile as the
    % current one ... and are NOT in the current data tile but close enough
    % to be included if we extend the current tile boundary a little.
    xYBoundryPolygonExt = polybuffer( ...
        xYBoundryPolygon, sqrt(2)/2*lidarRasterResolutionInM);
    assert(xYBoundryPolygonExt.NumRegions==1, ...
        'The extended data tile should have only one region!');

    % We only need to consider points not in the current tile.
    boolsCandidateProfileSampsNearCurTile = ...
        boolsProfileSampsWithCurTileAsClosest ...
        &(~boolsProfileSampsInCurTile);
    boolsProfileSampsNearCurTile = ...
        boolsCandidateProfileSampsNearCurTile;
    if any(boolsProfileSampsNearCurTile)
        try
            boolsProfileSampsNearCurTile( ...
                boolsCandidateProfileSampsNearCurTile) = InPolygon( ...
                profileSampLocs( ...
                boolsCandidateProfileSampsNearCurTile,1), ...
                profileSampLocs( ...
                boolsCandidateProfileSampsNearCurTile,2), ...
                xYBoundryPolygonExt.Vertices(:,1), ...
                xYBoundryPolygonExt.Vertices(:,2));
        catch err
            parsave('./debugInfoForLidarZExtrapolation.mat', ...
                xYBoundryPolygonExt, profileSampLocs, ...
                boolsCandidateProfileSampsNearCurTile);
            throw(err);
        end
    end

    % For these points' LiDAR z, use nearest neighbor extrapolation,
    % instead of linear interp2, because the latter will output NaN in this
    % case. Note that:
    %   indicesNearestPtInCurTile = dsearchn([lidarXs, lidarYs], ...
    %       profileSampLocs(boolsProfileSampsNearCurTile,1:2));
    % is slower than:
    %   indicesNearestPtInCurTile = knnsearch([lidarXs, lidarYs], ...
    %       profileSampLocs(boolsProfileSampsNearCurTile,1:2));
    % To speed things up, we will process the (lidarX, lidarY) pairs one by
    % one.
    numOfProfileSampsNearCurTile = sum(boolsProfileSampsNearCurTile);
    indicesProfileSampsNearCurTile = find(boolsProfileSampsNearCurTile);

    indicesNearestPtInCurTile = nan(numOfProfileSampsNearCurTile, 1);
    for idxSamp = 1:numOfProfileSampsNearCurTile
        curIdxProfileSampsNearCurTile ...
            = indicesProfileSampsNearCurTile(idxSamp);
        curProfileSampLoc ...
            = profileSampLocs(curIdxProfileSampsNearCurTile,1:2);
        % Filter out grid points that are too far away.
        boolsGridPtTooFar = abs(lidarXs-curProfileSampLoc(1))>deltaDist ...
            | abs(lidarYs-curProfileSampLoc(2))>deltaDist;
        curIndicesNearByGridPts = find(~boolsGridPtTooFar);

        curIndicesNearestPtInCurTile = knnsearch( ...
            [lidarXs(curIndicesNearByGridPts), ...
            lidarYs(curIndicesNearByGridPts)], ...
            curProfileSampLoc);
        indicesNearestPtInCurTile(idxSamp) ...
            = curIndicesNearByGridPts(curIndicesNearestPtInCurTile);
    end

    switch lower(terrainDataType)
        case 'elevation'
            % The terrain information saved should be able to cover all the
            % points.
            boolsEleProfSampsToUpdate ...
                = (boolsProfileSampsInCurTile ...
                | boolsProfileSampsNearCurTile) ...
                & isnan(profileTerrain);
            profileTerrain(boolsEleProfSampsToUpdate) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsEleProfSampsToUpdate,1), ...
                profileSampLocs(boolsEleProfSampsToUpdate,2));
        case 'lidar'
            boolsLidarProfSampsToUpdate ...
                = boolsProfileSampsInCurTile & isnan(profileLidar);
            profileLidar(boolsLidarProfSampsToUpdate) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsLidarProfSampsToUpdate,1), ...
                profileSampLocs(boolsLidarProfSampsToUpdate,2));
            candidateNearTileLidarZs = lidarZs(indicesNearestPtInCurTile);
            profileLidar( ...
                boolsProfileSampsNearCurTile & isnan(profileLidar)) ...
                = candidateNearTileLidarZs( ...
                isnan(profileLidar(boolsProfileSampsNearCurTile)));
        case 'both'
            boolsEleProfSampsToUpdate ...
                = (boolsProfileSampsInCurTile ...
                | boolsProfileSampsNearCurTile) ...
                & isnan(profileTerrain);
            profileTerrain(boolsEleProfSampsToUpdate) ...
                = getEleFromXYFct( ...
                profileSampLocs(boolsEleProfSampsToUpdate,1), ...
                profileSampLocs(boolsEleProfSampsToUpdate,2));

            boolsLidarProfSampsToUpdate ...
                = boolsProfileSampsInCurTile & isnan(profileLidar);
            profileLidar(boolsLidarProfSampsToUpdate) ...
                = getLiDarZFromXYFct( ...
                profileSampLocs(boolsLidarProfSampsToUpdate,1), ...
                profileSampLocs(boolsLidarProfSampsToUpdate,2));
            candidateNearTileLidarZs = lidarZs(indicesNearestPtInCurTile);
            profileLidar( ...
                boolsProfileSampsNearCurTile & isnan(profileLidar)) ...
                = candidateNearTileLidarZs( ...
                isnan(profileLidar(boolsProfileSampsNearCurTile)));
        otherwise
            error('Unsupported output data type!')
    end
end

switch lower(terrainDataType)
    case 'elevation'
        boolsNeedFallbackEle = isnan(profileTerrain);
    case 'lidar'
        boolsNeedFallbackEle = isnan(profileLidar);
    case 'both'
        boolsNeedFallbackEle = isnan(profileTerrain) | isnan(profileLidar);
end

if any(boolsNeedFallbackEle)
    [profileNanSampLats, profileNanSampLons] = utmToDegFct( ...
        profileSampLocs(boolsNeedFallbackEle,1), ...
        profileSampLocs(boolsNeedFallbackEle,2));

    eleForNanLocs(boolsNeedFallbackEle) ...
        = queryElevationPointsFromUsgsInChunks(profileNanSampLats, ...
        profileNanSampLons);
end

end
% EOF