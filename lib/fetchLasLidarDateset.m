function [lidarXYZ, hFigGoogleMap, hFigLidarData] ...
    = fetchLasLidarDateset(absPathToLas, absPathToLatLonRange, ...
    gridResolutionInM, figVisible)
%FETCHLASLIDARDATESET Load LiDAR data from the specified .las file and
%generate some illustration figures for it.
% Inputs:
%   - absPathToLas
%     The absolute path to the .las/.laz file of the LiDAR data.
%   - absPathToLatLonRange
%     The absolute path to the log file (.txt or .csv) of the (lat, lon)
%     range of the LiDAR dataset. This is necessary to properly determine
%     the UTM zone label and generate the Google map illustration figure.
%     The file should store the bottom-left and top-right boundary points
%     with the first column being lat and the second column being lon.
%   - gridResolutionInM
%     An optional scalar to specify the spatial resolution for the LiDAR
%     data. If specified, a grid will be constructed accordingly and linear
%     interpolation will be used to fetch LiDAR data on the grid.
%   - figVisible
%     A illustration plot for the resulting ground roughness will also be
%     generated. This flag constrols whether this plot will be visible.
%
% Outpus:
%   - lidarXYZ
%     A matrix for the fetched LiDAR data with each row being the (UTM x,
%     UTM y, LiDAR z) for a LiDAR sample.
%   - hFigGoogleMap, hFigLidarData
%     Optional outputs for the illustration figure handles. Note that, for
%     the LiDAR illustration figure, hFigLidarData, we will need
%     gridResolutionInM to be properly set so that a grid will be generated
%     to be used in the surf command.
%
% Yaguang Zhang, Purdue, 11/21/2019

if ~exist('gridResolutionInM', 'var')
    gridResolutionInM = nan;
end

% Cache results for speed.
[parentDirPathForLas, filenameForLas] = fileparts(absPathToLas);
absPathToCacheResults = fullfile(parentDirPathForLas, ...
    [filenameForLas, '_gridResolutionInM_', ...
    strrep(num2str(gridResolutionInM), '.', '_'), '.mat']);

if exist(absPathToCacheResults, 'file')
    historyResults = load(absPathToCacheResults);
    lidarXYZ = historyResults.lidarXYZ;
    utmZone = historyResults.utmZone;
else
    % Load LiDAR data boundary points.
    latLonRange = csvread(absPathToLatLonRange);
    minLat = min(latLonRange(:,1));
    maxLat = max(latLonRange(:,1));
    minLon = min(latLonRange(:,2));
    maxLon = max(latLonRange(:,2));
    
    [minX, minY, utmZone] = deg2utm(minLat, minLon);
    [maxX, maxY, utmZoneMax] = deg2utm(maxLat, maxLon);
    
    assert(strcmp(utmZone, utmZoneMax), ...
        'UTM zone labels are not consistent for the boundary points!');
    
    % Load LiDAR data.
    lidarData = lasdata(absPathToLas);
    lidarXYZ = [lidarData.x, lidarData.y, lidarData.z];
    
    if ~isnan(gridResolutionInM)
        gridXLabels = constructAxisGrid( ...
            mean([maxX, minX]), ...
            floor((maxX-minX)./gridResolutionInM), gridResolutionInM);
        gridYLabels = constructAxisGrid( ...
            mean([maxY, minY]), ...
            floor((maxY-minY)./gridResolutionInM), gridResolutionInM);
        [gridXs, gridYs] = meshgrid(gridXLabels,gridYLabels);
        
        getLiDarZFromXYFct ...
            = scatteredInterpolant(lidarData.x, lidarData.y, lidarData.z);
        gridZs = getLiDarZFromXYFct(gridXs, gridYs);
        lidarXYZ = [gridXs(:), gridYs(:), gridZs(:)];
    end
    
    save(absPathToCacheResults, 'lidarXYZ', 'utmZone');
end

if nargout>1
    % It is necessary to plot the Google map illustration. We first
    % construct a polygon for the area of interest.
    indicesLidarAreaConvHull = convhull(lidarXYZ(:,1), lidarXYZ(:,2));
    UTM_X_Y_BOUNDARY_OF_INTEREST ...
        = [lidarXYZ(indicesLidarAreaConvHull,1), ...
        lidarXYZ(indicesLidarAreaConvHull,2)];
    [numOfBoIPts, ~] = size(UTM_X_Y_BOUNDARY_OF_INTEREST);
    
    % Convert the UTM coordinates to GPS locations.
    [latsBoI, lonsBoI] = utm2deg( ...
        UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
        UTM_X_Y_BOUNDARY_OF_INTEREST(:,2), ...
        repmat(utmZone, numOfBoIPts, 1));
    
    hFigGoogleMap = figure('Visible', figVisible);
    h_BoI = plot(lonsBoI, latsBoI, '--w', 'LineWidth', 1.5);
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)')
    grid on; view(2); 
    plot_google_map('MapType', 'satellite');
    makescale;
    legend(h_BoI, 'Area of Interest');
    
    if nargout>2 && ~isnan(gridResolutionInM)
        % It is necessary to plot the LiDAR data illustration.
        xLabelToSet = 'UTM x (m)';
        yLabelToSet = 'UTM y (m)';
        cLabelToSet = 'LiDAR z (m)';
        
        hFigLidarData = figure('Visible', figVisible);
        plot3k(lidarXYZ, 'Labels', {'', xLabelToSet, yLabelToSet, '', ...
            cLabelToSet});
        view(2); axis equal; axis tight;
        xlabel('UTM x (m)'); ylabel('UTM y (m)');
    end
end

end
% EOF