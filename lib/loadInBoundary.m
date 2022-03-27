function [ inBoundaryLatLons, inBoundaryXYs, inBoundaryUtmZone ] ...
    = loadInBoundary(utmZone)
%LOADINBOUNDARY A helper function to generate and load the IN boundary.
%
% Input:
%   - utmZone
%     Optional. Default to '16 T'. The UTM zone to use for converting GPS
%     (lat, lon) to UTM (x, y).
%
% Outputs:
%   - inBoundaryLatLons
%     The boundary polygon for IN in terms of GPS (lats, lons) vertices.
%   - inBoundaryXYs
%     The boundary polygon for IN in terms of UTM (x, y) vertices.
%   - inBoundaryUtmZone
%     The UTM zone for inBoundaryXYs.
%
% Yaguang Zhang, Purdue, 03/26/2022

% Parameters
ABS_PATH_TO_SHARED_FOLDER = evalin('base', 'ABS_PATH_TO_SHARED_FOLDER');
dirToCachedResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries', 'IN');
fileToCachedResults = fullfile(dirToCachedResults, 'boundary.mat');

if exist(fileToCachedResults, 'file')
    cachedB = load(fileToCachedResults);
    inBoundaryLatLons = cachedB.boundary.LAT_LON_BOUNDARY_OF_INTEREST;
    inBoundaryXYs = cachedB.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    inBoundaryUtmZone = cachedB.boundary.UTM_ZONE;

    if exist('utmZone', 'var')
        assert(strcmp(utmZone, ...
            inBoundaryUtmZone), ['UTM zone mismatch! ', ...
            'Please delete the cached IN boundary to regenerate it!']);
    end
else
    if ~exist('utmZone', 'var')
        utmZone = '16 T';
    end

    STATEFP_IN = 18;
    dirToUsCountyShp = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'StateInformation', 'tl_2021_us_county', 'tl_2021_us_county.shp');

    % Load the location information.
    counties = shaperead(dirToUsCountyShp);
    countiesIn = counties(str2double({counties.STATEFP})==STATEFP_IN);

    lonLatPolyshapeIn = polyshape();
    for idxCnt = 1:length(countiesIn)
        lonLatPolyshapeIn = union(lonLatPolyshapeIn, ...
            polyshape(countiesIn(idxCnt).X, countiesIn(idxCnt).Y));
    end
    lonLatPolyshapeIn = rmholes(lonLatPolyshapeIn);

    inBoundaryLats = lonLatPolyshapeIn.Vertices(:,2);
    inBoundaryLons = lonLatPolyshapeIn.Vertices(:,1);

    % Close the polygon.
    inBoundaryLats(end+1) = inBoundaryLats(1);
    inBoundaryLons(end+1) = inBoundaryLons(1);

    inBoundaryLatLons = [inBoundaryLats, inBoundaryLons];
    inBoundaryUtmZone = utmZone;

    % For GPS and UTM conversions.
    [deg2utm_speZone, ~] ...
        = genUtmConvertersForFixedZone(inBoundaryUtmZone);

    [inBoundaryXs, inBoundaryYs] = deg2utm_speZone( ...
        inBoundaryLats, inBoundaryLons);
    inBoundaryXYs = [inBoundaryXs, inBoundaryYs];

    % Cache the restuls.
    boundary.LAT_LON_BOUNDARY_OF_INTEREST = inBoundaryLatLons;
    boundary.UTM_X_Y_BOUNDARY_OF_INTEREST = inBoundaryXYs;
    boundary.UTM_ZONE = inBoundaryUtmZone;

    % Create directories if necessary.
    if exist(dirToCachedResults, 'dir')~=7
        mkdir(dirToCachedResults);
    end
    save(fileToCachedResults, 'boundary');
end

end
% EOF