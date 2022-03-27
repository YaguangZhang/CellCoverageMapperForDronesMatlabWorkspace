function [ whinBoundaryLatLons, whinBoundaryXYs, whinBoundaryUtmZone ] ...
    = loadWhinBoundary(utmZone)
%LOADWHINBOUNDARY A helper function to generate and load the WHIN boundary.
%
% Input:
%   - utmZone
%     Optional. Default to '16 T'. The UTM zone to use for converting GPS
%     (lat, lon) to UTM (x, y).
%
% Outputs:
%   - whinBoundaryLatLons
%     The boundary polygon for WHIN in terms of GPS (lats, lons) vertices.
%   - whinBoundaryXYs
%     The boundary polygon for WHIN in terms of UTM (x, y) vertices.
%   - whinBoundaryUtmZone
%     The UTM zone for whinBoundaryXYs.
%
% Yaguang Zhang, Purdue, 11/23/2021

% Parameters
ABS_PATH_TO_SHARED_FOLDER = evalin('base', 'ABS_PATH_TO_SHARED_FOLDER');
dirToCachedResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries', 'WHIN');
fileToCachedResults = fullfile(dirToCachedResults, 'boundary.mat');

if exist(fileToCachedResults, 'file')
    cachedB = load(fileToCachedResults);
    whinBoundaryLatLons = cachedB.boundary.LAT_LON_BOUNDARY_OF_INTEREST;
    whinBoundaryXYs = cachedB.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    whinBoundaryUtmZone = cachedB.boundary.UTM_ZONE;

    if exist('utmZone', 'var')
        assert(strcmp(utmZone, ...
            whinBoundaryUtmZone), ['UTM zone mismatch! ', ...
            'Please delete the cached WHIN boundary to regenerate it!']);
    end
else
    if ~exist('utmZone', 'var')
        utmZone = '16 T';
    end

    STATEFP_IN = 18;
    COUNTY_NAMES_WHIN = {'Pulaski', 'White', 'Cass', 'Benton', ...
        'Carroll', 'Tippecanoe', 'Warren', 'Fountain', ...
        'Montgomery', 'Clinton'};
    dirToUsCountyShp = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'StateInformation', 'tl_2021_us_county', 'tl_2021_us_county.shp');

    % Load the location information.
    counties = shaperead(dirToUsCountyShp);
    countiesIn = counties(str2double({counties.STATEFP})==STATEFP_IN);
    countiesWhin = countiesIn(ismember({countiesIn.NAME}, ...
        COUNTY_NAMES_WHIN));

    lonLatPolyshapeWhin = polyshape();
    for idxCnt = 1:length(COUNTY_NAMES_WHIN)
        lonLatPolyshapeWhin = union(lonLatPolyshapeWhin, ...
            polyshape(countiesWhin(idxCnt).X, countiesWhin(idxCnt).Y));
    end
    lonLatPolyshapeWhin = rmholes(lonLatPolyshapeWhin);

    whinBoundaryLats = lonLatPolyshapeWhin.Vertices(:,2);
    whinBoundaryLons = lonLatPolyshapeWhin.Vertices(:,1);

    % Close the polygon.
    whinBoundaryLats(end+1) = whinBoundaryLats(1);
    whinBoundaryLons(end+1) = whinBoundaryLons(1);

    whinBoundaryLatLons = [whinBoundaryLats, whinBoundaryLons];
    whinBoundaryUtmZone = utmZone;

    % For GPS and UTM conversions.
    [deg2utm_speZone, ~] ...
        = genUtmConvertersForFixedZone(whinBoundaryUtmZone);

    [whinBoundaryXs, whinBoundaryYs] = deg2utm_speZone( ...
        whinBoundaryLats, whinBoundaryLons);
    whinBoundaryXYs = [whinBoundaryXs, whinBoundaryYs];

    % Cache the restuls.
    boundary.LAT_LON_BOUNDARY_OF_INTEREST = whinBoundaryLatLons;
    boundary.UTM_X_Y_BOUNDARY_OF_INTEREST = whinBoundaryXYs;
    boundary.UTM_ZONE = whinBoundaryUtmZone;

    % Create directories if necessary.
    if exist(dirToCachedResults, 'dir')~=7
        mkdir(dirToCachedResults);
    end
    save(fileToCachedResults, 'boundary');
end

end
% EOF