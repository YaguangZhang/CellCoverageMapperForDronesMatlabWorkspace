function [ tippBoundaryLatLons, tippBoundaryXYs, tippBoundaryUtmZone ] ...
    = loadTippBoundary(utmZone)
%LOADTIPPBOUNDARY A helper function to generate and load the Tippecanoe
%County (Tipp) boundary.
%
% Input:
%   - utmZone
%     Optional. Default to '16 T'. The UTM zone to use for converting GPS
%     (lat, lon) to UTM (x, y).
%
% Outputs:
%   - tippBoundaryLatLons
%     The boundary polygon for Tipp in terms of GPS (lats, lons) vertices.
%   - tippBoundaryXYs
%     The boundary polygon for Tipp in terms of UTM (x, y) vertices.
%   - tippBoundaryUtmZone
%     The UTM zone for tippBoundaryXYs.
%
% Yaguang Zhang, Purdue, 03/27/2022

% Parameters
ABS_PATH_TO_SHARED_FOLDER = evalin('base', 'ABS_PATH_TO_SHARED_FOLDER');
dirToCachedResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'StateInformation', 'IndianaBoundaries', 'Tipp');
fileToCachedResults = fullfile(dirToCachedResults, 'boundary.mat');

if exist(fileToCachedResults, 'file')
    cachedB = load(fileToCachedResults);
    tippBoundaryLatLons = cachedB.boundary.LAT_LON_BOUNDARY_OF_INTEREST;
    tippBoundaryXYs = cachedB.boundary.UTM_X_Y_BOUNDARY_OF_INTEREST;
    tippBoundaryUtmZone = cachedB.boundary.UTM_ZONE;

    if exist('utmZone', 'var')
        assert(strcmp(utmZone, ...
            tippBoundaryUtmZone), ['UTM zone mismatch! ', ...
            'Please delete the cached Tipp boundary to regenerate it!']);
    end
else
    if ~exist('utmZone', 'var')
        utmZone = '16 T';
    end

    STATEFP_IN = 18;
    COUNTY_NAMES_TIPP = {'Tippecanoe'};
    dirToUsCountyShp = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'StateInformation', 'tl_2021_us_county', 'tl_2021_us_county.shp');

    % Load the location information.
    counties = shaperead(dirToUsCountyShp);
    countiesIn = counties(str2double({counties.STATEFP})==STATEFP_IN);
    countiesTipp = countiesIn(ismember({countiesIn.NAME}, ...
        COUNTY_NAMES_TIPP));

    lonLatPolyshapeTipp = polyshape();
    for idxCnt = 1:length(COUNTY_NAMES_TIPP)
        lonLatPolyshapeTipp = union(lonLatPolyshapeTipp, ...
            polyshape(countiesTipp(idxCnt).X, countiesTipp(idxCnt).Y));
    end
    lonLatPolyshapeTipp = rmholes(lonLatPolyshapeTipp);

    tippBoundaryLats = lonLatPolyshapeTipp.Vertices(:,2);
    tippBoundaryLons = lonLatPolyshapeTipp.Vertices(:,1);

    % Close the polygon.
    tippBoundaryLats(end+1) = tippBoundaryLats(1);
    tippBoundaryLons(end+1) = tippBoundaryLons(1);

    tippBoundaryLatLons = [tippBoundaryLats, tippBoundaryLons];
    tippBoundaryUtmZone = utmZone;

    % For GPS and UTM conversions.
    [deg2utm_speZone, ~] ...
        = genUtmConvertersForFixedZone(tippBoundaryUtmZone);

    [tippBoundaryXs, tippBoundaryYs] = deg2utm_speZone( ...
        tippBoundaryLats, tippBoundaryLons);
    tippBoundaryXYs = [tippBoundaryXs, tippBoundaryYs];

    % Cache the restuls.
    boundary.LAT_LON_BOUNDARY_OF_INTEREST = tippBoundaryLatLons;
    boundary.UTM_X_Y_BOUNDARY_OF_INTEREST = tippBoundaryXYs;
    boundary.UTM_ZONE = tippBoundaryUtmZone;

    % Create directories if necessary.
    if exist(dirToCachedResults, 'dir')~=7
        mkdir(dirToCachedResults);
    end
    save(fileToCachedResults, 'boundary');
end

end
% EOF