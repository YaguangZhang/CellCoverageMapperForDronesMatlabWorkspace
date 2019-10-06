function [ eles ] = queryElevationPointsFromUsgs(lats, lons, language)
%QUERYELEVATIONPOINTSFROMUSGS Querry elevation data for input locations vis
%the National Map - Elevation Point Query Service.
%
% More information can be found at https://nationalmap.gov/epqs/
%
% Inputs:
%   - lats, lons
%     Column vectors for the GPS latitudes and longitudes of the locations
%     to query.
%   - language
%     Optional. Either 'Python' (default; faster) or 'Matlab'. One need to
%     make sure first Python works in the current Matlab instance, or
%     specify a valid Python path at the base workspace in the variable
%     ABS_PATH_TO_PYTHON (so that this function can automatically set up
%     Python).
% Output:
%   - The corresponding elevation in meters.
%
% Yaguang Zhang, Purdue, 09/18/2019

if ~exist('language', 'var')
    language = 'Python';
end

if strcmpi(language, 'Python')
    % Make sure Python is available.
    curPythonVersion = pyversion;
    if isempty(curPythonVersion)
        pyversion(evalin('base', 'ABS_PATH_TO_PYTHON'));
    end
    
    % Limit the number of concurrent HTTP requests to reduce the number of
    % failures.
    MAX_NUM_OF_CONCURRENT_REQUESTS = 100;
end

numToStrFormatter = '%.12f';
maxNumOfTrials = 10;
urlUsgsEpqs = 'https://nationalmap.gov/epqs/pqs.php';
usgsNanValue = -1000000;

numQueryPts = length(lats);
urls = arrayfun(@(idx) sprintf( ...
    [urlUsgsEpqs, '?x=%s&y=%s&units=Meters&output=json'], ...
    num2str(lons(idx), numToStrFormatter), ...
    num2str(lats(idx), numToStrFormatter)), 1:numQueryPts, ...
    'UniformOutput', false);

eles = nan(size(lats));
numPts = length(eles);
switch lower(language)
    case 'matlab'
        % Sequentially fetch the elevation data from USGS point by point.
        for idxPt = 1:numQueryPts
            ctnTrials = 0;
            while ctnTrials<maxNumOfTrials
                ctnTrials = ctnTrials+1;
                try
                    resp = webread(urls{idxPt});
                    fetchedEle = resp ...
                        .USGS_Elevation_Point_Query_Service ...
                        .Elevation_Query.Elevation;
                    % Only record valid elevation values.
                    if fetchedEle ~= usgsNanValue
                        eles(idxPt) = fetchedEle;
                    end
                    ctnTrials = inf;
                catch err
                    disp(['There was an error fetching data ', ...
                        'from USGS via Matlab webread (Pt #', ...
                        num2str(idxPt), '/', num2str(numQueryPts), ...
                        ' Trial #', num2str(ctnTrials), '/', ...
                        num2str(maxNumOfTrials), ')!'])
                    dispErr(err);
                end
            end
        end
    case 'python'
        ctnTrials = 0;
        while ctnTrials<maxNumOfTrials
            ctnTrials = ctnTrials+1;
            try
                % Concurrently fetch the elevation data from USGS via
                % Python with a maximum batch size.
                while(any(isnan(eles)))
                    curIdxStart = find(isnan(eles), 1);
                    if (numPts-curIdxStart+1) ...
                            <=MAX_NUM_OF_CONCURRENT_REQUESTS
                        curIdxEnd = numPts;
                    else
                        curIdxEnd ...
                            = curIdxStart+MAX_NUM_OF_CONCURRENT_REQUESTS-1;
                    end
                    fetchedPyJsons = py.ConcurrentWebreader ...
                        .concurrent_fetch_urls( ...
                        py.list(urls(curIdxStart:curIdxEnd)));
                    fetchedJsons = cellfun(@(jo) ...
                        jsondecode(native2unicode(jo)), ...
                        cell(fetchedPyJsons));
                    eles(curIdxStart:curIdxEnd) = arrayfun(@(resp) resp ...
                        .USGS_Elevation_Point_Query_Service ...
                        .Elevation_Query.Elevation, fetchedJsons');
                end
                
                eles(eles==usgsNanValue) = nan;
                ctnTrials = inf;
            catch err
                disp(['There was an error fetching data ', ...
                    'from USGS via Python (In total ', ...
                    num2str(numQueryPts), ' points', ...
                    ' Trial #', num2str(ctnTrials), '/', ...
                    num2str(maxNumOfTrials), ')!'])
                dispErr(err);
            end
        end
    otherwise
        error(['Unsupported language ', language, ...
            ' for fetching USGS elevation data!'])
end

end
% EOF