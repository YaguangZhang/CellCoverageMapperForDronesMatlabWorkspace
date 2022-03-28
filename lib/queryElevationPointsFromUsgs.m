function [ eles ] = queryElevationPointsFromUsgs(lats, lons, language)
%QUERYELEVATIONPOINTSFROMUSGS Querry elevation data for input locations via
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
% Test:
%   elesTest = queryElevationPointsFromUsgs(40:50, -80:-70);
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

% Update 20210415: USGS supports values with at most 10 digits after the
% decimal point now. For example:
%   https://nationalmap.gov/epqs/pqs.php?x=-87.58214445193&y=40.16903761691&units=Meters&output=json
% does not work (11 digits), but
%   https://nationalmap.gov/epqs/pqs.php?x=-87.5821444519&y=40.1690376169&units=Meters&output=json
% works (10 digits).
%
% Update 20220328: USGS uses a string to indicate unavailable data points
% now.
numToStrFormatter = '%.10f';
maxNumOfTrials = 100;
urlUsgsEpqs = 'https://nationalmap.gov/epqs/pqs.php';
usgsNanValue = '-1000000';

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
                    if ~ischar(fetchedEle)
                        eles(idxPt) = fetchedEle;
                    else
                        assert(strcmp(fetchedEle, usgsNanValue), ...
                            ['Unknown USGS ele value: ', fetchedEle, '!']);
                    end
                    ctnTrials = inf;
                catch err
                    disp(['There was an error fetching data ', ...
                        'from USGS via Matlab webread (Pt #', ...
                        num2str(idxPt), '/', num2str(numQueryPts), ...
                        ' Trial #', num2str(ctnTrials), '/', ...
                        num2str(maxNumOfTrials), ')!'])
                    dispErr(err);
                    if ctnTrials == maxNumOfTrials
                        save('errorWorkspaceInQuerryEle.mat');
                        error('Run out of allowed number of trials!');
                    end
                end
            end
        end
    case 'python'
        ctnTrials = 0;
        while ctnTrials<maxNumOfTrials
            ctnTrials = ctnTrials+1;
            try
                boolsEleFetched = false(numPts, 1);
                % Concurrently fetch the elevation data from USGS via
                % Python with a maximum batch size.
                while(any(~boolsEleFetched))
                    curIdxStart = find(~boolsEleFetched, 1);
                    if (numPts-curIdxStart+1) ...
                            <=MAX_NUM_OF_CONCURRENT_REQUESTS
                        curIdxEnd = numPts;
                    else
                        curIdxEnd ...
                            = curIdxStart+MAX_NUM_OF_CONCURRENT_REQUESTS-1;
                    end
                    curUrls = urls(curIdxStart:curIdxEnd);
                    fetchedPyJsons = py.ConcurrentWebreader ...
                        .concurrent_fetch_urls( ...
                        py.list(curUrls));
                    fetchedJsonsCell = cell(fetchedPyJsons);

                    boolsIsValidJson = cellfun( ...
                        @(json) ~isa(json, 'py.NoneType'), ...
                        fetchedJsonsCell);

                    try
                        fetchedJsons = cellfun(@(jo) ...
                            jsondecode(native2unicode(jo)), ...
                            fetchedJsonsCell(boolsIsValidJson));
                    catch err
                        % Save the fetched JSON object to a file for
                        % debugging.
                        fetchedPyJsonsChar = char(fetchedPyJsons);
                        warning( ...
                            ['Unable to process fetchedPyJsons: ', ...
                            fetchedPyJsonsChar]);
                        parsave('./degbugFetchedPyJsonsChar', ...
                            curUrls, fetchedPyJsonsChar);
                        throw(err);
                    end

                    indicesToUpdate = curIdxStart:curIdxEnd;
                    indicesToUpdate = indicesToUpdate(boolsIsValidJson);
                    fetchedElesCell = arrayfun(@(resp) resp ...
                        .USGS_Elevation_Point_Query_Service ...
                        .Elevation_Query.Elevation, fetchedJsons', ...
                        'UniformOutput', false);
                    boolsUsgsEleIsChar = cellfun( ...
                        @ (ele) ischar(ele), fetchedElesCell);
                    cellfun(@ (charEle) ...
                        assert(strcmp(charEle, usgsNanValue), ...
                        ['Unknown USGS ele value: ', charEle, '!']), ...
                        fetchedElesCell(boolsUsgsEleIsChar));

                    fetchedElesCell(boolsUsgsEleIsChar) = {nan};
                    eles(indicesToUpdate) = cell2mat(fetchedElesCell);
                    boolsEleFetched(indicesToUpdate) = true;
                end

                ctnTrials = inf;
            catch err
                disp(['There was an error fetching data ', ...
                    'from USGS via Python (In total ', ...
                    num2str(numQueryPts), ' points', ...
                    ' Trial #', num2str(ctnTrials), '/', ...
                    num2str(maxNumOfTrials), ')!'])
                dispErr(err);
                if ctnTrials == maxNumOfTrials
                    save('errorWorkspaceInQuerryEle.mat');
                    error('Run out of allowed number of trials!');
                end
            end
        end
    otherwise
        error(['Unsupported language ', language, ...
            ' for fetching USGS elevation data!'])
end

end
% EOF