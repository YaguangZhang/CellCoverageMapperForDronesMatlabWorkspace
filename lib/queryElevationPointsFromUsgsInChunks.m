function [ eles ] = queryElevationPointsFromUsgsInChunks( ...
    lats, lons, language, maxChunkSize)
%QUERYELEVATIONPOINTSFROMUSGSINCHUNKS Querry elevation data for input
%locations via the National Map - Elevation Point Query Service.
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
%   - maxChunkSize
%     The maximum number of locations to querry at once through the
%     function queryElevationPointsFromUsgs.
% Output:
%   - The corresponding elevation in meters.
%
% Yaguang Zhang, Purdue, 06/01/2020

if ~exist('language', 'var')
    language = 'Python';
end

if ~exist('maxChunkSize', 'var')
    maxChunkSize = 50;
end

numQueryPts = length(lats);
eles = nan(size(lats));
for idxChunkStart = 1:maxChunkSize:numQueryPts
    idxChunkEnd = idxChunkStart+maxChunkSize-1;
    if idxChunkEnd>numQueryPts
        idxChunkEnd = numQueryPts;
    end
    eles(idxChunkStart:idxChunkEnd) = queryElevationPointsFromUsgs( ...
        lats(idxChunkStart:idxChunkEnd), ...
        lons(idxChunkStart:idxChunkEnd), ...
        language);
end

end
% EOF