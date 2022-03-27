function [xYBoundryPolygon, lonLatBoundryPolygon, warnMsg] ...
    = loadCachedLidarDsmTile(curFullPathToSaveLidarResults)
%LOADCACHEDLIDARDSMTILE A helper snipet to try loading cached .mat LiDAR
%DSM tile.
%
% Input:
%   - curFullPathToSaveLidarResults
%     The absolute path to the .mat cache file.
%
% Outputs:
%   - xYBoundryPolygon, lonLatBoundryPolygon, warnMsg
%
% Yaguang Zhang, Purdue, 03/27/2022

% Clear last warning message.
lastwarn('');

% Make sure the file can be loaded properly.
try
    historyResult = load(curFullPathToSaveLidarResults, ...
        'xYBoundryPolygon', 'lonLatBoundryPolygon', ...
        'getLiDarZFromXYFct');
    xYBoundryPolygon = historyResult.xYBoundryPolygon;
    lonLatBoundryPolygon ...
        = historyResult.lonLatBoundryPolygon;

    % Make sure the function getLiDarZFromXYFct in the file has a valid
    % workspace. More specifically, the function handle
    % fctLonLatToLidarStatePlaneXY should be struct with field
    % 'spatialRef'.
    getLiDarZFromXYFctDetails = functions( ...
        historyResult.getLiDarZFromXYFct);
    fctLonLatToLidarStatePlaneXYDetails ...
        = functions(getLiDarZFromXYFctDetails.workspace{1} ...
        .fctLonLatToLidarStatePlaneXY);
    if ~isfield( ...
            fctLonLatToLidarStatePlaneXYDetails.workspace{1}, ...
            'spatialRef')
        warning('Field spatialRef not found!')
    end
catch err
    disp('            There was an error!')
    dispErr(err);

    % Output empty results.
    [xYBoundryPolygon, lonLatBoundryPolygon] = deal(polyshape);
    warning('The history result .mat file is invalid!');
end

% Check whether there is any warning in loading the desired data.
[warnMsg, ~] = lastwarn;

end
% EOF