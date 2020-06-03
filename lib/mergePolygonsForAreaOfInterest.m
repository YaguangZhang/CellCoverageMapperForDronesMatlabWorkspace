function [mergedPolygon] ...
    = mergePolygonsForAreaOfInterest(polygons, shrinkFactor)
%MERGEPOLYGONFORAREAOFINTEREST Compose one boundary polygon (polyshape) for
%the area covered by the input polygons (polyshapes).
%
% Input:
%   - polygons
%     A cell with all the 2D polygons to be considered.
%   - shrinkFactor
%     Optional. Shrink factor, specified as a scalar in the range of [0,1].
%     0 corresponds to the convex hull of the points. 1 corresponds to the
%     tightest single-region boundary around the points. Default to 0.5.
%
% Output:
%   - mergedPolygon
%     The resulting polygon from merging all input polygons.
%
% Yaguang Zhang, Purdue, 09/03/2019

if ~exist('shrinkFactor', 'var')
    shrinkFactor = 0.5;
end

inputPts = vertcat(cell2mat(cellfun(@(p) p.Vertices, polygons, ...
    'UniformOutput', false)));
boolsRowsToDiscard = isnan(inputPts(:,1)+inputPts(:,2));
inputPts(boolsRowsToDiscard, :) = [];
mergedPolygonVs = inputPts(boundary(inputPts, shrinkFactor), :);

% Remove duplicate vertices.
warning('off', 'MATLAB:polyshape:repairedBySimplify');
mergedPolygon = polyshape(mergedPolygonVs);
warning('on', 'MATLAB:polyshape:repairedBySimplify');

end
% EOF