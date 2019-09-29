function [centroids] = extractCentroidsFrom2DPolyCell(polyCell)
%EXTRACTCENTROIDSFROM2DPOLYCELL Extract the centroids for 2D polygons
%stored in a cell.
%
% Input:
%   - polyCell  
%     A vector cell containing polygons.
%
% Output:
%	- centroids
%     A matrix containing all the centroids for polygons in polyCell, with
%     each row corresponds to one polygon.
%
% Yaguang Zhang, Purdue, 09/09/2019

numOfPolys = length(polyCell);

if numOfPolys == 0
    centroids = [];
else
    [~, numOfDimensions] = size(polyCell{1}.Vertices);
    assert(numOfDimensions==2, ...
        'The dimension of the polygon vertices should be 2!')
    
    centroids = nan(numOfPolys, numOfDimensions);
    for idxP = 1:numOfPolys
        [centroids(idxP, 1), centroids(idxP, 2)] ...
            = centroid(polyCell{idxP});
    end
end

end
% EOF