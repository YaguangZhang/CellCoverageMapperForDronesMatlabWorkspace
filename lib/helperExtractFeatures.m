function [normals,curvatures,neighInds] = helperExtractFeatures(ptCloud,numNeighbors)
%helperExtractFeatures Computes features from point cloud using normals
%
%   This is an example helper function that is subject to change or removal
%   in future releases.

% Copyright 2020 MathWorks, Inc.
normals = pcnormals(ptCloud, numNeighbors);

% Estimate neigbor indices and curvature for each point
neighInds = nan(numNeighbors, ptCloud.Count);
curvatures = nan(ptCloud.Count, 1);
locs = ptCloud.Location;
for i=1:ptCloud.Count
    if ~any(isnan(locs(i,:)))
        ind = findNearestNeighbors(ptCloud, ptCloud.Location(i,:),numNeighbors);
        % Store neigbor indices
        neighInds(:,i) = ind;
        eigVal = eig(cov(ptCloud.Location(ind,:)));
        % Store curvature
        curvatures(i) = min(eigVal) / sum(eigVal);
    end
end
end