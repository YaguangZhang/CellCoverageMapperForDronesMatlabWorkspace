function labels = helperClassify(normals,curvatures,neighInds,norThres,curThres)
%helperClassify Classifies point clouds into high curvature and low
%   curvature segments
%
%     Class    |  Label id
%   Vegetation |   1 
%   Building   |   2
%
%   This is an example helper function that is subject to change or removal
%   in future releases.

% Copyright 2020 MathWorks, Inc.

% Sort curvatures
[~, sortedIds] = sort(curvatures);
minClusterSize = 100;
numPts  = size(normals,1);
numberOfNeighbors = size(neighInds,1);
notVisited = true(numPts, 1);
seedNotAdded = true(numPts, 1);
labels = ones(numPts, 1);
labels(isnan(curvatures)) = 0;
numValidPts = sum(~isnan(curvatures));
for i = 1 :numValidPts
    
    if notVisited(sortedIds(i))
        % Mark current idx as visited and take it as seed point
        notVisited(sortedIds(i)) = false;
        seedPtIdx =  sortedIds(i);
        regionPts = [];
        while ~isempty(seedPtIdx)
            curIdx = seedPtIdx(1);
            % Check whether the neighbors of seed point satisfy criteria
            for j = 1:numberOfNeighbors
                idx = neighInds(j, curIdx);
                if notVisited(idx)
                    normDiff = abs(dot(normals(curIdx, :), normals(idx, :)));
                    if normDiff >= norThres
                        notVisited(idx) = false;
                        regionPts(end+1) = idx;
                        
                    end
                    if curvatures(idx) <= curThres && seedNotAdded(idx)
                        seedNotAdded(idx) = false;
                        seedPtIdx(end+1) = idx;
                    end
                end
            end
            seedPtIdx(1) = [];
        end
        % Consider the regions that have atleast minClusterSize points
        if length(regionPts) > minClusterSize
            labels(regionPts) = 2;
        end
    end
end
end