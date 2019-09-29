function [ gridVs ] = constructAxisGrid(meanValue, ...
    numOfGridPts, gridResolution)
%CONSTRUCTAXISGRID Construct a grid for a one-dimentional range.
%
% Inputs:
%   - meanValue
%     The center of the range to be covered.
%   - numOfGridPts
%     The number of points to use.
%   - gridResolution
%     A positive scalar indicating the spacial resolution for the grid.
% Output:
%   - gridVs
%     A column vector with the resulting grid values.
%
% Yaguang Zhang, Purdue, 09/11/2019

gridVs = ((1:numOfGridPts)-1).*gridResolution;
gridVs = gridVs - (gridVs(end)-gridVs(1))/2 + meanValue;

end
% EOF