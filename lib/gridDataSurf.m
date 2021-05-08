function [hSurf] = gridDataSurf(XYZs, numOfPtsPerSide, varargin)
%GRIDDATASURF Surf with arbitrary (x, y) points by griddata.
%
% Ref:
%   https://www.mathworks.com/help/matlab/math/interpolating-scattered-data.html
%
% Inputs:
%   - XYZs
%     The matrix containing the reference (x,y,z) data.
%   - numOfPtsPerSide
%     We will create a new grid for the visulization.
%
% Output:
%   - hSurf
%
% Yaguang Zhang, Purdue, 03/03/2021

[xi,yi] = meshgrid( ...
    linspace(min(XYZs(:,1)), max(XYZs(:,1)), numOfPtsPerSide), ...
    linspace(min(XYZs(:,2)), max(XYZs(:,2)), numOfPtsPerSide));
zi = griddata(XYZs(:,1), XYZs(:,2), XYZs(:,3), xi, yi);
hSurf = surf(xi, yi, zi, varargin{:});

end
% EOF