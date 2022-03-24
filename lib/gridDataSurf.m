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
zi = griddata(XYZs(:,1), XYZs(:,2), XYZs(:,3), xi, yi, 'linear');
% Remove interpolation results out of the area covered by the input data.
% We will use a very big shrink factor to avoid (1) including to much
% unwanted area and (2) removing any wanted area.
indicesBoundPts = boundary(XYZs(:,1:2), 0.99);
isInvalidZi = ~InPolygon(xi(:), yi(:), ...
    XYZs(indicesBoundPts,1), XYZs(indicesBoundPts,2));
zi(isInvalidZi) = nan;
hSurf = surf(xi, yi, zi, varargin{:});

end
% EOF