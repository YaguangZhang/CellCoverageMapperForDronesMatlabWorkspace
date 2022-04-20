function [in, on] = InPolygon(xs, ys, xsBoundary, ysBoundary)
%INPOLYGON Replace the MEX function InPolygon with inpoly2.
%
% All inputs needs to be column vectors. Please refer to the built-in
% Matlab function inpolygon, the 3rd-party MEX implementation InPolygon,
% and the 3rd-party Matlab implementation inpoly2 for more details.
%
% Yaguang Zhang, Purdue, 04/19/2022

[in, on] = inpoly2([xs, ys], [xsBoundary, ysBoundary]);

end
% EOF