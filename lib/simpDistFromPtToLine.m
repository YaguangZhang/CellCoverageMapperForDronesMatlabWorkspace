function [ds] = simpDistFromPtToLine(pts, v1, v2)
%SIMPDISTFROMPTTOLINE A simple and fast function to calculate the 2D
%distance(s) between the input point(s), pts, and the line defined by v1
%and v2.
%
% Inputs:
%   - pts, v1, v2
%     Input query points and two different points along the line of
%     interest. Each row being a (x,y) pair.
%
% A test:
%   simpDistFromPtToLine([0,0;0,1;2,2], [0,0], [1,0])
% should output [0; 1; 2].
%
% A test:
%   simpDistFromPtToLine([-1,-1;0,0;0,1;2,2], [0,0], [1,1])
% should output [0; 0; 0.7071; 0].
%
% Yaguang Zhang, Purdue, 05/02/2022

a = v1 - v2;
bs = pts - v2;

if a(1) == 0
    boolsForceZeroArea = bs(:,1)==0;
elseif a(2) == 0
    boolsForceZeroArea = bs(:,2)==0;
else
    boolsToCheck = bs(:,2)~=0;
    ratio = a(1)./a(2);
    boolsForceZeroArea = false(length(boolsToCheck), 1);
    boolsForceZeroArea(boolsToCheck) = ...
        bs(boolsToCheck,2)./bs(boolsToCheck,1) == ratio;
end

if any(boolsForceZeroArea)
    bs(boolsForceZeroArea, :) = repmat([0,0], sum(boolsForceZeroArea), 1);
end

ds = abs(a(1).*bs(:,2) - a(2).*bs(:,1)) ./ norm(a);

end
% EOF