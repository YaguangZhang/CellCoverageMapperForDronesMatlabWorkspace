function [newAxis] = extendAxisByFactor(curAxis, factorExtend)
%EXTENDAXISBYFACTOR Extend the axis area by the specified factor (for each
%side).
%
% Yaguang Zhang, Purdue, 10/17/2019

xMin = curAxis(1);
xMax = curAxis(2);
yMin = curAxis(3);
yMax = curAxis(4);

xMean = (xMin+xMax)/2;
xDelta = (xMax - xMin)*(1+factorExtend)/2;
xMin = xMean - xDelta;
xMax = xMean + xDelta;

yMean = (yMin+yMax)/2;
yDelta = (yMax - yMin)*(1+factorExtend)/2;
yMin = yMean - yDelta;
yMax = yMean + yDelta;

newAxis = [xMin xMax yMin yMax];

end
% EOF