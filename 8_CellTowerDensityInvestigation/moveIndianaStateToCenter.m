%MOVEINDIANASTATETOCENTER A helper script to center the Indiana State on a
%map.
%
% This also works for other area of interest, as long as
% boundOfInterestLats and boundOfInterestLons are properly set.
%
% Yaguang Zhang, Purdue, 01/07/2021

curAxis = axis;
halfSideLenX = (curAxis(2)- curAxis(1))/2;
halfSideLenY = (curAxis(4)- curAxis(3))/2;

bOILonCenter = mean([max(boundOfInterestLons), min(boundOfInterestLons)]);
bOILatCenter = mean([max(boundOfInterestLats), min(boundOfInterestLats)]);

minX = bOILonCenter-halfSideLenX;
maxX = bOILonCenter+halfSideLenX;
minY = bOILatCenter-halfSideLenY;
maxY = bOILatCenter+halfSideLenY;

axis([minX, maxX, minY, maxY]);
% EOF