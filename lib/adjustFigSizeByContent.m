function [ ] = adjustFigSizeByContent(hFig, minAxis)
%ADJUSTFIGSIZEBYCONTENT Adjust the size of the figure automatically
%according to its content.
%
% Inputs:
%   - hFig
%     The handle to the target figure. 
%   - minAxis
%     The tighest desired axis range to show.
%
% Yaguang Zhang, Purdue, 10/17/2019

% Weight the width a little more for more natrual figures.
WIDTH_WEIGHT_FACTOR = 1.2;

if ~exist('minAxis', 'var')
    gcf(hFig);
    axis tight;
    minAxis = axis;
end

minAxisWidth = minAxis(2)-minAxis(1);
minAxisHeight = minAxis(4)-minAxis(3);

set(hFig, 'Unit', 'pixels');
figPos = get(hFig, 'Position');
figPosWidth = figPos(3);
figPosHeight = figPos(4);

% We will adjust the figure according to the longer axis side.
if minAxisWidth>minAxisHeight
    figPosHeight = ceil(figPosWidth./minAxisWidth.*minAxisHeight ...
        ./WIDTH_WEIGHT_FACTOR);
else
    figPosWidth = ceil(figPosHeight./minAxisHeight.*minAxisWidth ...
        .*WIDTH_WEIGHT_FACTOR);
end
set(hFig, 'Position', [figPos(1:2) figPosWidth figPosHeight]);

axis auto;
axis(minAxis);

end
% EOF