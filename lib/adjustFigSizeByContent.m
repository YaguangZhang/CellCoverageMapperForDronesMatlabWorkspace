function [ ] = adjustFigSizeByContent(hFig, minAxis, refSide, ...
    weightForWidth)
%ADJUSTFIGSIZEBYCONTENT Adjust the size of the figure automatically
%according to its content.
%
% Inputs:
%   - hFig
%     The handle to the target figure.
%   - minAxis
%     The tighest desired axis range to show.
%   - refSide
%     An optional string ('width' or 'height') to control which side of the
%     figure dose not change size.
%   - weightForWidth
%     An optional float (default 1) to indicate the weight of the width for
%     further adjust the figure size (weightForHeight is assumed to be 1).
%     For example, set this to 2 will set the width over height ratio twice
%     as that used in the default case.
%
% Yaguang Zhang, Purdue, 10/17/2019

% Weight the width a little more for more natrual figures.
WIDTH_WEIGHT_FACTOR = 1.2;

if (~exist('minAxis', 'var')) || (length(minAxis)~=4)
    figure(hFig);
    axis tight;
    minAxis = axis;
end

if ~exist('weightForWidth', 'var')
    weightForWidth = 1;
end

minAxisWidth = minAxis(2)-minAxis(1);
minAxisHeight = minAxis(4)-minAxis(3);

set(hFig, 'Unit', 'pixels');
figPos = get(hFig, 'Position');
figPosWidth = figPos(3);
figPosHeight = figPos(4);

refWidth = minAxisWidth;
if exist('refSide', 'var')
    switch lower(refSide)
        case 'width'
            refWidth = inf;
        case 'height'
            refWidth = -inf;
        otherwise
            error(['Unsupported refSide ', refSide, '!']);
    end
end

% We will adjust the figure according to the longer axis side.
if refWidth>minAxisHeight
    figPosHeight = ceil(figPosWidth./minAxisWidth.*minAxisHeight ...
        ./WIDTH_WEIGHT_FACTOR)./weightForWidth;
else
    figPosWidth = ceil(figPosHeight./minAxisHeight.*minAxisWidth ...
        .*WIDTH_WEIGHT_FACTOR).*weightForWidth;
end
set(hFig, 'Position', [figPos(1:2) figPosWidth figPosHeight]);

axis auto;
axis(minAxis);

end
% EOF