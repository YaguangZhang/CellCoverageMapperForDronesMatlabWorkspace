%RECREATEOPENSIGNALMAP A snippet to recreate the map canvas captured in the
%reference Opensignal measurement result screenshot for IN.2021.
%
% Yaguang Zhang, Purdue, 03/23/2022

% Ref: the size of the Opensignal figure is 1971x3006.
targetFigHInPixel = 600;
targetFigSize = [targetFigHInPixel/3006*1971, targetFigHInPixel];
greyOutAlpha = 0.33;

targetLatRange = [37.644445, 41.899804];
targetLonRange = [-88.235131 -84.631344];
figure('Position', [0,0,targetFigSize]);
set(gca, 'InnerPosition', [0,0,1,1]);
hold on;
plot(lonsBoundShrinkedIN, latsBoundShrinkedIN, 'k-');
xlim(targetLonRange); ylim(targetLatRange);
hGM = plot_google_map;
xticklabels(''); yticklabels('');
axis manual;

% Greyout regions out of the AoI.
polyshapeSrhinkedIN = polyshape(lonsBoundShrinkedIN, latsBoundShrinkedIN);
polyshapeGreyOut = polyshape( ...
    [targetLonRange(1), targetLonRange(1), ...
    targetLonRange(2), targetLonRange(2)], ...
    [targetLatRange(1), targetLatRange(2), ...
    targetLatRange(2),targetLatRange(1)]);
polyshapeGreyOut = subtract(polyshapeGreyOut, polyshapeSrhinkedIN);
hGrey = plot(polyshapeGreyOut, ...
    'FaceColor', 'k', 'LineStyle', 'none', 'FaceAlpha', greyOutAlpha);

% EOF