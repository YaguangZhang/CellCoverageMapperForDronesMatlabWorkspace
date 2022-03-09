function [hFig] = surfWithContour(XYZs, labels, figVisible)
%SURFWITHCONTOUR Similar to plot3k, but use surf and contour for show data
%with color.
%
% We will plot the data on a new grid with z values obtained by linear
% interpolation.
%
% Inputs:
%   - XYZs
%     The matrix for the input data with each row being (x, y, z) for one
%     sample.
%   - labels
%     A cell for
%           {title, xlabel, ylabel, zlabel, clabel, ctitle}.
%   - figVisible
%     A illustration plot for the resulting ground roughness will also be
%     generated. This flag constrols whether this plot will be visible.
%
% Output:
%   - hFig
%     The handle to the resultant figure.
%
% Yaguang Zhang, Purdue, 11/22/2019

%% Inputs.

xs = XYZs(:, 1);
ys = XYZs(:, 2);
zs = XYZs(:, 3);

colorRange = [min(zs), max(zs)];

titleToSet = labels{1};
xLabelToSet = labels{2};
yLabelToSet = labels{3};
zLabelToSet = labels{4};
cLabelToSet = labels{5};
cTitleToSet = labels{6};

upSampFactor = 10;

%% Create meshgrid for surf.

% We first construct a polygon for the area of interest.
indicesLidarAreaConvHull = convhull(XYZs(:,1), XYZs(:,2));
UTM_X_Y_BOUNDARY_OF_INTEREST ...
    = [XYZs(indicesLidarAreaConvHull,1), ...
    XYZs(indicesLidarAreaConvHull,2)];

xsBoI = UTM_X_Y_BOUNDARY_OF_INTEREST(:, 1);
ysBoI = UTM_X_Y_BOUNDARY_OF_INTEREST(:, 2);

xMinBoI = min(xsBoI);
xMaxBoI = max(xsBoI);
yMinBoI = min(ysBoI);
yMaxBoI = max(ysBoI);

distsToFirstSamp = vecnorm((XYZs(:,1:2)-XYZs(1,1:2))');
dataSpacialResoluation = min(distsToFirstSamp(distsToFirstSamp>0));
sufNumPtsPerSide = max([xMaxBoI-xMinBoI; yMaxBoI-yMinBoI]) ...
    ./dataSpacialResoluation...
    .*upSampFactor;

[xsNew, ysNew] = meshgrid( ...
    linspace(xMinBoI, xMaxBoI, sufNumPtsPerSide), ...
    linspace(yMinBoI, yMaxBoI, sufNumPtsPerSide));
zsNew = griddata(xs, ys, zs, xsNew, ysNew);

% Ignore points out of the area of interest by seting the z values for them
% to NaN.
boolsInOrOnPoly = InPolygon(xsNew(:), ysNew(:), xsBoI, ysBoI);
boolsPtsToIgnore = ~boolsInOrOnPoly;
if any(boolsPtsToIgnore)
    zsNew(boolsPtsToIgnore) = nan;
end

%% Plot.

hFig = figure('Visible', figVisible);
hold on;
surf(xsNew, ysNew, zsNew, ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
caxis(colorRange); 
title(titleToSet); 
xlabel(xLabelToSet); ylabel(yLabelToSet); zlabel(zLabelToSet);
hCb = colorbar; ylabel(hCb, cLabelToSet);
title(hCb, cTitleToSet);
[c,h] = contour(xsNew, ysNew, zsNew);
clabel(c,h); view(2); axis equal; axis tight;
xticks([]); yticks([]);

end
% EOF