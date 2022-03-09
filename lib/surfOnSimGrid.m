function [hSurf] = surfOnSimGrid(gridLonLatZs, simConfigs, ...
    contourLevelVs, faceAlpha, fontSize, lineWidth, clabelTextColor)
%SURFONSIMGRID A customized surf function for the simulation grid.
%
% Yaguang Zhang, Purdue, 01/07/2021

if ~exist('faceAlpha', 'var')
    faceAlpha = 0.5;
end

if ~exist('fontSize', 'var')
    fontSize = 10;
end

if ~exist('lineWidth', 'var')
    fontSize = 0.5;
end

if ~exist('clabelTextColor', 'var')
    clabelTextColor = 'k';
end

% Create meshgrid for surf.
upSampFactor = 10;
sufNumPtsPerSide = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE ...
    .*upSampFactor;
lons = gridLonLatZs(:, 1);
lats = gridLonLatZs(:, 2);
zs = gridLonLatZs(:, 3);

% Find the ranges for the boundary of interet (BoI) to build a new grid for
% showing the results.
[latsBoI, lonsBoI] = simConfigs.utm2deg_speZone( ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,1), ...
    simConfigs.UTM_X_Y_BOUNDARY_OF_INTEREST(:,2));
lonMinBoI = min(lonsBoI);
lonMaxBoI = max(lonsBoI);
latMinBoI = min(latsBoI);
latMaxBoI = max(latsBoI);

[lonsNew, latsNew] = meshgrid( ...
    linspace(lonMinBoI, lonMaxBoI, sufNumPtsPerSide), ...
    linspace(latMinBoI, latMaxBoI, sufNumPtsPerSide));
zsNew = griddata(lons,lats,zs,lonsNew,latsNew);

% Replace the resulting NaN grid data points with nearest neighbor's value
% in the original dataset.
ziNearest = griddata(lons,lats,zs,lonsNew,latsNew,'Nearest');
boolsZsNewIsNan = isnan(zsNew);
zsNew(boolsZsNewIsNan) = ziNearest(boolsZsNewIsNan);

% For positions with NaN values in ziNearest, zsNew should also have NaN.
zsNew(isnan(ziNearest)) = nan;

% Ignore points out of the area of interest by seting the z values for them
% to NaN.
boolsInOrOnPoly = InPolygon(lonsNew(:), latsNew(:), lonsBoI, latsBoI);
boolsPtsToIgnore = ~boolsInOrOnPoly;
if any(boolsPtsToIgnore)
    zsNew(boolsPtsToIgnore) = nan;
end

hSurf = surf(lonsNew,latsNew,zsNew, ...
    'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
caxis([min(zs) max(zs)]);
if exist('contourLevelVs', 'var')
    [c,h] = contour(lonsNew,latsNew,zsNew,contourLevelVs);
else
    [c,h] = contour(lonsNew,latsNew,zsNew);
end
if exist('contourLevelVs', 'var')
    h.LineWidth = lineWidth;
end
h.ContourZLevel = max(zs)+1;
if exist('contourLevelVs', 'var')
    clabel(c, h, contourLevelVs, ...
        'FontSize', fontSize, 'Color', clabelTextColor, ...
        'LineWidth', lineWidth);
else
    clabel(c, h, ...
        'FontSize', fontSize, 'Color', clabelTextColor, ...
        'LineWidth', lineWidth);
end
end
% EOF