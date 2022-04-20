function [ hS, hCb ] = overlayOpenSigStyleMap(simConfigs, ...
    mapLatLonPts, mapVs, colorRange)
%OVERLAYOPENSIGSTYLEMAP Plot a simulation map on the current figure
%following the OpenSignal color scheme.
%
% Yaguang Zhang, Purdue, 03/23/2022

% A very small positive value. We will shift the input mapVs up by this
% amount to make sure everything is shown properly after add the Google
% Maps layer.
deltaZ = 0;

% Create meshgrid for surf.
upSampFactor = 10;
sufNumPtsPerSide = simConfigs.NUM_OF_PIXELS_FOR_LONGER_SIDE ...
    .*upSampFactor;
lons = mapLatLonPts(:, 2);
lats = mapLatLonPts(:, 1);

zs = mapVs+deltaZ;

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
[in,on] = inpoly2([lonsNew(:), latsNew(:)], [lonsBoI, latsBoI]);
boolsPtsToIgnore = ~(in|on);
if any(boolsPtsToIgnore)
    zsNew(boolsPtsToIgnore) = nan;
end

hS = surf(lonsNew, latsNew, zsNew, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5);
caxis(colorRange); zlim([0, max(mapVs)]);

% The OpenSignal app uses green (0, 175, 15) to red (1, 0, 0) => good to
% bad. Note that we have the higher the input value, the worse the
% situation.
numOfColorPts = 1000;
% Green to red.
openSigColorMap = [(linspace(0,1,numOfColorPts))', ...
    (linspace(175,0,numOfColorPts)./255)', ...
    (linspace(15,0,numOfColorPts)./255)'];
colormap(openSigColorMap)

% Add a color bar and modify the max label.
hCb =colorbar('Location', 'west', 'FontSize', 10, 'FontWeight', 'bold');
hCb.TickLabels{end} = ['≥', hCb.TickLabels{end}];
hCb.Position(4) = 0.9;

end
% EOF