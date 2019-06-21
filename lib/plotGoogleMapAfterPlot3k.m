function [outputArg1,outputArg2] = plotGoogleMapAfterPlot3k(hFig, mapType)
%PLOTGOOGLEMAPAFTERPLOT3K Add a Google map background on a plot3k figure.
%
% The command plot_google_map messes up the color legend of plot3k. We
% will take care of this issue here.
%
% Yaguang Zhang, Purdue, 06/19/2019

plot_google_map('MapType', mapType);

hCb = findall( allchild(hFig), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap)+1,length(hCb.TickLabels));

end
% EOF