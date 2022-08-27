function kmlStr = genKmlStrForGpsTrackAs3DLines( ...
    lons, lats, eles, ...
    datetimeStamps, maxAllowedTimeGapInS, ...
    name, lineColor, lineWidth, extrude)
%GENKMLSTRFORGPSTRACKAS3DLINES Generate .kml string for a GPS track.
%
%   The track will be broken into and plotted as separte lines if any time
%   gap between two adjacent samples is too big.
%
% Inputs:
%   - lons, lats, eles
%     Column vectores. The GPS coordinates (lats, lons) in degree with
%     elevations in meter for the track.
%   - datetimeStamps
%     GPS time stamps stored as a Matlab datetime (column) array.
%   - maxAllowedTimeGapInS
%     The maximum allowed GPS time gps in second within one line.
%   - name, lineColor, lineWidth, extrude
%     Plotting options. Please refer to ge_plot3 under
%     lib\ext\googleearth_matlab\googleearth for more details.
%
% Yaguang Zhang, Purdue, 08/24/2022

numOfPts = length(lons);
timeGapsInS = seconds(datetimeStamps(2:end) - datetimeStamps(1:(end-1)));
indicesToBreakTrack = 2:numOfPts;
indicesToBreakTrack = indicesToBreakTrack( ...
    abs(timeGapsInS)>maxAllowedTimeGapInS);

lonsToPlot = lons;
latsToPlot = lats;
elesToPlot = eles;

if ~isempty(indicesToBreakTrack)
    numOfBreaks = length(indicesToBreakTrack);
    tempIndicesToBreak = indicesToBreakTrack;

    for idxBreakPt = 1:numOfBreaks
        curIdxToBreak = tempIndicesToBreak(idxBreakPt);

        lonsToPlot = [lonsToPlot(1:(curIdxToBreak-1)); nan; ...
            lonsToPlot(curIdxToBreak:end)];
        latsToPlot = [latsToPlot(1:(curIdxToBreak-1)); nan; ...
            latsToPlot(curIdxToBreak:end)];
        elesToPlot = [elesToPlot(1:(curIdxToBreak-1)); nan; ...
            elesToPlot(curIdxToBreak:end)];

        tempIndicesToBreak = tempIndicesToBreak+1;
    end

    assert(length(lonsToPlot) == numOfPts + numOfBreaks, ...
        'Unexpected total number of points to plot!');
end

kmlStr = ge_plot3( ...
    lonsToPlot, latsToPlot, elesToPlot, ...
    'name', name, ...
    'lineColor', lineColor, ...
    'lineWidth', lineWidth, ...
    'extrude', extrude);
end

% EOF