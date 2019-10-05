function [ hProcTimeFig, hProcTimeHistFig ] = plotProcessingTime( ...
    simState, plotType, flagVisible)
%PLOTPROCESSINGTIME Plot the processing time for all pixels of the specific
%type of path loss map.
%
% Inputs:
%   - simState
%     The struct for the simulation results. We need field:
%       - TimeUsedInSForEachPixel
%         For debugging and evaluating computation performance of the
%         simulator, we also record the time used in second for each map
%         pixel. For example,
%               simState.TimeUsedInSForEachPixel{idxCell}{idxDroneHeight}
%         is the processing time vector (for locations in
%         simState.mapGridXYPts) for the idxCell-th cellular tower that has
%         effect on the simulation, at the idxDroneHeight-th drone height
%         that needs to be inspected.
%
%         We expect the time used for a pixel at its first drone height
%         will be way longer than for other heights, because for the other
%         heights, we will reuse the terrain profiles generated for the
%         first height.
%   - plotType
%     A case-insensitive string controling what type of plot to generate.
%     Currently support: 'pixels' for pixel-wise processing time;
%     'cellTowers' for cellular-tower-wise processing time.
%   - flagVisible
%     An optional boolean for whether to show the resultant figure or not.
%     Default to true.
%
% Output:
%   - hProcTimeFig
%     The handle to the resultant processing time figure.
%   - hProcTimeHistFig
%     The handle to the resultant processing time histogram figure.
%
% Yaguang Zhang, Purdue, 10/04/2019

% Number of bins to use in the hist plot.
NUM_OF_HIST_BINS = 50;

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

switch lower(plotType)
    case 'pixels'
        timeUsedForAllMapSets = [simState.TimeUsedInSForEachPixel{:}];
        % Time used for all pixels.
        ys = [timeUsedForAllMapSets{:}]';
    case 'mapsets'
        % Time used for all mapsets.
        ys = cell2mat(cellfun(@(timeVector) sum(timeVector), ...
            [simState.TimeUsedInSForEachPixel{:}], ...
            'UniformOutput', false)');
    case 'celltowers'
        % Time used for all cellular towers.
        ys = cell2mat(cellfun(@(timeVector) sum(cell2mat(timeVector)), ...
            simState.TimeUsedInSForEachPixel, ...
            'UniformOutput', false));
    otherwise
        error(['Unsupported plot type ', plotType, '!']);
end

xs = 1:length(ys);

titleToSet = ['Total time: ', seconds2human(sum(ys))];

maxY = max(ys);
secondsInAMinute = 60;
secondsInAnHour = secondsInAMinute*60;
secondsInADay = secondsInAnHour*24;
secondsInAWeek = secondsInADay*7;
if maxY > secondsInAWeek
    % Week.
    ys = ys/secondsInAWeek;
    ysUnit = 'Week';
elseif maxY > secondsInADay
    % Day.
    ys = ys/secondsInADay;
    ysUnit = 'Day';
elseif maxY > secondsInAnHour
    % Hour.
    ys = ys/secondsInAnHour;
    ysUnit = 'Hour';
elseif maxY > secondsInAMinute
    % Minute.
    ys = ys/secondsInAMinute;
    ysUnit = 'min';
else
    % Second.
    ysUnit = 's';
end

hProcTimeFig = figure('Visible', flagVisible); hold on;
hYs = plot(xs, ys, '.');
hAverage = plot([0; length(ys)], ones(2,1).*mean(ys));
legend([hYs, hAverage], 'Records', 'Average');
xlabel('Index'); ylabel(['Processing Time (', ysUnit, ')']); ...
    grid on; grid minor;
axis tight; title(titleToSet);

hProcTimeHistFig = figure('Visible', flagVisible);
hist(ys, NUM_OF_HIST_BINS);
xlabel(['Processing Time (', ysUnit, ')']); ylabel('Occurrence');
grid on; grid minor;
title(titleToSet);

end
% EOF