function [hFig, covRatios] = plotCovRatioVsInspectedHeight( ...
    mapsCell, heightsInM, flagVisible)
%PLOTCOVRATIOVSINSPECTEDHEIGHT Plot the coverage ratio over inspected
%height.
%
% Inputs:
%   - mapsCell
%     A cell of path loss maps for the input heights. Each of mapsCell's
%     element is a vector of path losses in dB. Blocage is indicated by NaN
%     values.
%   - heightsInM
%     The inspected heights.
%   - flagVisible
%     An optional boolean for whether to show the resultant figure or not.
%     Default to true.
%
% Outputs:
%   - hFig
%     The handle to the resultant figure.
%   - covRatios
%     The coverage ratios for the input maps.
%
% Yaguang Zhang, Purdue, 10/17/2019

% Set an appropriate figure size for publication.
desiredFigSizeInPixel = [500, 180];

if ~exist('flagVisible', 'var')
    flagVisible = true;
end

covRatios = cellfun(@(m) sum(~isnan(m))./length(m), mapsCell);
if ~iscolumn(covRatios)
    covRatios = covRatios';
end

% Plot.
hFig= figure('Visible', flagVisible, ...
    'Position', [0, 0, desiredFigSizeInPixel]);
hold on; set(gca, 'fontWeight', 'bold');
plot(heightsInM, covRatios, '*--', 'MarkerSize', 6, 'LineWidth', 1);
xlabel('Height (m)'); ylabel('LoS Coverage Ratio'); grid on; grid minor;

end
% EOF