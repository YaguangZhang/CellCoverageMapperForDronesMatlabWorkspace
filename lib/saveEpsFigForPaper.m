function [ ] = saveEpsFigForPaper(hFig, fullPathToSaveFig, flagFixTicks)
%SAVEEPSFIGFORPAPER Save an .eps version for figure hFig at
%fullPathToSaveFig.
%
% Note: In order to use export_fig, one needs to install ghostscript at
%       http://www.ghostscript.com
% and manually locate pdftops (if asked) at
%       NistMeasurementCampaignCode/lib/ext/xpdf-tools-win-4.00
%
% Inputs:
%   - hFig
%     The handle for the figure to export.
%   - fullPathToSaveFig
%     The full path (the folder + file name) to save the figure.
%   - flagFixTicks
%     Optional. Default to true (sometimes this makes the tick labels look
%     better).
%
% Outputs:
%   - An eps file and a png file for the figure at the desired directory.
%
% Yaguang Zhang, Purdue, 07/16/2018

if ~exist('flagFixTicks', 'var')
    flagFixTicks = true;
end

[dirToSave, figName, fileExt] = fileparts(fullPathToSaveFig);
% Create directories if necessary.
if exist(dirToSave, 'dir')~=7
    mkdir(dirToSave);
end
if (~isempty(fileExt)) && (~strcmpi(fileExt, '.eps'))
    warning( ...
        'The file extention specified is not .eps and will be ignored.')
end
epsFullPathToSave = fullfile(dirToSave, [figName, '.eps']);
% Also save a .png copy for easy preview.
pngFullPathToSave = fullfile(dirToSave, [figName, '.png']);

% Set background to (non-transparent) white.
curFigure = gcf;
set(0, 'currentfigure', hFig);
curFigureColor = get(gcf,'Color');
set(gcf, 'Color', 'white');

if flagFixTicks
    % Fix ticks.
    curXTicks = xticks;
    curYTicks = yticks;
    xticks('manual'); xticks(curXTicks);
    yticks('manual'); yticks(curYTicks);
end

try
    export_fig(epsFullPathToSave, '-eps');
catch
    saveas(hFig, epsFullPathToSave, 'epsc');
end

try
    export_fig(pngFullPathToSave, '-png', '-transparent');
catch
    saveas(hFig, pngFullPathToSave);
end

set(gcf, 'Color', curFigureColor);
set(0, 'currentfigure', curFigure);
end

% EOF