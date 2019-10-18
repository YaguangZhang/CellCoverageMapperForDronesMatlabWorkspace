function [ ] = saveEpsFigForPaper(hFig, fullPathToSaveFig)
%SAVEEPSFIGFORPAPER Save an .eps version for figure hFig at
%fullPathToSaveFig.
%
% Note: In order to use export_fig, one needs to install ghostscript at
%       http://www.ghostscript.com 
% and manually locate pdftops (if asked) at 
%       NistMeasurementCampaignCode/lib/ext/xpdf-tools-win-4.00
%
% Yaguang Zhang, Purdue, 07/16/2018

[dirToSave, figName, fileExt] = fileparts(fullPathToSaveFig);
% Create directories if necessary.
if exist(dirToSave, 'dir')~=7
    mkdir(dirToSave);
end
if ~strcmpi(fileExt, '.eps')
    warning( ...
        'The file extention specified is not .eps and will be ignored.')
end

epsFullPathToSave = fullfile(dirToSave, [figName, '.eps']);
% Also save a .png copy for easy preview.
pngFullPathToSave = fullfile(dirToSave, [figName, '.png']);

curFigure = gcf;
set(0, 'currentfigure', hFig);
curFigureColor = get(gcf,'Color');
set(gcf, 'Color', 'white');

export_fig(epsFullPathToSave, '-eps'); 
export_fig(pngFullPathToSave, '-png', '-transparent');

set(gcf, 'Color', curFigureColor);
set(0, 'currentfigure', curFigure);
end

% EOF