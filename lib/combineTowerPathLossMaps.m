function [coverageMap, coverageMapXLabels, coverageMapYLabels] ...
    = combineTowerPathLossMaps(towerPathLossMaps, ...
    towerPathLossMapsXLabels, towerPathLossMapsYLabels, ...
    allowedPathLossRangeInDB)
%COMBINETOWERPATHLOSSMAPS Combine tower path loss maps into one signal path
%loss map.
%
% Inputs:
%   - towerPathLossMaps
%     A cell with path loss matrices for towers.
%   - towerPathLossMapsXLabels, towerPathLossMapsYLabels
%     A cell with vectors for x labels and y labels, repectively, of the
%     tower path loss maps.
%   - allowedPathLossRangeInDB
%     The allowed path loss range in dB in terms of [minPathLossAllowed,
%     maxPathLossAllowed]. Any path loss out of this range will be
%     considered as inf.
%
% Outputs:
%   - coverageMap
%     The combined path loss matrix for the coverage map.
%   - coverageMapXLabels, coverageMapYLabels
%     Two vectors for x labels and y labels, repectively, of the resulting
%     coverage map.
%
% Yaguang Zhang, Purdue, 06/19/2019

coverageMapXLabels = towerPathLossMapsXLabels{1};
coverageMapYLabels = towerPathLossMapsYLabels{1};
coverageMap = towerPathLossMaps{1};

% Make sure the input tower path loss maps are covering the save area.
numTowerPathLossMaps = length(towerPathLossMaps);
for idxTowerPLM = 2:numTowerPathLossMaps
    assert( all(towerPathLossMapsXLabels{idxTowerPLM}(:) ...
        == coverageMapXLabels(:)) ...
        && all(towerPathLossMapsYLabels{idxTowerPLM}(:) ...
        == coverageMapYLabels(:)), ...
        'The input tower path loss maps are not covering the same area!');
    
    % Combine the input maps by keeping the minimum path loss value for
    % each location of interest.
    coverageMap = min(coverageMap, towerPathLossMaps{idxTowerPLM});
end

% Clip the result to the limits.
boolsSetToNan = (coverageMap(:)<allowedPathLossRangeInDB(1)) ...
    | (coverageMap(:)>allowedPathLossRangeInDB(2));
coverageMap(boolsSetToNan) = nan;

end
% EOF