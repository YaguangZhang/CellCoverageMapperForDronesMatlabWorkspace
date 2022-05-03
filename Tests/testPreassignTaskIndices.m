% TESTPREASSIGNTASKINDICES A snippet for testing different approach of
% assigning locations.
%
% Yaguang Zhang, Purdue, 05/02/2022

cd(fileparts(mfilename('fullpath')));

addpath('.'); addpath(fullfile('..', 'lib'));
load('testPreassignTaskIndices.mat');

curFigPos = [0,0,420,560];

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers);
figure('Position', curFigPos); axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Normal');
saveas(gcf, 'testPreassignTaskIndices_1_Normal.jpg');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers, 'randomized');
figure('Position', curFigPos); axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Randomized');
saveas(gcf, 'testPreassignTaskIndices_2_Randomized.jpg');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers, ...
    'rectangle', gridInfo);
figure('Position', curFigPos); axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Rectangle');
saveas(gcf, 'testPreassignTaskIndices_3_Rectangle.jpg');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers, ...
    'tile', gridInfo);
figure('Position', curFigPos); axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Tile');
saveas(gcf, 'testPreassignTaskIndices_4_Tile.jpg');

% EOF