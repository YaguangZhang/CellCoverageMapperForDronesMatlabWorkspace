% TESTPREASSIGNTASKINDICES A snippet for testing different approach of
% assigning locations.
%
% Yaguang Zhang, Purdue, 05/02/2022

cd(fileparts(mfilename('fullpath')));

addpath('.'); addpath(fullfile('..', 'lib'));
load('testPreassignTaskIndices.mat');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers);
figure; axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Normal');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers, 'randomized');
figure; axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Randomized');

locIndicesForAllWorkers = preassignTaskIndicesToWorkers( ...
    numOfTasks, numOfAvailableWorkers, ...
    'rectangle', gridInfo);
figure; axis equal; hold on;
plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
for idxW = 1:numOfAvailableWorkers
    plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
        gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
end
title('Rectangle');
% EOF