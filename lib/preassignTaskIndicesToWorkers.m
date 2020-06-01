function [ locIndicesForAllWorkers ] ...
    = preassignTaskIndicesToWorkers(numOfTasks, numOfAvailableWorkers)
%PREASSIGNTASKINDICESTOWORKERS Preassign 1:numOfTasks to available workers.
%
% Input:
%   - numOfTasks
%     The number of different tasks to consider. Essentially, values
%     1:numOfTasks will be assigned to available workers (roughly) evenly.
%   - numOfAvailableWorkers
%     Optional. The number of available workers. If not specified, a local
%     parallel computing pool will be started and the number of workers
%     available there will be used.
%
% Output:
%   - locIndicesForAllWorkers
%     A cell containing the task indices (as a row vector) that need to be
%     processed by each worker.
%
% Yaguang Zhang, Purdue, 09/18/2019

if ~exist('numOfAvailableWorkers', 'var')
    localCluster = gcp;
    numOfAvailableWorkers = localCluster.NumWorkers;
end

locIndicesForAllWorkers = cell(numOfAvailableWorkers, 1);

% Pre-allocate the indices to workers. We will use a worker-based approach
% to avoid frequently changing cell size.
for idxWorker = 1:numOfAvailableWorkers
    % Essentially, we need numbers within the range 1:numOfTasks that have
    % a remainder of idxWorker after being divided by
    % numOfAvailableWorkers.
    locIndicesForAllWorkers{idxWorker} ...
        = idxWorker:numOfAvailableWorkers:numOfTasks;
end

% % An index-based approach, which acheives the same goal but is slower.
% for idxLoc = 1:numOfTasks
%     idxWorkerToAssign = mod(idxLoc-1, numOfAvailableWorkers)+1;
%     locIndicesForAllWorkers{idxWorkerToAssign} = ...
%         [locIndicesForAllWorkers{idxWorkerToAssign}, idxLoc];
% end

end
% EOF