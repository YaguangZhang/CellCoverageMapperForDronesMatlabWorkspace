function [ locIndicesForAllWorkers ] ...
    = preassignTaskIndicesToWorkers(numOfTasks, numOfAvailableWorkers, ...
    FLAG_RANDOMIZE_TASK)
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
%   - FLAG_RANDOMIZE_TASK
%     Optional. Default to true. Set FLAG_RANDOMIZE_TASK to false to stop
%     randomizing the output indices.
%
% Output:
%   - locIndicesForAllWorkers
%     A cell containing the task indices (as a row vector) that need to be
%     processed by each worker.
%
% Yaguang Zhang, Purdue, 09/18/2019

if ~exist('FLAG_RANDOMIZE_TASK', 'var')
    FLAG_RANDOMIZE_TASK = true;
end

if ~exist('numOfAvailableWorkers', 'var')
    localCluster = gcp;
    numOfAvailableWorkers = localCluster.NumWorkers;
end

locIndicesForAllWorkers = cell(numOfAvailableWorkers, 1);

if FLAG_RANDOMIZE_TASK
    rng(1);
    randIntegerTaskIds = randperm(numOfTasks);
end

% Pre-allocate the indices to workers. We will use a worker-based approach
% to avoid frequently changing cell size.
for idxWorker = 1:numOfAvailableWorkers
    if FLAG_RANDOMIZE_TASK
        locIndicesForAllWorkers{idxWorker} ...
            = randIntegerTaskIds( ...
            idxWorker:numOfAvailableWorkers:numOfTasks);
    else
        % Essentially, we need numbers within the range 1:numOfTasks that
        % have a remainder of idxWorker after being divided by
        % numOfAvailableWorkers.
        locIndicesForAllWorkers{idxWorker} ...
            = idxWorker:numOfAvailableWorkers:numOfTasks;
    end
end

end
% EOF