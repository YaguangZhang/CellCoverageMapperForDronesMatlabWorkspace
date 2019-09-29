function [ locIndicesForAllWorkers ] ...
    = preassignTaskIndicesToWorkers(numOfTasks)
%PREASSIGNTASKINDICESTOWORKERS Preassign 1:numOfTasks to available workers.
%
% Input:
%   - numOfTasks
%     The number of different tasks to consider. Essentially, values
%     1:numOfTasks will be assigned to available workers (roughly) evenly.
%
% Output:
%   - locIndicesForAllWorkers
%     A cell containing the task indices (as a row vector) that need to be
%     processed by each worker.
%
% Yaguang Zhang, Purdue, 09/18/2019

localCluster = gcp;
numOfAvailableWorkers = localCluster.NumWorkers;

locIndicesForAllWorkers = cell(numOfAvailableWorkers, 1);

for idxLoc = 1:numOfTasks
    idxWorkerToAssign = mod(idxLoc-1, numOfAvailableWorkers)+1;
    locIndicesForAllWorkers{idxWorkerToAssign} = ...
        [locIndicesForAllWorkers{idxWorkerToAssign}, idxLoc];
end

end
% EOF