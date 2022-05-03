function [ locIndicesForAllWorkers ] ...
    = preassignTaskIndicesToWorkers(numOfTasks, numOfAvailableWorkers, ...
    CHUNK_STYLE, gridInfo)
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
%   - CHUNK_STYLE
%     Optional. Default to "normal". Supported styles:
%       - 'normal'
%         1:numOfTasks will be divided into continous chunks (almost)
%         evenly.
%       - 'randomized'
%         1:numOfTasks will be permutated randomly first before being
%         divided into continous chunks (almost) evenly.
%       - 'rectangle'
%         The locations will be grouped as rectangle chunks, each with
%         almost the same number of elements. This option will require the
%         input gridInfo for more details about the input locations:
%           - gridInfo.xys
%             A matrix of the locations of interest, each row being (x, y)
%             coordinates of one location, cooresponding to one task. Note
%             that the grid could have an arbitrary shape/outline.
%           - gridInfo.resolution
%             The spatial resoltuion of the grid.
%
% Output:
%   - locIndicesForAllWorkers
%     A cell containing the task indices (as a row vector) that need to be
%     processed by each worker.
%
% Debugging figure:
%     figure; axis equal; hold on;
%      plot(gridInfo.xys(:, 1), gridInfo.xys(:, 2), 'k.');
%     for idxW = 1:numOfAvailableWorkers
%         plot(gridInfo.xys(locIndicesForAllWorkers{idxW}, 1), ...
%             gridInfo.xys(locIndicesForAllWorkers{idxW}, 2), 'x');
%     end
%
% Yaguang Zhang, Purdue, 05/02/2022

if numOfTasks == 0
    locIndicesForAllWorkers = {};
    return;
end

if ~exist('CHUNK_STYLE', 'var')
    CHUNK_STYLE = 'normal';
end

if (~exist('numOfAvailableWorkers', 'var')) || isnan(numOfAvailableWorkers)
    localCluster = gcp;
    numOfAvailableWorkers = localCluster.NumWorkers;
end

locIndicesForAllWorkers = cell(numOfAvailableWorkers, 1);

switch lower(CHUNK_STYLE)
    case {'normal', 'randomized'}
        if strcmpi(CHUNK_STYLE, 'normal')
            FLAG_RANDOMIZE_TASK = false;
        else
            FLAG_RANDOMIZE_TASK = true;
        end

        if FLAG_RANDOMIZE_TASK
            rng(1);
            randIntegerTaskIds = randperm(numOfTasks);
        end

        % Pre-allocate the indices to workers. We will use a worker-based
        % approach to avoid frequently changing cell size.
        for idxWorker = 1:numOfAvailableWorkers
            if FLAG_RANDOMIZE_TASK
                locIndicesForAllWorkers{idxWorker} ...
                    = randIntegerTaskIds( ...
                    idxWorker:numOfAvailableWorkers:numOfTasks);
            else
                % Essentially, we need numbers within the range
                % 1:numOfTasks that have a remainder of idxWorker after
                % being divided by numOfAvailableWorkers.
                locIndicesForAllWorkers{idxWorker} ...
                    = idxWorker:numOfAvailableWorkers:numOfTasks;
            end
        end
    case 'rectangle'
        assert(exist('gridInfo', 'var'), ...
            ['Input gridInfo is required for task assignment pattern ', ...
            '"rectangle"!']);
        assert(numOfTasks==size(gridInfo.xys,1), ...
            ['Input numOfTasks does not agree with ', ...
            'the row number of gridInfo.xys!']);

        numOfTasksPerWorker = ceil(numOfTasks/numOfAvailableWorkers);
        integerTaskIds = 1:numOfTasks;

        % Find the center rectangle chunk. To speed things up, we will
        % start with a rectangle which will contain enough locs if it is
        % pacted fully with grid points based on the given resolution.
        centerXY = mean(gridInfo.xys, 1);
        rectSide = sqrt(numOfTasksPerWorker)*gridInfo.resolution;
        centerChunkBoundLeftX = centerXY(1) - rectSide/2;
        centerChunkBoundRightX = centerXY(1) + rectSide/2;
        centerChunkBoundTopY = centerXY(2) + rectSide/2;
        centerChunkBoundBottomY = centerXY(2) - rectSide/2;

        boolsLocsInChunk = (gridInfo.xys(:,1)>=centerChunkBoundLeftX ...
            & gridInfo.xys(:,1)<=centerChunkBoundRightX ...
            & gridInfo.xys(:,2)<=centerChunkBoundTopY ...
            & gridInfo.xys(:,2)>=centerChunkBoundBottomY)';
        extensionCnt = 0;
        while sum(boolsLocsInChunk)<numOfTasksPerWorker
            switch mod(extensionCnt, 4)
                case 0
                    centerChunkBoundLeftX = centerChunkBoundLeftX ...
                        - gridInfo.resolution;
                case 1
                    centerChunkBoundRightX = centerChunkBoundRightX ...
                        + gridInfo.resolution;
                case 2
                    centerChunkBoundTopY = centerChunkBoundTopY ...
                        + gridInfo.resolution;
                case 3
                    centerChunkBoundBottomY = centerChunkBoundBottomY ...
                        - gridInfo.resolution;
            end

            extensionCnt = extensionCnt+1;

            boolsLocsInChunk = ...
                (gridInfo.xys(:,1)>=centerChunkBoundLeftX ...
                & gridInfo.xys(:,1)<=centerChunkBoundRightX ...
                & gridInfo.xys(:,2)<=centerChunkBoundTopY ...
                & gridInfo.xys(:,2)>=centerChunkBoundBottomY)';

            if all(boolsLocsInChunk)
                break;
            end
        end

        boolsAllTasksAssigned = boolsLocsInChunk;
        locIndicesForAllWorkers{1} = integerTaskIds(boolsLocsInChunk);

        if all(boolsAllTasksAssigned)
            return;
        else
            if numOfAvailableWorkers<13 %9
                % warning(['At least 9 workers are required for ', ...
                %     'task assignment pattern "rectangle"!', ...
                %      ' Currently we have ', ...
                %     num2str(numOfAvailableWorkers), ' free workers.']);

                % We will just divide the region along the long side.
                minX = min(gridInfo.xys(:,1));
                maxX = max(gridInfo.xys(:,1));
                rangeX = maxX - minX;

                minY = min(gridInfo.xys(:,2));
                maxY = max(gridInfo.xys(:,2));
                rangeY = maxY - minY;

                if rangeX>rangeY
                    startX = minX;
                    deltaX = (rangeX + gridInfo.resolution) ...
                        /numOfAvailableWorkers;

                    for idxW = 1:numOfAvailableWorkers
                        endX = startX + deltaX;

                        locIndicesForAllWorkers{idxW} = integerTaskIds( ...
                            gridInfo.xys(:,1)>=startX ...
                            & gridInfo.xys(:,1)<endX);

                        startX = startX + deltaX;
                    end
                else
                    startY = minY;
                    deltaY = (rangeY + gridInfo.resolution) ...
                        /numOfAvailableWorkers;

                    for idxW = 1:numOfAvailableWorkers
                        endY = startY + deltaY;

                        locIndicesForAllWorkers{idxW} = integerTaskIds( ...
                            gridInfo.xys(:,2)>=startY ...
                            & gridInfo.xys(:,2)<endY);

                        startY = startY + deltaY;
                    end
                end
                return;
            end
        end

        % Extend left, right, up, and down, to find more chunks.
        curWorkerCnt = 2;
        remainingXys = gridInfo.xys(~boolsAllTasksAssigned, :);
        remainingLocsIndices = integerTaskIds(~boolsAllTasksAssigned);

        % Left.
        leftChunkBoundLeftX = centerChunkBoundLeftX ...
            - gridInfo.resolution;
        leftChunkBoundRightX = centerChunkBoundLeftX;
        leftChunkBoundTopY = centerChunkBoundTopY;
        leftChunkBoundBottomY = centerChunkBoundBottomY;

        boolsLocsInChunk = false(1, numOfTasks);
        boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
            remainingXys(:,1)>=leftChunkBoundLeftX ...
            & remainingXys(:,1)<=leftChunkBoundRightX ...
            & remainingXys(:,2)<=leftChunkBoundTopY ...
            & remainingXys(:,2)>=leftChunkBoundBottomY))) = true;
        while any(remainingXys(:,1)<=leftChunkBoundRightX)
            while (sum(boolsLocsInChunk)<numOfTasksPerWorker) ...
                    && any(remainingXys(:,1)<leftChunkBoundLeftX)
                leftChunkBoundLeftX = leftChunkBoundLeftX ...
                    - gridInfo.resolution;

                boolsLocsInChunk = false(1, numOfTasks);
                boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                    remainingXys(:,1)>=leftChunkBoundLeftX ...
                    & remainingXys(:,1)<=leftChunkBoundRightX ...
                    & remainingXys(:,2)<=leftChunkBoundTopY ...
                    & remainingXys(:,2)>=leftChunkBoundBottomY))) = true;
            end

            locIndicesForAllWorkers{curWorkerCnt} ...
                = integerTaskIds(boolsLocsInChunk);
            curWorkerCnt = curWorkerCnt+1;

            assert(~any(boolsAllTasksAssigned & boolsLocsInChunk), ...
                'Loc index/indices got assigned more than once!');
            boolsAllTasksAssigned = ...
                boolsAllTasksAssigned | boolsLocsInChunk;
            remainingXys = gridInfo.xys(~boolsAllTasksAssigned, :);
            remainingLocsIndices = integerTaskIds(~boolsAllTasksAssigned);

            leftChunkBoundRightX = leftChunkBoundLeftX;
            leftChunkBoundLeftX = leftChunkBoundRightX ...
                - gridInfo.resolution;

            boolsLocsInChunk = false(1, numOfTasks);
            boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                remainingXys(:,1)>=leftChunkBoundLeftX ...
                & remainingXys(:,1)<=leftChunkBoundRightX ...
                & remainingXys(:,2)<=leftChunkBoundTopY ...
                & remainingXys(:,2)>=leftChunkBoundBottomY))) = true;
        end

        % Right.
        rightChunkBoundLeftX = centerChunkBoundRightX;
        rightChunkBoundRightX = centerChunkBoundRightX ...
            + gridInfo.resolution;
        rightChunkBoundTopY = centerChunkBoundTopY;
        rightChunkBoundBottomY = centerChunkBoundBottomY;

        boolsLocsInChunk = false(1, numOfTasks);
        boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
            remainingXys(:,1)>=rightChunkBoundLeftX ...
            & remainingXys(:,1)<=rightChunkBoundRightX ...
            & remainingXys(:,2)<=rightChunkBoundTopY ...
            & remainingXys(:,2)>=rightChunkBoundBottomY))) = true;
        while any(remainingXys(:,1)>=rightChunkBoundLeftX)
            while (sum(boolsLocsInChunk)<numOfTasksPerWorker) ...
                    && any(remainingXys(:,1)>rightChunkBoundRightX)
                rightChunkBoundRightX = rightChunkBoundRightX ...
                    + gridInfo.resolution;

                boolsLocsInChunk = false(1, numOfTasks);
                boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                    remainingXys(:,1)>=rightChunkBoundLeftX ...
                    & remainingXys(:,1)<=rightChunkBoundRightX ...
                    & remainingXys(:,2)<=rightChunkBoundTopY ...
                    & remainingXys(:,2)>=rightChunkBoundBottomY))) = true;
            end

            locIndicesForAllWorkers{curWorkerCnt} ...
                = integerTaskIds(boolsLocsInChunk);
            curWorkerCnt = curWorkerCnt+1;

            assert(~any(boolsAllTasksAssigned & boolsLocsInChunk), ...
                'Loc index/indices got assigned more than once!');
            boolsAllTasksAssigned = ...
                boolsAllTasksAssigned | boolsLocsInChunk;
            remainingXys = gridInfo.xys(~boolsAllTasksAssigned, :);
            remainingLocsIndices = integerTaskIds(~boolsAllTasksAssigned);

            rightChunkBoundLeftX = rightChunkBoundRightX;
            rightChunkBoundRightX = rightChunkBoundLeftX ...
                + gridInfo.resolution;

            boolsLocsInChunk = false(1, numOfTasks);
            boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                remainingXys(:,1)>=rightChunkBoundLeftX ...
                & remainingXys(:,1)<=rightChunkBoundRightX ...
                & remainingXys(:,2)<=rightChunkBoundTopY ...
                & remainingXys(:,2)>=rightChunkBoundBottomY))) = true;
        end

        % Top.
        topChunkBoundLeftX = centerChunkBoundLeftX;
        topChunkBoundRightX = centerChunkBoundRightX;
        topChunkBoundTopY = centerChunkBoundTopY + gridInfo.resolution;
        topChunkBoundBottomY = centerChunkBoundBottomY;

        boolsLocsInChunk = false(1, numOfTasks);
        boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
            remainingXys(:,1)>=topChunkBoundLeftX ...
            & remainingXys(:,1)<=topChunkBoundRightX ...
            & remainingXys(:,2)<=topChunkBoundTopY ...
            & remainingXys(:,2)>=topChunkBoundBottomY))) = true;
        while any(remainingXys(:,2)>=topChunkBoundBottomY)
            while (sum(boolsLocsInChunk)<numOfTasksPerWorker) ...
                    && any(remainingXys(:,2)>topChunkBoundTopY)
                topChunkBoundTopY = topChunkBoundTopY ...
                    + gridInfo.resolution;

                boolsLocsInChunk = false(1, numOfTasks);
                boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                    remainingXys(:,1)>=topChunkBoundLeftX ...
                    & remainingXys(:,1)<=topChunkBoundRightX ...
                    & remainingXys(:,2)<=topChunkBoundTopY ...
                    & remainingXys(:,2)>=topChunkBoundBottomY))) = true;
            end

            locIndicesForAllWorkers{curWorkerCnt} ...
                = integerTaskIds(boolsLocsInChunk);
            curWorkerCnt = curWorkerCnt+1;

            assert(~any(boolsAllTasksAssigned & boolsLocsInChunk), ...
                'Loc index/indices got assigned more than once!');
            boolsAllTasksAssigned = ...
                boolsAllTasksAssigned | boolsLocsInChunk;
            remainingXys = gridInfo.xys(~boolsAllTasksAssigned, :);
            remainingLocsIndices = integerTaskIds(~boolsAllTasksAssigned);

            topChunkBoundBottomY = topChunkBoundTopY;
            topChunkBoundTopY = topChunkBoundBottomY ...
                + gridInfo.resolution;

            boolsLocsInChunk = false(1, numOfTasks);
            boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                remainingXys(:,1)>=topChunkBoundLeftX ...
                & remainingXys(:,1)<=topChunkBoundRightX ...
                & remainingXys(:,2)<=topChunkBoundTopY ...
                & remainingXys(:,2)>=topChunkBoundBottomY))) = true;
        end

        % Bottom.
        bottomChunkBoundLeftX = centerChunkBoundLeftX;
        bottomChunkBoundRightX = centerChunkBoundRightX;
        bottomChunkBoundTopY = centerChunkBoundBottomY;
        bottomChunkBoundBottomY = centerChunkBoundBottomY ...
            - gridInfo.resolution;

        boolsLocsInChunk = false(1, numOfTasks);
        boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
            remainingXys(:,1)>=bottomChunkBoundLeftX ...
            & remainingXys(:,1)<=bottomChunkBoundRightX ...
            & remainingXys(:,2)<=bottomChunkBoundTopY ...
            & remainingXys(:,2)>=bottomChunkBoundBottomY))) = true;
        while any(remainingXys(:,2)<=bottomChunkBoundTopY)
            while (sum(boolsLocsInChunk)<numOfTasksPerWorker) ...
                    && any(remainingXys(:,2)<bottomChunkBoundBottomY)
                bottomChunkBoundBottomY = bottomChunkBoundBottomY ...
                    - gridInfo.resolution;

                boolsLocsInChunk = false(1, numOfTasks);
                boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                    remainingXys(:,1)>=bottomChunkBoundLeftX ...
                    & remainingXys(:,1)<=bottomChunkBoundRightX ...
                    & remainingXys(:,2)<=bottomChunkBoundTopY ...
                    & remainingXys(:,2)>=bottomChunkBoundBottomY))) = true;
            end

            locIndicesForAllWorkers{curWorkerCnt} ...
                = integerTaskIds(boolsLocsInChunk);
            curWorkerCnt = curWorkerCnt+1;

            assert(~any(boolsAllTasksAssigned & boolsLocsInChunk), ...
                'Loc index/indices got assigned more than once!');
            boolsAllTasksAssigned = ...
                boolsAllTasksAssigned | boolsLocsInChunk;
            remainingXys = gridInfo.xys(~boolsAllTasksAssigned, :);
            remainingLocsIndices = integerTaskIds(~boolsAllTasksAssigned);

            bottomChunkBoundTopY = bottomChunkBoundBottomY;
            bottomChunkBoundBottomY = bottomChunkBoundTopY ...
                - gridInfo.resolution;

            boolsLocsInChunk = false(1, numOfTasks);
            boolsLocsInChunk(integerTaskIds(remainingLocsIndices( ...
                remainingXys(:,1)>=bottomChunkBoundLeftX ...
                & remainingXys(:,1)<=bottomChunkBoundRightX ...
                & remainingXys(:,2)<=bottomChunkBoundTopY ...
                & remainingXys(:,2)>=bottomChunkBoundBottomY))) = true;
        end

        if all(boolsAllTasksAssigned)
            return;
        end

        % Recursived find the remaining chunks.
        boolsRemainingTopLeft = false(1, numOfTasks);
        boolsRemainingTopLeft(remainingLocsIndices( ...
            (remainingXys(:,1) < centerChunkBoundLeftX) ...
            & (remainingXys(:,2) > centerChunkBoundTopY))) = true;

        boolsRemainingBottomLeft = false(1, numOfTasks);
        boolsRemainingBottomLeft(remainingLocsIndices( ...
            (remainingXys(:,1) < centerChunkBoundLeftX) ...
            & (remainingXys(:,2) < centerChunkBoundBottomY))) = true;

        boolsRemainingBottomRight = false(1, numOfTasks);
        boolsRemainingBottomRight(remainingLocsIndices( ...
            (remainingXys(:,1) > centerChunkBoundRightX) ...
            & (remainingXys(:,2) < centerChunkBoundBottomY))) = true;

        boolsRemainingTopRight = false(1, numOfTasks);
        boolsRemainingTopRight(remainingLocsIndices( ...
            (remainingXys(:,1) > centerChunkBoundRightX) ...
            & (remainingXys(:,2) > centerChunkBoundTopY))) = true;

        % Assign workers.
        numOfCurFreeWorkers = numOfAvailableWorkers - curWorkerCnt + 1;
        boolsRemainingLocsForEachRegion = {boolsRemainingTopLeft, ...
            boolsRemainingBottomLeft, boolsRemainingBottomRight, ...
            boolsRemainingTopRight};
        numOfTasksForEachRegion = [sum(boolsRemainingTopLeft), ...
            sum(boolsRemainingBottomLeft), ...
            sum(boolsRemainingBottomRight), ...
            sum(boolsRemainingTopRight)];
        numOfWorkersForEachRegion = ceil(numOfCurFreeWorkers ...
            /sum(numOfTasksForEachRegion).*numOfTasksForEachRegion);
        while sum(numOfWorkersForEachRegion) > numOfCurFreeWorkers
            [maxNumOfAssignedWorkers, idxMaxNumOfAssignedWorkers] ...
                = max(numOfWorkersForEachRegion);
            numOfWorkersForEachRegion(idxMaxNumOfAssignedWorkers) ...
                = maxNumOfAssignedWorkers - 1;
        end

        % Top left.
        for idxRegion = 1:length(numOfWorkersForEachRegion)
            curBoolsRemainingLocs ...
                = boolsRemainingLocsForEachRegion{idxRegion};
            curGridInfo.xys = gridInfo.xys(curBoolsRemainingLocs, :);
            curGridInfo.resolution = gridInfo.resolution;
            curIndicesForWorkers = preassignTaskIndicesToWorkers( ...
                numOfTasksForEachRegion(idxRegion), ...
                numOfWorkersForEachRegion(idxRegion), ...
                CHUNK_STYLE, curGridInfo);

            % Convert the indices in the new grid to the bigger/old one.
            curLocIndicesRemainingLocs ...
                = integerTaskIds(curBoolsRemainingLocs);
            curIndicesForWorkers = cellfun( ...
                @(indices) curLocIndicesRemainingLocs(indices), ...
                curIndicesForWorkers, 'UniformOutput', false);

            curNumOfWorkers = length(curIndicesForWorkers);
            locIndicesForAllWorkers( ...
                curWorkerCnt:(curWorkerCnt+curNumOfWorkers-1)) ...
                = curIndicesForWorkers;
            curWorkerCnt = curWorkerCnt+curNumOfWorkers;

            boolsLocsInRegion = false(1, numOfTasks);
            boolsLocsInRegion([curIndicesForWorkers{:}]) = true;

            assert(sum(any( ...
                boolsAllTasksAssigned & boolsLocsInRegion)) ==0, ...
                'Loc index/indices got assigned more than once!');
            boolsAllTasksAssigned = ...
                boolsAllTasksAssigned | boolsLocsInRegion;
        end

        assert(sum(boolsAllTasksAssigned)==numOfTasks, ...
            'Not all locations are assigned!');
    otherwise
        error(['Unsupported task assignment pattern: ', CHUNK_STYLE, '!']);
end

end
% EOF