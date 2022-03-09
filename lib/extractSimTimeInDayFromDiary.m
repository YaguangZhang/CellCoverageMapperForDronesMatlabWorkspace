function [ simTimeInDTotal, simTimesInDForEffeTs, simTimeInDTotalMore ] ...
    = extractSimTimeInDayFromDiary(pathToSimDiary, numOfEffeTowers)
%EXTRACTSIMTIMEINDAYFROMDIARY Extract the total simulation time (in days)
%from a simulation diary file. We will find the processing time for each
%effective tower and sum them up for the total simulation time.
%
% Note: simTimeInDTotal only considers the simulation time (start to end
% for each effective tower); simTimeInDTotalMore considers some
% intermittent time between simulations, too (start of the first tower to
% end of the last tower).
%
% Yaguang Zhang, Purdue, 01/22/2022

fid = fopen(pathToSimDiary, 'r');

% We are looking for start and end lines for each effective tower
% processed. They look like:
%   [2021/05/11 22:25:43] Effective cellular tower #6/50 ...
% and
%   [2021/05/11 22:36:13] Finished cellular tower #6/50
timePattern = '\[(\d{4}\/\d{2}\/\d{2}\ \d{2}:\d{2}:\d{2})\]';
buildStartPatten = @(idx) [timePattern, ' Effective cellular tower #', ...
    num2str(idx), '\/', num2str(numOfEffeTowers)];
buildEndPatten = @(idx) [timePattern, ' Finished cellular tower #', ...
    num2str(idx), '\/', num2str(numOfEffeTowers)];

simTimesInDForEffeTs = nan(numOfEffeTowers, 1);
simTimeInDTotalMoreTimestrs = cell(2, 1);
try
    preTowerIdx = 0;
    % Continue reading the file until all sim times for effective towers
    % are extracted.
    while preTowerIdx<numOfEffeTowers
        curpreTowerIdx = preTowerIdx+1;
        curStartPattern = buildStartPatten(curpreTowerIdx);
        curEndPattern = buildEndPatten(curpreTowerIdx);

        curDateStrStart = '';
        curDateStrEnd = '';
        % Keep getting new lines until the end time str is found for this
        % tower.
        while isempty(curDateStrEnd)
            newline = fgetl(fid);

            matchedStrStaCell = regexp(newline, curStartPattern, 'tokens');
            if ~isempty(matchedStrStaCell)
                assert(length(matchedStrStaCell)==1 ...
                    && length(matchedStrStaCell{1})==1, ...
                    'More than one match found!');
                curDateStrStart = matchedStrStaCell{1}{1};
                if curpreTowerIdx == 1
                    simTimeInDTotalMoreTimestrs{1} = curDateStrStart;
                end
            end

            matchedStrEndCell = regexp(newline, curEndPattern, 'tokens');
            if ~isempty(matchedStrEndCell)
                assert(length(matchedStrEndCell)==1 ...
                    && length(matchedStrEndCell{1})==1, ...
                    'More than one match found!');
                curDateStrEnd = matchedStrEndCell{1}{1};
                if curpreTowerIdx == numOfEffeTowers
                    simTimeInDTotalMoreTimestrs{2} = curDateStrStart;
                end
            end
        end

        % The start time str for this tower is expected to be found, too,
        % by now.
        assert(~isempty(curDateStrStart), 'Start time str is missing!');

        simTimesInDForEffeTs(curpreTowerIdx) ...
            = days(datestrDiff(curDateStrStart, curDateStrEnd));

        preTowerIdx = preTowerIdx+1;
    end
catch err
    fclose(fid);
    rethrow(err);
end

fclose(fid);
simTimeInDTotal = sum(simTimesInDForEffeTs);
simTimeInDTotalMore = days(datestrDiff( ...
    simTimeInDTotalMoreTimestrs{1}, simTimeInDTotalMoreTimestrs{2}));
end
% EOF