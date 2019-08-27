% LOGEXECTIMEFORCOVMAPGEN Log the execution time for coverage map
% generation.
%
% Yaguang Zhang, Purdue, 08/14/2019

if ~exist('flagProcessInterrupted', 'var')
    flagProcessInterrupted = false;
end

numOfHsInspected = length(towerPathLossMapsEHata);
numOfCellAntsObserved = length(towerPathLossMapsEHata{1});

numOfMapsGenerated = numOfHsInspected.*numOfCellAntsObserved;

numOfPixelsCovered = 0;
numOfPixelsCoveredForHs = zeros(numOfHsInspected, 1);
for idxH = 1:numOfHsInspected
    assert(numOfCellAntsObserved==length(towerPathLossMapsEHata{idxH}), ...
        'The same number of cell towers should be observed for all maps!');
    for idxC = 1:numOfCellAntsObserved
        % Cound the number of pixels in each side of the map.
        [pixM, pixN] = size(towerPathLossMapsEHata{idxH}{idxC});
        curNumOfPixels = pixM.*pixN;
        numOfPixelsCoveredForHs(idxH) = numOfPixelsCoveredForHs(idxH) ...
            + curNumOfPixels;
        numOfPixelsCovered = numOfPixelsCovered + curNumOfPixels;
    end
end

pathToSaveExecTime = fullfile(pathToSaveResults, ...
    'executionTime.txt');
% Log the information to a txt file.
logFileId = fopen(pathToSaveExecTime,'wt');

fprintf(logFileId, '%s\n', fileNameHintRuler);
fprintf(logFileId, '%s\n', ['  ', curFileName]);
fprintf(logFileId, '%s\n', fileNameHintRuler);

fprintf(logFileId, '%s\n', ['Started at: ', dataTimeStrStart]);
fprintf(logFileId, '%s\n', ['Host name: ', curHostname]);
fprintf(logFileId, '%s\n', ['Num of workers: ', ...
    num2str(numOfAvailableWorkers)]);

fprintf(logFileId, '\n');
fprintf(logFileId, '%s\n', ['Num of TX towers involved: ', ...
    num2str(numOfCellAntsObserved)]);
fprintf(logFileId, '%s\n', ['Num of RX heights inspected: ', ...
    num2str(numOfHsInspected)]);

fprintf(logFileId, '\n');
fprintf(logFileId, '%s\n', ...
    ['Num of individual maps generated: ', ...
    num2str(numOfMapsGenerated)]);
fprintf(logFileId, '%s\n', ...
    ['Num of pixels (path loss values) evaluated: ', ...
    num2str(numOfPixelsCovered)]);
fprintf(logFileId, '\n');

if flagProcessInterrupted
    fprintf(logFileId, '%s\n', ...
        'Current session was resumed from history progress data.');
else
    % Collect execution time information.
    execTimeInSec = toc(timerValueStart);
    aveTimePerMapInSec = execTimeInSec./numOfMapsGenerated;
    aveTimePerPixInSec = execTimeInSec./numOfPixelsCovered;
    
    
    fprintf(logFileId, '%s\n', ['Total time used: ', ...
        num2str(execTimeInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['    (', ...
        seconds2human(execTimeInSec), ')']);
    fprintf(logFileId, '%s\n', ['Average time used per map: ', ...
        num2str(aveTimePerMapInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['    (', ...
        seconds2human(aveTimePerMapInSec), ')']);
    fprintf(logFileId, '%s\n', ['Average time used per pixel: ', ...
        num2str(aveTimePerPixInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['    (', ...
        seconds2human(aveTimePerPixInSec), ')']);
end

for idxH = 1:numOfHsInspected
    curHExecTimeInSec = execTimeInSecForAllHs(idxH);
    curHAveTimePerMapInSec = curHExecTimeInSec./numOfCellAntsObserved;
    curHAveTimePerPixInSec ...
        = curHExecTimeInSec./numOfPixelsCoveredForHs(idxH);
    
    fprintf(logFileId, '\n%s\n', ['For RX height # ', num2str(idxH), ':']);
    fprintf(logFileId, '%s\n', ['    Total time used: ', ...
        num2str(curHExecTimeInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['        (', ...
        seconds2human(curHExecTimeInSec), ')']);
    fprintf(logFileId, '%s\n', ['    Average time used per map: ', ...
        num2str(curHAveTimePerMapInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['        (', ...
        seconds2human(curHAveTimePerMapInSec), ')']);
    fprintf(logFileId, '%s\n', ['    Average time used per pixel: ', ...
        num2str(curHAveTimePerPixInSec), ' seconds']);
    fprintf(logFileId, '%s\n', ['        (', ...
        seconds2human(curHAveTimePerPixInSec), ')']);
end

fprintf(logFileId, '%s\n', fileNameHintRuler);
fclose(logFileId);

% EOF