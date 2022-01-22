function [ flagSuccess ] = exportSimTimeToCsv(simTime, ...
    presets, carrierFreqsInMHz, pathToSaveCsv, roundNumOfDigit)
%EXPORTSIMTIMETOCSV Save simulation time values to a human-readable .csv
%file.
%
% Example content of the .csv file:
%
%   Carrier\Scenario	Tipp	WHIN	Indiana
%    1900MHz
%   3700MHz
%    4700MHz
%   Ave Time
%
% Inputs:
%   - simTime
%     A matrix for the simulation time values. Each row corresponds to the
%     time values for all presets with one carrier frequency.
%   - presets, carrierFreqsInMHz
%     Row cells of the presets and carrier frequencies, respectively.
%   - pathToSaveCsv
%     Full path to save the .csv file. Please make sure the parent folder
%     already exists before calling exportSimTimeToCsv.
%   - roundNumOfDigit
%     Optional. If present, the output values will be rounded to
%     roundNumOfDigit digits. More details can be found in function round.
%
% Yaguang Zhang, Purdue, 01/22/2022

flagSuccess = false; %#ok<NASGU>

aveSimTime = mean(simTime);
if exist('roundNumOfDigit', 'var')
    simTime = round(simTime, roundNumOfDigit);
    aveSimTime = round(aveSimTime, roundNumOfDigit);
end

% Start with the first column.
firstColTable = cell2table( ...
    (cellfun(@(f) [num2str(f), 'MHz'], carrierFreqsInMHz, ...
    'UniformOutput', false))', ...
    'VariableNames', {'Carrier\Scenario'});

% Append the table with simulation time values.
simTimeTable = [firstColTable, ...
    array2table(simTime, 'VariableNames', presets)];

% Export the table to a .csv file.
writetable(simTimeTable, pathToSaveCsv);

% Append the .csv file with an average-time row.
writecell(['Ave Time', ...
    arrayfun(@(t) num2str(t), aveSimTime, ...
    'UniformOutput', false)], ...
    pathToSaveCsv, 'WriteMode', 'append');

flagSuccess = true;
end
% EOF