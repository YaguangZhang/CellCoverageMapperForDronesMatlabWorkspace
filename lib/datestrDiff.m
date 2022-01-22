function [deltaDatetime] = datestrDiff(dateStr1, dateStr2, inputFormat)
%DATESTRDIFF Compare the input date strings and output the difference
%(dateStr2 - dateStr1) as a datetime object.
%
% By default, inputFormat = 'yyyy/MM/dd HH:mm:ss', corresponding to the
% output format datetimeFormat = 'yyyy/mm/dd HH:MM:ss' for function datastr
% defined in prepareSimulationEnv. This function could be very helpful in
% getting the simulation time from the diary logs.
%
% Yaguang Zhang, Purdue, 01/22/2022

if ~exist('inputFormat', 'var')
    inputFormat = 'yyyy/MM/dd HH:mm:ss';
end

deltaDatetime = datetime(dateStr2, 'InputFormat', inputFormat) ...
    - datetime(dateStr1, 'InputFormat', inputFormat);
end
% EOF