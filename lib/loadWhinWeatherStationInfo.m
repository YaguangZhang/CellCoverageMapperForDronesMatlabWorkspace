function [ lats, lons, intIds, names ] = loadWhinWeatherStationInfo()
%LOADWHINWEATHERSTATIONINFO A helper function to extract the WHIN weather
%station information from the raw log file.
%
% Outputs:
%   - lats, lons
%     GPS (lats, lons) for the weather stations.
%   - intIds
%     Integer IDs of the weather stations.
%   - names
%     A column cell for the string names of the weather stations.
%
% Yaguang Zhang, Purdue, 11/18/2021

FLAG_SORT_RECORDS = false;

% Parameters
ABS_PATH_TO_SHARED_FOLDER = evalin('base', 'ABS_PATH_TO_SHARED_FOLDER');
dirToWhinWeatherStationLog = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'CellTowerInfo', 'WhinWeatherStations', ...
    'weather_station_locations.csv');

% Load the location information.
whinWeatherStationTb = readtable(dirToWhinWeatherStationLog);

if FLAG_SORT_RECORDS
    % Sort the records by station ID.
    whinWeatherStationTb = sortrows(whinWeatherStationTb, ...
        'id'); %#ok<UNRCH>
    assert(all(whinWeatherStationTb.id ...
        ==(1:length(whinWeatherStationTb.id))'), ...
        'Weather station data set is not complete!');
end

lats = whinWeatherStationTb.latitude;
lons = whinWeatherStationTb.longitude;
intIds = whinWeatherStationTb.id;
names = whinWeatherStationTb.name;

end
% EOF