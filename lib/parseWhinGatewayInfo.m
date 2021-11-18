function [ pathToSaveWhinGatewayLocs ] = parseWhinGatewayInfo( ...
    ABS_PATH_TO_SHARED_FOLDER, dirToWhinGatewayLog)
%PARSEWHINGATEWAYINFO A helper function to extract the WHIN gateway
%information from the raw log file and convert it to a .csv file compliant
%with the cellular coverage simulator.
%
% Inputs:
%   - ABS_PATH_TO_SHARED_FOLDER
%     Full path to the cellular coverage analysis project. Please refer to
%     setPath.m for more details.
%   - dirToWhinGatewayLog
%     The full path to the WHIN gateway information file. If empty, a
%     default path will be provided.
%
% Output:
%   - pathToSaveWhinGatewayLocs
%     The full path to the resultant .csv file.
%
% Yaguang Zhang, Purdue, 11/17/2021

% Parameters
if (~exist('dirToWhinGatewayLog', 'var')) ...
        || isempty(dirToWhinGatewayLog)
    dirToWhinGatewayLog = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
        'CellTowerInfo', 'WhinWeatherStations', ...
        'lora_tower_gateway_positions.csv');
end
[folderToSaveWhinGatewayLocs, ~, ~] = fileparts(dirToWhinGatewayLog);
pathToSaveWhinGatewayLocs = fullfile(folderToSaveWhinGatewayLocs, ...
    'WhinLoraWanGateways.csv');

% Only process the raw log if it is not yet done.
if ~exist(pathToSaveWhinGatewayLocs, 'file')
    whinGatewayTb = readtable(dirToWhinGatewayLog);

    Lat = whinGatewayTb.latitude;
    Long = whinGatewayTb.longitude;
    HeightInM = whinGatewayTb.altitude_above_ground_meters;
    Site = (1:length(Lat))';

    writetable(table(Site, Lat, Long, HeightInM), ...
        pathToSaveWhinGatewayLocs)
end

% EOF