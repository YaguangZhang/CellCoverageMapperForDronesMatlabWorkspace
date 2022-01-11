function [ flagSuccess ] = exportMaps( ...
    pathToSaveResults, simConfigs, simState)
%EXPORTMAPS A helper function to export the raw data from selected maps
%(stored in simState) to .csv files.
%
% Yaguang Zhang, Purdue, 01/11/2022

disp(' ')
disp('    Exporting raw data for selected maps ...')

flagSuccess = false; %#ok<NASGU>

folderToSaveCsvs = fullfile(pathToSaveResults, 'MapData');
if ~exist(folderToSaveCsvs, 'dir')
    mkdir(folderToSaveCsvs);
end
rxAntHs = simConfigs.RX_ANT_HEIGHTS_TO_INSPECT_IN_M;

% We only need the 1.5 m case for now.
%   for curIdxH = 1:length(rxAntHs)
for curIdxH = 1
    rxAntH = rxAntHs(curIdxH);
    curPathToSaveCsv = fullfile(folderToSaveCsvs, ...
        ['latLonPathLossWithVeg_RxHeight_', num2str(rxAntH), '.csv']);

    if ~exist(curPathToSaveCsv, 'file')
        Lat = simState.mapGridLatLonPts(:,1);
        Lon = simState.mapGridLatLonPts(:,2);
        PathLossWithVegInDb = simState.pathLossWithVegMaps{curIdxH};
    end
    writetable(table(Lat, Lon, PathLossWithVegInDb), ...
        curPathToSaveCsv)
end
flagSuccess = true;

disp('    Done!')
end
% EOF