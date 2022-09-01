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

        writetable(table(Lat, Lon, PathLossWithVegInDb), ...
            curPathToSaveCsv);
    end
end

% More records for ACRE LoRaWAN simulations: lon, lat, install_height_m,
% ehata_path_loss_db, accumulate_blockage_dist_m.
for curIdxH = 1:length(rxAntHs)
    rxAntH = rxAntHs(curIdxH);
    curPathToSaveCsv = fullfile(folderToSaveCsvs, ...
        ['latLonHPathLossBlockDist_RxHeight_', num2str(rxAntH), '.csv']);

    if ~exist(curPathToSaveCsv, 'file')
        lat = simState.mapGridLatLonPts(:,1);
        lon = simState.mapGridLatLonPts(:,2);
        install_height_m = ones(size(lat)).*rxAntH;
        ehata_path_loss_db = simState.coverageMaps{curIdxH};
        accumulate_blockage_dist_m = simState.blockageDistMaps{curIdxH};

        if isfield(simState, 'coverageMapsItm')
            itm_path_loss_db = simState.coverageMapsItm{curIdxH};

            writetable(table(lon, lat, install_height_m, ...
                ehata_path_loss_db, itm_path_loss_db, ...
                accumulate_blockage_dist_m), ...
                curPathToSaveCsv);
        else
            writetable(table(lon, lat, install_height_m, ...
                ehata_path_loss_db, accumulate_blockage_dist_m), ...
                curPathToSaveCsv);
        end
    end
end

% Export the geo info for the grid if that is available.
if isfield(simState, 'mapGridGroundEleInM')
    curPathToSaveCsv = fullfile(folderToSaveCsvs, ...
        'mapGridGeoInfo.csv');

    ele_usgs_m = simState.mapGridGroundEleInM;
    lidar_z_m = simState.mapGridLiarZInM;

    writetable(table(lon, lat, ele_usgs_m, lidar_z_m), curPathToSaveCsv);
end

flagSuccess = true;

disp('    Done!')
end
% EOF