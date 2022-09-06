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

if simConfigs.FLAG_EVAL_LOSS_THROUGH_VEG
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
end

% More records for ACRE LoRaWAN simulations: lon, lat, install_height_m,
% ehata_path_loss_db, accumulate_blockage_dist_m.
lat = simState.mapGridLatLonPts(:,1);
lon = simState.mapGridLatLonPts(:,2);
for curIdxH = 1:length(rxAntHs)
    ehata_path_loss_db = simState.coverageMaps{curIdxH};
    accumulate_blockage_dist_m = simState.blockageDistMaps{curIdxH};
    if isfield(simState, 'coverageItmMaps')
        itm_path_loss_db = simState.coverageItmMaps{curIdxH};
    end

    rxAntH = rxAntHs(curIdxH);
    curPathToSaveCsv = fullfile(folderToSaveCsvs, ...
        ['latLonHPathLossBlockDist_RxHeight_', num2str(rxAntH), '.csv']);

    if ~exist(curPathToSaveCsv, 'file')
        install_height_m = ones(size(lat)).*rxAntH;

        if isfield(simState, 'coverageItmMaps')
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

    % Generate a few overview figures.
    curPathToSaveFig = fullfile(folderToSaveCsvs, ...
        ['accumulate_blockage_dist_m-', num2str(rxAntH), '.csv']);
    if ~exist([curPathToSaveFig, '.jpg'], 'file')
        figure; plot3k([lon, lat, accumulate_blockage_dist_m], ...
            'Labels', {'Accumulate Blockage Distance (m)', ...
            'Longitude (degree)', 'Latitude (degree)', '', ''});
        view(2); plot_google_map('MapType', 'hybrid');
        saveas(gcf, [curPathToSaveFig, '.fig']);
        saveas(gcf, [curPathToSaveFig, '.jpg']);
    end

    curPathToSaveFig = fullfile(folderToSaveCsvs, ...
        ['ehata_path_loss_db', num2str(rxAntH), '.csv']);
    if ~exist([curPathToSaveFig, '.jpg'], 'file')
        figure; plot3k([lon, lat, ehata_path_loss_db], 'Labels', ...
            {'eHata Path Loss (dB)', ...
            'Longitude (degree)', 'Latitude (degree)', '', ''});
        view(2); plot_google_map('MapType', 'hybrid');
        saveas(gcf, [curPathToSaveFig, '.fig']);
        saveas(gcf, [curPathToSaveFig, '.jpg']);
    end

    if isfield(simState, 'coverageItmMaps')
        curPathToSaveFig = fullfile(folderToSaveCsvs, ...
            ['itm_path_loss_db', num2str(rxAntH), '.csv']);
        if ~exist([curPathToSaveFig, '.jpg'], 'file')
            figure; plot3k([lon, lat, itm_path_loss_db], 'Labels', ...
                {'ITM Path Loss (dB)', ...
                'Longitude (degree)', 'Latitude (degree)', '', ''});
            view(2); plot_google_map('MapType', 'hybrid');
            saveas(gcf, [curPathToSaveFig, '.fig']);
            saveas(gcf, [curPathToSaveFig, '.jpg']);
        end
    end
end

% Export the geo info for the grid if that is available.
if isfield(simState, 'mapGridGroundEleInM')
    curPathToSaveCsv = fullfile(folderToSaveCsvs, ...
        'mapGridGeoInfo.csv');

    ele_usgs_m = simState.mapGridGroundEleInM;
    lidar_z_m = simState.mapGridLiarZInM;

    if ~exist(curPathToSaveCsv, 'file')
        writetable(table(lon, lat, ele_usgs_m, lidar_z_m), ...
            curPathToSaveCsv);
    end

    % Generate a few overview figures.
    curPathToSaveFig = fullfile(folderToSaveCsvs, 'ele_usgs_m');
    if ~exist([curPathToSaveFig, '.jpg'], 'file')
        figure; plot3k([lon, lat, ele_usgs_m], 'Labels', ...
            {'Bare-Earth Ground Elevation (m)', ...
            'Longitude (degree)', 'Latitude (degree)', '', ''});
        view(2); plot_google_map('MapType', 'hybrid');
        saveas(gcf, [curPathToSaveFig, '.fig']);
        saveas(gcf, [curPathToSaveFig, '.jpg']);
    end

    curPathToSaveFig = fullfile(folderToSaveCsvs, 'lidar_z_m');
    if ~exist([curPathToSaveFig, '.jpg'], 'file')
        figure; plot3k([lon, lat, lidar_z_m], 'Labels', ...
            {'LiDAR z (m)', ...
            'Longitude (degree)', 'Latitude (degree)', '', ''});
        view(2); plot_google_map('MapType', 'hybrid');
        saveas(gcf, [curPathToSaveFig, '.fig']);
        saveas(gcf, [curPathToSaveFig, '.jpg']);
    end

    curPathToSaveFig = fullfile(folderToSaveCsvs, 'lidar_z_m-ele_usgs_m');
    if ~exist([curPathToSaveFig, '.jpg'], 'file')
        figure; plot3k([lon, lat, lidar_z_m-ele_usgs_m], 'Labels', ...
            {'LiDAR z - Bare-Earth Ground Elevation (m)', ...
            'Longitude (degree)', 'Latitude (degree)', '', ''});
        view(2); plot_google_map('MapType', 'hybrid');
        saveas(gcf, [curPathToSaveFig, '.fig']);
        saveas(gcf, [curPathToSaveFig, '.jpg']);
    end
end

close all;
flagSuccess = true;

disp('    Done!')
end
% EOF