% INVESTIGATEGROUNDROUGHNESS Investigate different approaches of testing
% ground roughness.
%
% Yaguang Zhang, Purdue, 11/21/2019

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..'); addpath('lib');
curFileName = mfilename;

prepareSimulationEnv;

%% Script Parameters

ABS_PATH_TO_PURDUE_FARM_LIDAR = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'Lidar', 'PurdueFarms');

% For interpolating the LiDAR data.
lidarGridResultionInM = 1;

% For the sample density of the ground roughness results.
numberOfPixelsOnLongerSide = 100;
regionalCircleRadiaToInspectInM = [5, 10, 50, 100];

% For plotting.
flagGenFigSilently = true;

pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', '7_GroundRoughnessInvestigation');

% Create directories if necessary.
if exist(pathToSaveResults, 'dir')~=7
    mkdir(pathToSaveResults);
end

%% Locate LiDAR Data

disp(' ')
disp('    Locating LiDAR datasets ...')

% Locate the LiDAR data folders for the Purdue research farms.
purdueFramLidarSets = dir(ABS_PATH_TO_PURDUE_FARM_LIDAR);
purdueFramLidarSets = purdueFramLidarSets([purdueFramLidarSets.isdir]);
% Remove '.' and '..'.
purdueFramLidarSets = purdueFramLidarSets( ...
    ~ismember({purdueFramLidarSets.name}, {'.','..'}));
numOfLidarSets = length(purdueFramLidarSets);

disp('    Done!')

%% Investigate the Local RMSD Approach

disp(' ')
disp('    Computing local RMSDs ...')

for idxLidarSet = 1:numOfLidarSets
    % Load LiDAR data.
    curDataSetLabel = purdueFramLidarSets(idxLidarSet).name;
    disp(' ')
    disp(['        Loading dataset for ', curDataSetLabel, ' ...'])
    
    curDataSetDir = fullfile(purdueFramLidarSets(idxLidarSet).folder, ...
        curDataSetLabel, 'points.las');
    curDataSetLatLonRangeDir ...
        = fullfile(purdueFramLidarSets(idxLidarSet).folder, ...
        curDataSetLabel, 'LatLonRangeUsed.txt');
    
    [lidarXYZ, hFigGoogleMap, hFigLidarData] ...
        = fetchLasLidarDateset(curDataSetDir, curDataSetLatLonRangeDir, ...
        lidarGridResultionInM, ~flagGenFigSilently);
    saveas(hFigGoogleMap, fullfile(pathToSaveResults, ...
        [curDataSetLabel, '_MapWithAreaOfInterest.png']));
    saveas(hFigLidarData, fullfile(pathToSaveResults, ...
        [curDataSetLabel, '_Lidar.png']));
    
    disp('        Done!')
    
    disp(' ')
    disp('        Computing local RMSDs ...')
    for regionalCircleRadiusInM = regionalCircleRadiaToInspectInM
        disp(['            For regionalCircleRadiusInM = ', ...
            num2str(regionalCircleRadiusInM), ' ...'])
        
        % Cache the results for speed.
        curAbsPathToCacheRegRmsd = fullfile(pathToSaveResults, ...
            [curDataSetLabel, '_RegRmsd_regionalCircleRadiusInM_', ...
            strrep(num2str(regionalCircleRadiusInM), '.', '_'), '.mat']);
        if exist(curAbsPathToCacheRegRmsd, 'file')
            load(curAbsPathToCacheRegRmsd);
            hFigRegRmsdWrtMean = surfWithContour( ...
                [gridXYPts, regRmsdsWrtMean], ...
                {'', '', '', '', 'Regional RMSD (m)', ''}, ...
                ~flagGenFigSilently);
        else
            [hFigRegRmsdWrtMean, gridXYPts, regRmsdsWrtMean] ...
                = plotRegionalRmsdWrtMean(...
                lidarXYZ, regionalCircleRadiusInM, ...
                numberOfPixelsOnLongerSide, ~flagGenFigSilently);
            
            save(curAbsPathToCacheRegRmsd, 'gridXYPts', 'regRmsdsWrtMean');
        end
        
        curFileName = [curDataSetLabel, '_RegCirRInM_', ...
            num2str(regionalCircleRadiusInM)];
        saveas(hFigRegRmsdWrtMean, fullfile(pathToSaveResults, ...
            [curFileName, '_RegRmsdWrtMean.png']));
        
        disp('            Done!')
    end
    
    if flagGenFigSilently
        close all;
    end
    
    disp('        Done!')
end

disp('    Done!')

%% Empirical CDFs

disp(' ')
disp('    Plotting empirical CDFs ...')

markersToUse = {''}; % {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', 'p'};
linesToUse = {'-', '--', ':', '-.'};
numOfMarkers = length(markersToUse);
numOfLines = length(linesToUse);

for regionalCircleRadiusInM = regionalCircleRadiaToInspectInM
    hFigEmpCdf = figure('Visible', ~flagGenFigSilently);
    hold on; legendHs = []; legendsToSet = '';
    
    for idxLidarSet = 1:numOfLidarSets
        curDataSetLabel = purdueFramLidarSets(idxLidarSet).name;
        
        % Loading the cached regional RMSDs.
        curAbsPathToCacheRegRmsd = fullfile(pathToSaveResults, ...
            [curDataSetLabel, '_RegRmsd_regionalCircleRadiusInM_', ...
            strrep(num2str(regionalCircleRadiusInM), '.', '_'), ...
            '.mat']);
        
        histRegRmsdResults = load(curAbsPathToCacheRegRmsd);
        curRegRmsds = histRegRmsdResults.regRmsdsWrtMean;
        [curCdfValues, curCdfXs] = ecdf(curRegRmsds(:));
        
        legendsToSet = strcat(legendsToSet, ...
            ", '", curDataSetLabel, "'");
        markerToUse = [markersToUse{mod(idxLidarSet,numOfMarkers)+1}, ...
            linesToUse{mod(idxLidarSet,numOfLines)+1}];
        curHCdf = plot(curCdfXs, curCdfValues, markerToUse, ...
            'LineWidth', 1);
        legendHs(end+1) = curHCdf; %#ok<SAGROW>
    end
    
    % Remove the starting ', '.
    legendsToSet = char(legendsToSet);
    legendsToSet = legendsToSet(3:end);
    axis tight; grid on; grid minor;
    set(gca, 'XScale', 'log');
    eval(['legend(legendHs,', legendsToSet, ', "Location", "SouthEast");']);
    xlabel('Regional RMSD (m)'); ylabel('Empirical CDF');
    
    curFileName = ['EmpCdfForRegRmsd_RegCirRInM_', ...
        num2str(regionalCircleRadiusInM)];
    saveas(hFigEmpCdf, fullfile(pathToSaveResults, ...
        [curFileName, '_RegRmsdWrtMean.png']));
end

if flagGenFigSilently
    close all;
end

disp('    Done!')

% EOF