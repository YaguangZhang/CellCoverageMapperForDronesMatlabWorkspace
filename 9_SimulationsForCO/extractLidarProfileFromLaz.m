%EXTRACTLIDARPROFILEFROMLAZ Snippet to extract a LiDAR profile for a long
%link (~14.5 km) in Colorado based on .laz file(s).
%
% Developed with Matlab R2022a.
%
% Yaguang Zhang, Purdue, 10/19/2022

clear; clc; close all; dbstop if error;

% Locate the Matlab workspace and save the current filename.
cd(fileparts(mfilename('fullpath'))); cd('..');
addpath('lib'); addpath('.');
curFileName = mfilename;

prepareSimulationEnv;

%% Configurations

% Path to the folder holding the LiDAR laz file(s) to use. Currently, we
% have datasets:
%   - 'OneTileExample'
%      One example tile from Anemo fetched via the USGS TNM-LPC API.
%   - 'Complete'
%      A complete dataset covering the area needed for the long link, with
%      some redundancy: there are tiles covering the same region, but
%      probably are from different datasets.
%   - 'Selected'
%      Selected tiles from the complete set which are enough for the long
%      link.
lidarPreset = 'Selected';

switch lidarPreset
    case {'OneTileExample', 'Complete', 'Selected'}
        pathToLazFolder = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
            'Lidar', 'CO', 'YZ_CA_CO_LiDAR', lidarPreset);
    otherwise
        error(['Unsupported LiDAR dataset: ', lidarPreset, '!'])
end

% Link of interest:
latLonStart = [39.991794, -105.274711];
latLonEnd = [40.121142, -105.250972];

% The absolute path to the folder for saving the results.
pathToSaveResults = fullfile(ABS_PATH_TO_SHARED_FOLDER, ...
    'PostProcessingResults', ...
    'CO_LongLink_LidarProfileForProfA', lidarPreset);

if ~exist(pathToSaveResults, 'dir')
    mkdir(pathToSaveResults)
end

maxAllowedResInM = 1;
minNumSamps = 10;

%% Preprocess LiDAR Data
% We will generate for each .laz file a .mat meta file with the same file
% name containing fields:
%   - lazPath
%     Absolute path to the corresponding .laz file.
%   - lonLatBoundary
%     A 2-column matrix with each row containing the (longitude, latitude)
%     of a vertex of the boundary polygon. The first row and the last row
%     should be the same to indicate a closed polygon.
%   - crs
%     The coordinate reference system object for projfwd and projinv.
%   - xYBoundary
%     The (x, y) boundary from converting lonLatBoundary with crs.

% Find all laz files.
dirsLaz = rdir(fullfile(pathToLazFolder, '*.laz'));
numOfLazs = length(dirsLaz);

% Extract the coverage boundaries.
lonLatBounds = cell(numOfLazs, 1);
crses = cell(numOfLazs, 1);
% parfor may cause use too much resource because of the LiDAR plots.
for idxLaz = 1:numOfLazs
    curDirLaz = dirsLaz(idxLaz);
    [~, lazFileName] =  fileparts(curDirLaz.name);

    curPathToLazMeta = fullfile(curDirLaz.folder, [lazFileName, '.mat']);
    if exist(curPathToLazMeta, 'file')
        curMeta = load(curPathToLazMeta, 'lonLatBoundary', 'crs');
        lonLatBounds{idxLaz} = curMeta.lonLatBoundary;
        crses{idxLaz} = curMeta.crs;
    else
        lazPath = curDirLaz.name;

        lasReader = lasFileReader(lazPath);
        ptCloud = readPointCloud(lasReader);
        crs = readCRS(lasReader);

        % For simplicity, we will use a rectangle as the boundary.
        %     minX = min(ptCloud.Location(:, 1));
        %      maxX = max(ptCloud.Location(:, 1));
        %     minY = min(ptCloud.Location(:, 2));
        %      maxY = max(ptCloud.Location(:, 2));

        % We will expand the boundary by half of the resolution.
        spatialRes = min( ...
            min(abs( diff(unique(ptCloud.Location(:, 1))) )), ...
            min(abs( diff(unique(ptCloud.Location(:, 2))) )));
        distToExp = spatialRes/2;

        minX = ptCloud.XLimits(1) - distToExp;
        maxX = ptCloud.XLimits(2) + distToExp;
        minY = ptCloud.YLimits(1) - distToExp;
        maxY = ptCloud.YLimits(2) + distToExp;
        xYBoundary = [minX, minY; minX, maxY; maxX, maxY; maxX, minY; ...
            minX, minY];

        % Convert UTM (x, y) to (lat, lon).
        [latsBound, lonsBound] = projinv(crs, ...
            xYBoundary(:, 1), xYBoundary(:, 2));
        lonLatBoundary = [lonsBound, latsBound];

        % Also generate an overview plot of the LiDAR data.
        [lidarLats, lidarLons] = projinv(crs, ...
            ptCloud.Location(:, 1), ptCloud.Location(:, 2));

        curFigNamePrefix = fullfile(curDirLaz.folder, lazFileName);

        hOverviewTile = figure;
        plot3k([lidarLons, lidarLats, ptCloud.Location(:, 3)], ...
            'Labels', {'Anemo LiDAR Data Overview', ...
            'Longitude (degree)', 'Latitude (degree)', 'LiDAR z (m)', ...
            'LiDAR z (m)'});
        view(2);
        plot_google_map('MapType', 'hybrid');

        saveas(hOverviewTile, [curFigNamePrefix, '.jpg']);
        close(hOverviewTile);

        hOverviewTile = figure;
        maxLidarZInM = max(ptCloud.Location(:, 3));
        hLink = plot3([latLonStart(2), latLonEnd(2)], ...
            [latLonStart(1), latLonEnd(1)], ...
            [maxLidarZInM, maxLidarZInM], '-r'); % #ok<PFBNS>
        view(2);
        plot_google_map('MapType', 'hybrid');
        plot3k([lidarLons, lidarLats, ptCloud.Location(:, 3)], ...
            'Labels', {'Anemo LiDAR Data Overview', ...
            'Longitude (degree)', 'Latitude (degree)', 'LiDAR z (m)', ...
            'LiDAR z (m)'});
        view(2);
        legend(hLink, 'Link of Interest');

        saveas(hOverviewTile, [curFigNamePrefix, '_WithExampleLink.jpg']);
        close(hOverviewTile);

        % Save the results.
        parsave(curPathToLazMeta, lazPath, lonLatBoundary, ...
            crs, xYBoundary);

        lonLatBounds{idxLaz} = lonLatBoundary;
        crses{idxLaz} = crs;
    end
end

% Choose a common crs to use by voting.
uniqueCrses = {};
uniqueCrsesCnts = [];
for idxCrs = 1:numOfLazs
    curCrs = crses{idxCrs};

    idxMatchingUniCrs = nan;
    if ~isempty(uniqueCrses)
        idxMatchingUniCrs = nan;
        for idxUniCrs = 1:length(uniqueCrses)
            if strcmp(uniqueCrses{idxUniCrs}.Name, curCrs.Name)
                % CRS already registered.
                idxMatchingUniCrs = idxUniCrs;
                break;
            else
                [~, d1, d2] = comp_struct( ...
                    uniqueCrses{idxUniCrs}.ProjectionParameters, ...
                    curCrs.ProjectionParameters);
                if isempty(d1) && isempty(d2) ...
                        && strcmp( ...
                        uniqueCrses{idxUniCrs}.ProjectionMethod, ...
                        curCrs.ProjectionMethod) ...
                        && strcmp(uniqueCrses{idxUniCrs}.LengthUnit, ...
                        curCrs.LengthUnit)
                    % CRS already registered.
                    idxMatchingUniCrs = idxUniCrs;
                    break;
                end
            end
        end
    end

    if isnan(idxMatchingUniCrs)
        uniqueCrses{end+1, 1} = curCrs; %#ok<SAGROW>
        uniqueCrsesCnts(end+1, 1) = 1;  %#ok<SAGROW>
    else
        uniqueCrsesCnts(idxMatchingUniCrs) ...
            = uniqueCrsesCnts(idxMatchingUniCrs) + 1; %#ok<SAGROW>
    end
end

[~, maxCntIdx] = max(uniqueCrsesCnts);
bestCrs = uniqueCrses{maxCntIdx};

%% Construct Link of Interest in UTM

[utmXStart, utmYStart] = projfwd(bestCrs, latLonStart(1), latLonStart(2));
[utmXEnd, utmYEnd] = projfwd(bestCrs, latLonEnd(1), latLonEnd(2));

profLocsXY = generateTerrainProfileSampLocs( ...
    [utmXStart, utmYStart], [utmXEnd, utmYEnd], ...
    maxAllowedResInM, minNumSamps);

[profLocsLat, profLocsLon] = projinv(bestCrs, ...
    profLocsXY(:, 1), profLocsXY(:, 2));
profLocsLonLat = [profLocsLon, profLocsLat];

%% Overview Figure of the LiDAR Coverage

hOverviewLidarCov = figure; hold on;
for idxTile = 1:numOfLazs
    hCov = patch('XData', lonLatBounds{idxTile}(:, 1), ...
        'YData', lonLatBounds{idxTile}(:, 2), ...
        'FaceColor', 'b', 'FaceAlpha',.5);
end
hLink = plot([latLonStart(2), latLonEnd(2)], ...
    [latLonStart(1), latLonEnd(1)], '-r');
hProfLocs = plot(profLocsLonLat(:, 1), profLocsLonLat(:, 2), '.r');
plot_google_map('MapType', 'hybrid');
legend([hLink, hCov], 'Link of Interest', 'LiDAR Coverage');

curFigFilename = fullfile(pathToSaveResults, ...
    'LidarOverview_Anemo');
saveas(hOverviewLidarCov, [curFigFilename, '.jpg']);

% Zoomed-in versions.
axis([-105.2585011053211, -105.2574053433766, ...
    40.0826508226302, 40.0833120621302]);
saveas(hOverviewLidarCov, [curFigFilename, '_Zoom_1.jpg']);
axis([-105.2579914369130, -105.2579637175118, ...
    40.0829797815678, 40.0829965088883]);
saveas(hOverviewLidarCov, [curFigFilename, '_Zoom_2.jpg']);
axis([-105.2579747017336, -105.2579741493096, ...
    40.0829873664333, 40.0829877085021]);
saveas(hOverviewLidarCov, [curFigFilename, '_Zoom_3.jpg']);

close(hOverviewLidarCov);

%% Extract Profile

% Find which tile is needed for each sampled location along the link.
numOfProfSamps = size(profLocsLonLat, 1);
sampIndicesForTile = nan(numOfLazs, 1);
tileIndicesForSamp = nan(numOfProfSamps, 1);
profZs = nan(numOfProfSamps, 1);

for idxTile = 1:numOfLazs    
    crs = crses{idxTile};
    [curProfLocsX, curProfLocsY] = projfwd(crs, ...
        profLocsLonLat(:, 2), profLocsLonLat(:, 1));
    [curBoundX, curBoundY] = projfwd(crs, ...
        lonLatBounds{idxLaz}(:, 2), lonLatBounds{idxLaz}(:, 1));

    curSampsInTile = inpoly2( ...
        [curProfLocsX, curProfLocsY], [curBoundX, curBoundY]);

    if ~isempty(curSampsInTile)
        assert(all(isnan(tileIndicesForSamp(curSampsInTile))), ...
            'Tile already found!');

        % Load current tile.
        curDirLaz = dirsLaz(idxLaz);
        lazPath = curDirLaz.name;

        lasReader = lasFileReader(lazPath);
        ptCloud = readPointCloud(lasReader);

        % Fetch LiDAR z by nearest neighbor.
        %     curProfZs = ptCloud.Location(dsearchn( ...
        %         [ptCloud.Location(:, 1), ptCloud.Location(:, 2)], ...
        %          [curProfLocsX(curSampsInTile), ...
        %         curProfLocsY(curSampsInTile)] ...
        %          ), 3);
        %     curProfZs = interp2( ...
        %         ptCloud.Location(:, 1), ptCloud.Location(:, 2), ...
        %          ptCloud.Location(:, 3), ...
        %         curProfLocsX(curSampsInTile), ...
        %          curProfLocsY(curSampsInTile), ...
        %         'nearest');
        F = scatteredInterpolant( ...
            ptCloud.Location(:, 1), ptCloud.Location(:, 2), ...
            ptCloud.Location(:, 3), 'nearest');

        curProfZs = F(curProfLocsX(curSampsInTile), ...
            curProfLocsY(curSampsInTile));

        assert(all(~isnan(curProfZs)), 'Invalid LiDAR z!');

        profZs(curSampsInTile) = curProfZs;

        % Overview figure of this profile segment.
        hOverviewProfSeg = figure;
        plot3k([profLocsLonLat(curSampsInTile, :), curProfZs]);
        view(2);
        plot_google_map('MapType', 'hybrid');

        curFigFilename = fullfile(pathToSaveResults, ...
            ['ProfSeg_LidarTileIdx_', num2str(idxTile)]);
        saveas(hOverviewProfSeg, [curFigFilename, '.jpg']);
    end
end

%% Save Results

save(fullfile(pathToSaveResults, 'LidarProfile.mat'), ...
    'profLocsLonLat', 'profZs');

%% Overview Figure of the Profile

% EOF