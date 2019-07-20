function [ medianBPLMap, medianBPLMapXLabels,  medianBPLMapYLabels] ...
    = genMedianBasicPropLossMaps(fsMHz, ...
    baseAntXY, baseAntHeightInM, ...
    rxLocXs, rxLocYs, mobileAntHeightInM, ...
    terrainXYZs, lidarXYZs, region, ...
    eleProfileResForEHataInM, libraryToUse, NTIA_EHATA_RELIABILITY, ...
    numOfWorkers, EHATA_VALID_DIST_RANGE_IN_M)
%GENMEDIANBASICPROPLOSSMAPS Compute the median basic propogation loss via
%the Extended Hata model (when appropriate) and the free-space path loss
%model for a grid to form a map.
%
% Inputs:
%   - fsMHz
%     The carrier frequency in MHz.
%   - baseAntXY, baseAntHeightInM
%     The basestation location (x,y) in the UTM system and its antenna's
%     height in meter, respectively.
%   - rxLocXs, rxLocYs
%     The x and y labels in the UTM system, respectively, for the grid to
%     form the map of the area of interest. All combinations of their
%     elements are the receiver's (i.e. the mobile station's) locations.
%   - mobileAntHeightInM
%     The antenna height of the mobile station in meter.
%   - terrainXYZs, lidarXYZs
%     The terrain and LiDAR information for the area of interest,
%     respectively. Each of them should be:
%       - A matrix representing the terrain elevation/LiDAR information in
%         the form of [xs, ys, zs], where xs and ys are column vectors
%         reprenting the UTM locations where the elevation information is
%         available, and zs is a column vector storing the corresponding
%         elevations/LiDAR zs in meter. Or,
%       - A function to fetch the elevation/LiDAR z for given locations,
%         e.g. an instance 'scatteredInterpolant' generated by
%         scatteredInterpolant from the terrain matrix above. In this case,
%         for example, the elevation z values can be gotten by
%         terrainXYZs(xs, ys).
%   - region
%     The type of the terrain.
%       - For the NTIA C++ eHata libary, it is an integer for the NLCD
%         environment code.
%       - For the Matlab version of eHATA library by Thao Nguyen, it is a
%         string ('DenseUrban', 'Urban', or 'Suburban').
%   - eleProfileResForEHataInM
%     The elevation profile resolution in meter for the extended Hata
%     model. The smaller this value is, the more points between the
%     transmitter and the receiver will be considered to generate the
%     elevation profile for the path loss evaluation via the extended Hata
%     model.
%   - libraryToUse
%     A string specifying which eHATA model library to use:
%       - 'CPlusPlus'
%         The NITA C++ version of the eHATA library will be used. Source
%         code available at https://github.com/NTIA/ehata.
%       - 'Matlab'
%         The Matlab version of eHATA library by Thao Nguyen will be used.
%         Source code available at https://github.com/Thao-Nguyen/eHATA.
%   - NTIA_EHATA_RELIABILITY
%     Only required when the NTIA C++ eHata libary is used. It is a double
%     value from [0, 1] specifying the quantile percent not exceeded of the
%     signal.
%   - numOfWorkers
%     Additional. Default: 0. The number of parallel computing works to be
%     used (e.g. 0 for non parellel and inf for as many as system allowed).
%   - EHATA_VALID_DIST_RANGE_IN_M
%     Distance range to use eHata instead of the LoS FSPL model.
%
% Outputs:
%   - medianBPLMap
%     A matrix representing the median basic transmission loss map. The map
%     can be easily visualized via the command: image(medianBPLMap).
%   - medianBPLMapXLabels,  medianBPLMapYLabels
%     The x and y labels in the UTM system, respectively, for the output
%     map. Note that the Xs are increaseing and Ys are decreasing so that
%     the (i,j) element of medianBPLMap corresponds to the median basic
%     transmission loss at location [Xs(j), Ys(i)].
%
% Note that all the UTM locations should be in the same zone.
%
% Yaguang Zhang, Purdue, 06/12/2019

if ~exist('numOfWorkers', 'var')
    numOfWorkers = 0;
end

if ~exist('EHATA_VALID_DIST_RANGE_IN_M', 'var')
    % According to the Matlab eHata library.
    EHATA_VALID_DIST_RANGE_IN_M = [1000, 100000];
end

% According to Prof. Anderson.
EHATA_MIN_BASE_ANT_H_IN_M = 20;
% Any path loss above this will be considered as inf.
MAX_MEDIAN_BPL_ALLOWED_IN_DB = inf;

% Carrier wavelength.
lambdaInM = physconst('LightSpeed')/fsMHz/1e6;

medianBPLMapXLabels = sort(rxLocXs);
medianBPLMapYLabels = sort(rxLocYs, 'descend');

numRxXs = length(medianBPLMapXLabels);
numRxYs = length(medianBPLMapYLabels);

if isa(terrainXYZs, 'scatteredInterpolant') ...
        || isa(terrainXYZs, 'function_handle')
    fetchEles = terrainXYZs;
elseif isa(terrainXYZs, 'double')
    % This approach is slower.
    terrainXYZs = @(xs, ys) griddata( ...
        terrainXYZs(:,1), terrainXYZs(:,2), ...
        terrainXYZs(:,3), ...
        xs, ys);
else
    error(['Unsupported terrain information type: ', ...
        class(terrainXYZs), '!']);
end

if isa(lidarXYZs, 'scatteredInterpolant') ...
        || isa(lidarXYZs, 'function_handle')
    fetchLidarZs = lidarXYZs;
elseif isa(lidarXYZs, 'double')
    % This approach is slower.
    fetchLidarZs = @(xs, ys) griddata( ...
        lidarXYZs(:,1), lidarXYZs(:,2), ...
        lidarXYZs(:,3), ...
        xs, ys);
else
    error(['Unsupported lidar information type: ', ...
        class(lidarXYZs), '!']);
end

% For progress feedback.
ratioToReportProgress = 0.05;

medianBPLMap = nan(numRxYs, numRxXs);

if strcmpi(libraryToUse, 'cplusplus') ...
        && (EHATA_VALID_DIST_RANGE_IN_M(2)>=EHATA_VALID_DIST_RANGE_IN_M(1))
    ExtendedHata_PropLoss_CPP = @(elePro, fsMHz, baseAntHeightInM, ...
        mobileAntHeightInM, regionInt8, reliability) ...
        calllib('ehata', 'ExtendedHata', ...
        elePro, fsMHz, baseAntHeightInM, ...
        mobileAntHeightInM, regionInt8, reliability);
else
    ExtendedHata_PropLoss_CPP = nan;
end

disp('    ------');

% Pre-assign pixels to workers to avoid unnecessary data copying.
localCluster = gcp;
maxNumOfWorkers = min(localCluster.NumWorkers, numOfWorkers);

indicesRxXsYsForAllWorkers = cell(maxNumOfWorkers, 1);

pixelCnt = 0;
for idxRxX= 1:numRxXs
    for idxRxY = 1:numRxYs
        indicesRxXsYsForAllWorkers{mod(pixelCnt, maxNumOfWorkers)+1} = ...
            [indicesRxXsYsForAllWorkers{ ...
            mod(pixelCnt, maxNumOfWorkers)+1}; ...
            idxRxX idxRxY];
        pixelCnt = pixelCnt+1;
    end
end

allMedianBPLs = cell(maxNumOfWorkers,1);
parfor (idxWorker = 1:maxNumOfWorkers, maxNumOfWorkers)
    % Load the NTIA eHata library first, if necessary, to avoid the "unable
    % to find ehata" error.
    if strcmpi(libraryToUse, 'cplusplus') && (~libisloaded('ehata'))
        loadlibrary('ehata');
    end
    
    [curNumOfPixs,~] = size(indicesRxXsYsForAllWorkers{idxWorker});
    
    numPixsProcessed = 0;
    numPixsToReportProgress = ceil(curNumOfPixs.*ratioToReportProgress);
    
    curMedianBPLs = inf(curNumOfPixs, 1);
    for idxPixel = 1:curNumOfPixs
        if mod(numPixsProcessed, numPixsToReportProgress)==0
            disp(['    Worker #', num2str(idxWorker), ' (', ...
                num2str(numPixsProcessed/curNumOfPixs*100, ...
                '%.2f'), '%) ...']);
        end
        
        curIdxRxX = indicesRxXsYsForAllWorkers{idxWorker}(idxPixel, 1);
        curIdxRxY = indicesRxXsYsForAllWorkers{idxWorker}(idxPixel, 2);
        
        curRxX = medianBPLMapXLabels(curIdxRxX); %#ok<PFBNS>
        curRxY = medianBPLMapYLabels(curIdxRxY); %#ok<PFBNS>
        
        % Generate the elevation profile needed by the extended Hata model
        % according to our terrain information.
        %   - elev
        %     An array containing elevation profile between Tx & Rx, where:
        %       - elev(1) = numPoints - 1 (for Matlab eHata lib) or
        %       numPoints + 1 (for C++ eHata lib)
        %         (note, numPoints is the number of points between Tx & Rx)
        %       - elev(2) = distance between points (in meters).
        %         (thus, elev(1)-1)*elev(2)=distance between Tx & Rx)
        %       - elev(3) = Tx elevation
        %         (in meters)
        %       - elev(numPoints+2) = Rx elevation
        %         (in meters)
        distTxToRx = norm(baseAntXY-[curRxX curRxY]);
        
        % For generating the elevation/LiDAR z profile.
        numPoints = ceil(distTxToRx./eleProfileResForEHataInM);
        eleProfXs = linspace(baseAntXY(1), curRxX, numPoints)';
        eleProfYs = linspace(baseAntXY(2), curRxY, numPoints)';
        
        % For the extended Hata model, the distance between transmitter and
        % receiver should be in the valid distance range.
        if (distTxToRx>=EHATA_VALID_DIST_RANGE_IN_M(1)) ...
                &&(distTxToRx<=EHATA_VALID_DIST_RANGE_IN_M(2)) ...
                &&(baseAntHeightInM>=EHATA_MIN_BASE_ANT_H_IN_M) %#ok<PFBNS>
            curEleProfile = nan(numPoints+2, 1);
            curEleProfile(1) = numPoints-1;
            curEleProfile(2) = distTxToRx/(numPoints-1);
            curEleProfile(3:end) ...
                = fetchEles(eleProfXs, eleProfYs); %#ok<PFBNS>
            
            switch lower(libraryToUse)
                case 'cplusplus'
                    curMedianBPL = ExtendedHata_PropLoss_CPP( ...
                        curEleProfile, fsMHz, ...
                        baseAntHeightInM, mobileAntHeightInM, ...
                        int8(region), NTIA_EHATA_RELIABILITY);
                case 'matlab'
                    curMedianBPL = ExtendedHata_PropLoss( ...
                        fsMHz, baseAntHeightInM, mobileAntHeightInM, ...
                        region, curEleProfile);
                otherwise
                    error(['Unsupported library: ', model, '!'])
            end
            
            if curMedianBPL <= MAX_MEDIAN_BPL_ALLOWED_IN_DB
                curMedianBPLs(idxPixel) = curMedianBPL;
            end
        else
            % Use the free-space path loss model, instead.
            eleProfDists = linspace(0, distTxToRx, numPoints)';
            
            curLidarProfileZs ...
                = fetchLidarZs(eleProfXs, eleProfYs); %#ok<PFBNS>
            parsLoSPath = polyfit([0; distTxToRx], ...
                [baseAntHeightInM+fetchEles(baseAntXY(1),baseAntXY(2)); ...
                mobileAntHeightInM+fetchEles(curRxX, curRxY)], 1);
            curLosPathHs = polyval(parsLoSPath, eleProfDists);
            
            % Make sure there is no obstacles blocking the LoS path
            % according to the LiDAR data.
            if all(curLosPathHs>=curLidarProfileZs)
                curMedianBPLs(idxPixel) = fspl(distTxToRx, lambdaInM);
            end
        end
        numPixsProcessed = numPixsProcessed+1;
    end
    
    allMedianBPLs{idxWorker} = curMedianBPLs;
    
    disp(['    Worker #', num2str(idxWorker), ' (', ...
        num2str(100, '%.2f'), '%) ...']);
end

% Output the results.
for idxW = 1:maxNumOfWorkers
    % Note the order for indices (row, col) is different from that for UTM
    % labels (x, y).
    medianBPLMap(sub2ind(size(medianBPLMap), ...
        indicesRxXsYsForAllWorkers{idxW}(:, 2), ...
        indicesRxXsYsForAllWorkers{idxW}(:, 1))) ...
        = allMedianBPLs{idxW};
end

disp('    ------');

end
% EOF