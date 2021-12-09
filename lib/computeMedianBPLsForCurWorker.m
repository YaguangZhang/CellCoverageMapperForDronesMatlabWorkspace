function curMedianBPLs = computeMedianBPLsForCurWorker(idxWorker, ...
    indicesRxXsYsForAllWorkers, ...
    medianBPLMapXLabels, medianBPLMapYLabels, baseAntXY, ...
    eleProfileResForEHataInM, MIN_NUM_ELE_SAMPLES_PER_PROFILE, ...
    EHATA_VALID_DIST_RANGE_IN_M, baseAntHeightInM, ...
    mobileAntHeightInM, libraryToUse, ...
    ExtendedHata_PropLoss_CPP, ...
    fetchEles, fetchLidarZs, fsMHz, region, ...
    NTIA_EHATA_RELIABILITY, EHATA_MIN_BASE_ANT_H_IN_M, ...
    MAX_MEDIAN_BPL_ALLOWED_IN_DB, lambdaInM)
%COMPUTEMEDIANBPLSFORCURWORKER Compute the median basic path loss for one
%worker in the parfor loop.
%
% This is just a helper function encapsulating the needed script to
% facilitate debugging.
%
% Yaguang Zhang, Purdue, 08/06/2019

% For progress feedback.
ratioToReportProgress = 0.05; % One feedback per 5% progress.

% The clearance ratio of the first Fresnel zone for a LoS path: at least
% this ratio of the first Fresnel zone needs to be clear for a path to be
% considered as "line of sight" (LoS); we expect it to be larger or equal
% to 50%.
LOS_FIRST_FRES_CLEAR_RATIO = 0.6;

[curNumOfPixs,~] = size(indicesRxXsYsForAllWorkers{idxWorker});

numPixsProcessed = 0;
numPixsToReportProgress ...
    = ceil(curNumOfPixs.*ratioToReportProgress);

curMedianBPLs = inf(curNumOfPixs, 1);

for idxPixel = 1:curNumOfPixs
    if mod(numPixsProcessed, numPixsToReportProgress)==0
        disp(['    Worker #', num2str(idxWorker), ' (', ...
            num2str(numPixsProcessed/curNumOfPixs*100, ...
            '%.2f'), '%) ...']);
    end
    try
        curIdxRxX = indicesRxXsYsForAllWorkers{idxWorker}(idxPixel, 1);
        curIdxRxY = indicesRxXsYsForAllWorkers{idxWorker}(idxPixel, 2);
        
        curRxX = medianBPLMapXLabels(curIdxRxX);
        curRxY = medianBPLMapYLabels(curIdxRxY);
        
        % Generate the elevation profile needed by the extended Hata model
        % according to our terrain information.
        %   - elev
        %     An array containing elevation profile between Tx & Rx, where:
        %       - elev(1) = numPoints - 1
        %         (for both Matlab eHata lib and C++ eHata lib; note,
        %         numPoints is the number of points between Tx & Rx)
        %       - elev(2) = distance between points (in meters).
        %         (thus, elev(1)*elev(2)=distance between Tx & Rx)
        %       - elev(3) = Tx elevation
        %         (in meters)
        %       - elev(numPoints+2) = Rx elevation
        %         (in meters)
        distTxToRx = norm(baseAntXY-[curRxX curRxY]);
        
        % For generating the elevation/LiDAR z profile.
        numPoints = max( ...
            ceil(distTxToRx./eleProfileResForEHataInM), ...
            MIN_NUM_ELE_SAMPLES_PER_PROFILE);
        
        % Only evaluate the path loss if the distance is not too big.
        if distTxToRx<=EHATA_VALID_DIST_RANGE_IN_M(2)
            % Here we need to construct the following parameters for the
            % eHata model:
            % 	- txHeightInM, rxHeightInM
            %     Elevation + antenna height for the TX and RX,
            %     respectively.
            %   - curEleProfile
            %     The terrain profile between the TX and the RX.
            eleProfXs = linspace(baseAntXY(1), curRxX, numPoints)';
            eleProfYs = linspace(baseAntXY(2), curRxY, numPoints)';
            
            % Path loss evaluation via eHata.
            curEleProfile = nan(numPoints+2, 1);
            curEleProfile(1) = numPoints-1;
            curEleProfile(2) = distTxToRx/(numPoints-1);
            
            baseAntElePlusHeight ...
                = baseAntHeightInM ...
                + fetchEles(baseAntXY(1),baseAntXY(2));
            mobileAntElePlusHeight ...
                = mobileAntHeightInM + fetchEles(curRxX, curRxY);
            % Inputs for the eHata model.
            if baseAntElePlusHeight>=mobileAntElePlusHeight
                txHeightInM = baseAntHeightInM;
                rxHeightInM = mobileAntHeightInM;
                curEleProfile(3:end) ...
                    = fetchEles(eleProfXs, eleProfYs);
                effectiveBaseAntHInM = baseAntHeightInM;
            else
                % Assuming reciprocal channels, we will flip TX and RX for
                % cases where the mobile antenna is higher than the
                % cellular tower antenna, because eHata is designed for
                % higher TXs.
                txHeightInM = mobileAntHeightInM;
                rxHeightInM = baseAntHeightInM;
                curEleProfile(3:end) ...
                    = fetchEles(eleProfXs(end:-1:1), ...
                    eleProfYs(end:-1:1));
                effectiveBaseAntHInM = mobileAntHeightInM;
            end

            % We will consider the 3D distance in FSPL computation.
            distTxToRx3D = norm([distTxToRx, txHeightInM-rxHeightInM]);
            curFsplMedianBPL = fspl(distTxToRx3D, lambdaInM);

            switch lower(libraryToUse)
                case 'cplusplus'    % C++ eHata
                    curEHataMedianBPL = ExtendedHata_PropLoss_CPP( ...
                        curEleProfile, fsMHz, ...
                        txHeightInM, rxHeightInM, ...
                        int8(region), NTIA_EHATA_RELIABILITY);
                case 'matlab'       % Matlab eHata
                    curEHataMedianBPL = ExtendedHata_PropLoss( ...
                        fsMHz, txHeightInM, rxHeightInM, ...
                        region, curEleProfile);
                case 'fspl'         % LoS FSPL
                    % Make sure there is no obstacles blocking the LoS path
                    % too much according to the LiDAR data. Note that here
                    % we always conider the profiles with the cellular
                    % tower as the start point and the mobile device as the
                    % end point.
                    
                    % We only need to compare the points between the
                    % cellular tower and the mobile device.
                    curLidarProfileZs ...
                        = fetchLidarZs(eleProfXs(2:(end-1)), ...
                        eleProfYs(2:(end-1)));
                    
                    lidarProfDists = linspace(0, distTxToRx, numPoints)';
                    lidarProfDists = lidarProfDists(2:(end-1));
                    
                    parsLoSPath = polyfit([0; distTxToRx], ...
                        [baseAntElePlusHeight; mobileAntElePlusHeight], 1);
                    curLosPathHs = polyval(parsLoSPath, lidarProfDists);
                    
                    if all(curLosPathHs>=curLidarProfileZs)
                        % A 50-percent clearance is now ensured. We need to
                        % consider the first Fresnel zone here to more
                        % accurately determine the LoS paths.
                        
                        % Distance between the celluar tower and the lidar
                        % z locations.
                        d1s = vecnorm( ...
                            [lidarProfDists(:)'; ...
                            curLidarProfileZs(:)'-baseAntElePlusHeight]);
                        % Distance between the lidar z locations and the
                        % mobile device.
                        d2s = vecnorm( ...
                            [distTxToRx - lidarProfDists(:)'; ...
                            curLidarProfileZs(:)'-mobileAntElePlusHeight]);
                        % First Fresnel zone radii for the lidar z
                        % locations.
                        firstFresRadii ...
                            = (sqrt( ...
                            (lambdaInM .* d1s .* d2s)./(d1s + d2s)))';
                        
                        % The distances between the lidar z locations and
                        % the TX-RX direct path.
                        distsToDirectPath ...
                            = point_to_line_distance( ...
                            [lidarProfDists(:), curLidarProfileZs(:)], ...
                            [0, baseAntElePlusHeight], ...
                            [distTxToRx, mobileAntElePlusHeight]);
                        
                        % The extra clearance ratio respective to the first
                        % Fresnel zone radius for LoS paths.
                        extraClearRatioVsR ...
                            = (LOS_FIRST_FRES_CLEAR_RATIO-0.5)*2;
                        
                        if all(distsToDirectPath ...
                                -extraClearRatioVsR*firstFresRadii>0)
                            curMedianBPLs(idxPixel) = curFsplMedianBPL;
                        end
                    end
                otherwise
                    error(['Unsupported library: ', libraryToUse, '!'])
            end
            
            if strcmpi(libraryToUse, 'cplusplus') ...
                    || strcmpi(libraryToUse, 'matlab')
                % For the extended Hata model, the distance between
                % transmitter and receiver should be in the valid distance
                % range.
                if (distTxToRx>=EHATA_VALID_DIST_RANGE_IN_M(1)) ...
                        &&(effectiveBaseAntHInM>=EHATA_MIN_BASE_ANT_H_IN_M)
                    % No need to consider the FSPL model.
                    if curEHataMedianBPL <= MAX_MEDIAN_BPL_ALLOWED_IN_DB
                        curMedianBPLs(idxPixel) = curEHataMedianBPL;
                    end
                else
                    % Consider the FSPL model.
                    ratioForEHata ...
                        = distTxToRx./EHATA_VALID_DIST_RANGE_IN_M(1);
                    curMedianBPLs(idxPixel) ...
                        = ratioForEHata.*curEHataMedianBPL ...
                        +(1-ratioForEHata).*curFsplMedianBPL;
                end
            end
        end
        
        numPixsProcessed = numPixsProcessed+1;
    catch err
        disp(['Error occurred in parfor idxWorker = ', ...
            num2str(idxWorker), ' idxPixel = ', num2str(idxPixel)]);
        rethrow(err);
    end
end

disp(['    Worker #', num2str(idxWorker), ' (', ...
    num2str(100, '%.2f'), '%) ...']);

end
% EOF