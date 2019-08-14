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
ratioToReportProgress = 0.05;

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
            
            curFsplMedianBPL = fspl(distTxToRx, lambdaInM);
            
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
                    eleProfDists = linspace(0, distTxToRx, numPoints)';
                    
                    curLidarProfileZs ...
                        = fetchLidarZs(eleProfXs, eleProfYs);
                    
                    parsLoSPath = polyfit([0; distTxToRx], ...
                        [baseAntElePlusHeight; mobileAntElePlusHeight], 1);
                    curLosPathHs = polyval(parsLoSPath, eleProfDists);
                    
                    % Make sure there is no obstacles blocking the LoS path
                    % according to the LiDAR data.
                    if all(curLosPathHs>=curLidarProfileZs)
                        curMedianBPLs(idxPixel) = curFsplMedianBPL;
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
                    % Consider the LoS FSPL model and weigh it in when its
                    % evaluation is valid.                    
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