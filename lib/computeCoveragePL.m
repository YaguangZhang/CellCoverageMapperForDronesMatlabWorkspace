function coveragePL = computeCoveragePL(txXYH, rxXYH, ...
    elePro, simConfigs, lidarProfile)
%COMPUTECOVERAGEPL Compute the path loss (in dB) for a pixel of the
%coverage map.
%
% For the coverage map, we will
%   (1) For drone locations with distances over (and including) 1km from
%   the cell tower
%       We use NTIA eHata model to evaluate the path loss.
%   (2) For drone locations with distances from the cell tower strictly
%   under 1km
%       If the path loss from the blockage map is valid for this drone
%       location, we will weight that value with the eHata value (please
%       see the code below for more details); otherwise, NaN will be
%       assigned for the path loss of that location.
%
% Inputs:
%   - txXYH, rxXYH
%     The UTM (x,y) with height (elvation = elevation + height) for the TX
%     and the RX, respectively.
%   - elePro
%     The terrain pofile defined by the eHata model. It is essentially an
%     array containing elevation profile between Tx & Rx (including Tx &
%     Rx, i.e., with them as end points), where:
%       - elev(1) = numPoints - 1
%         (for both Matlab eHata lib and C++ eHata lib; note, numPoints is
%         the number of points between Tx & Rx)
%       - elev(2) = distance between points (in meters).
%         (thus, elev(1)*elev(2)=distance between Tx & Rx)
%       - elev(3) = Tx elevation
%         (in meters)
%       - elev(numPoints+2) = Rx elevation
%         (in meters)
%   - simConfigs
%     A structure holding basic parameters of the simulation; we need here:
%       - simConfigs.NTIA_EHATA_RELIABILITY
%         For the NTIA eHata library: the quantile percent not exceeded of
%         the signal. Limits: 0 < reliability < 1.
%       - simConfigs.NTIA_EHATA_ENVIRO_CODE
%   	  For the NTIA eHata library: the NLCD environment code.
%       - simConfigs.CARRIER_WAVELENGTH_IN_M
%         Carrier wavelength in meters.
%   - lidarProfile
%     The LiDAR pofile (LiDAR z in meters) between (and including) the TX
%     and the RX. Only needed when distTxToRx<1km.
%
% Output:
%   - blockagePL
%     The path loss for the blockage map. If the link is blocked, a NaN
%     will be returned.
%
% Yaguang Zhang, Purdue, 09/18/2019

% Set this to true to ignore obstacles for FSPL computation within the
% distance threshold MIN_TX_TO_RX_DIST_FOR_EHATA_IN_M. Otherwise, we will
% output NaN if there is obstacles blocking the 1st Fresnel zone.
FLAG_IGNORE_OBSTACLES_FOR_FSPL = true;

% The minimum and maximum TX-to-RX distance for which the eHata model is
% valid.
MIN_TX_TO_RX_DIST_FOR_EHATA_IN_M = 1000;
MAX_TX_TO_RX_DIST_FOR_EHATA_IN_M = 100000;
regionInt8 = int8(simConfigs.NTIA_EHATA_ENVIRO_CODE);

if ~libisloaded('ehata')
    loadlibrary('ehata');
end

distTxToRx = norm(txXYH(1:2)-rxXYH(1:2));
coveragePL = nan;

% eHata is designed for the case where TX is higher than RX.
txAlt = elePro(3)+txXYH(3);
rxAlt = elePro(end)+rxXYH(3);
assert(txAlt>=rxAlt, 'TX is lower than RX!');

if distTxToRx<=MAX_TX_TO_RX_DIST_FOR_EHATA_IN_M
    curEHataPL = calllib('ehata', 'ExtendedHata', ...
        elePro, simConfigs.CARRIER_FREQUENCY_IN_MHZ, txXYH(3), ...
        rxXYH(3), regionInt8, simConfigs.NTIA_EHATA_RELIABILITY);
end

if distTxToRx<MIN_TX_TO_RX_DIST_FOR_EHATA_IN_M
    if FLAG_IGNORE_OBSTACLES_FOR_FSPL
        % We will consider the 3D distance in FSPL computation.
        distTxToRx3D = norm([txXYH(1:2), txAlt] - [rxXYH(1:2), rxAlt]);
        curFSPL = fspl(distTxToRx3D, simConfigs.CARRIER_WAVELENGTH_IN_M);
    else
        % Consider the non-blocked FSPL model.
        curFSPL = computeBlockagePL([txXYH(1:2), txAlt], ...
            [rxXYH(1:2), rxAlt], ...
            lidarProfile, simConfigs); %#ok<UNRCH>
    end

    ratioForEHata ...
        = distTxToRx./MIN_TX_TO_RX_DIST_FOR_EHATA_IN_M;
    coveragePL ...
        = ratioForEHata.*curEHataPL ...
        +(1-ratioForEHata).*curFSPL;
elseif (distTxToRx>=MIN_TX_TO_RX_DIST_FOR_EHATA_IN_M) ...
        && (distTxToRx<=MAX_TX_TO_RX_DIST_FOR_EHATA_IN_M)
    coveragePL = curEHataPL;
end

end

% EOF