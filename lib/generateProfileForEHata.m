function [elev] = generateProfileForEHata(txXYH, rxXYH, terrainProfile)
%GENERATEPROFILEFOREHATA Encapsulate terrain information to a profile for
%eHata.
% Inputs:
%   - txXYH, rxXYH
%     The (x, y, heightInM) coordinates for the Tx and the Rx,
%     respectively.
%   - terrainProfile
%     The terrain z values for points between Tx and Rx (including Tx & Rx,
%     i.e., with them as end points).
%
% Output:
%   - elev
%     An array containing elevation profile between Tx & Rx (including Tx &
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
%
% Yaguang Zhang, Purdue, 02/11/2022

numTerrianPoints = length(terrainProfile);

elev = nan(numTerrianPoints+2, 1);
elev(1) = numTerrianPoints-1;
distTxToRx = norm(txXYH(1:2)-rxXYH(1:2));
elev(2) = distTxToRx/(numTerrianPoints-1);
elev(3:end) = terrainProfile;

end
% EOF