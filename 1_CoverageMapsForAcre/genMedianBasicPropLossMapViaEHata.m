function [ MedianBPLMap, MedianBPLMapXLabels,  MedianBPLMapYLabels] ...
    = genMedianBasicPropLossMapViaEHata(fsMHz, ...
    baseAntXY, baseAntHeightInM, ...
    rxLocXs, rxLocYs, mobileAntHeightInM, ...
    terrainXYZs, region, ...
    eleProfileResForEHataInM, libraryToUse, NTIA_EHATA_RELIABILITY)
%GENMEDIANBASICPROPLOSSMAPVIAEHATA Compute the median basic propogation
%loss via the Extended Hata model for a grid to form a map.
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
%   - terrainXYZs
%     The terrain information for the area of interest. It can be:
%       - A matrix representing the terrain elevation information in the
%         form of [xs, ys, zs], where xs and ys are column vectors
%         reprenting the UTM locations where the elevation information is
%         available, and zs is a column vector storing the corresponding
%         elevations in meter. Or,
%       - A function to fetch the elevation for given locations, e.g. an
%         instance 'scatteredInterpolant' generated by scatteredInterpolant
%         from the terrain matrix above. In this case, the z values can be
%         gotten by terrainXYZs(xs, ys).
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
% Outputs:
%   - MedianBPLMap
%     A matrix representing the median basic transmission loss map. The map
%     can be easily visualized via the command: image(MedianBPLMap).
%   - MedianBPLMapXLabels,  MedianBPLMapYLabels
%     The x and y labels in the UTM system, respectively, for the output
%     map. Note that the Xs are increaseing and Ys are decreasing so that
%     the (i,j) element of MedianBPLMap corresponds to the median basic
%     transmission loss at location [Xs(j), Ys(i)].
%
% Note that all the UTM locations should be in the same zone.
%
% Yaguang Zhang, Purdue, 06/12/2019

MedianBPLMapXLabels = sort(rxLocXs);
MedianBPLMapYLabels = sort(rxLocYs, 'descend');

numRxXs = length(MedianBPLMapXLabels);
numRxYx = length(MedianBPLMapYLabels);

if isa(terrainXYZs, 'scatteredInterpolant') ...
        || isa(terrainXYZs, 'function_handle')
    fetchZs = terrainXYZs;
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

MedianBPLMap = inf(numRxYx, numRxXs);
for idxRxX= 1:numRxXs
    for idxRxY = 1:numRxYx
        curRxX = MedianBPLMapXLabels(idxRxX);
        curRxY = MedianBPLMapYLabels(idxRxY);
        
        % Generate the elevation profile needed by ExtendedHata_PropLoss
        % according to our terrain information.
        %   - elev
        %     An array containing elevation profile between Tx & Rx, where:
        %       - elev(1) = numPoints - 1
        %         (note, numPoints is the number of points between Tx & Rx)
        %       - elev(2) = distance between points (in meters).
        %         (thus, elev(1)-1)*elev(2)=distance between Tx & Rx)
        %       - elev(3) = Tx elevation
        %         (in meters)
        %       - elev(numPoints+2) = Rx elevation
        %         (in meters)
        distTxToRx = norm(baseAntXY-[curRxX curRxY]);
        
        % For the extended Hata model, the distance between transmitter and
        % receiver should be in the range of [1, 100] km
        if (distTxToRx>=1000) && (distTxToRx<=100000)
            numPoints = ceil(distTxToRx./eleProfileResForEHataInM);
            eleProfXs = linspace(baseAntXY(1), curRxX, numPoints)';
            eleProfYs = linspace(baseAntXY(2), curRxY, numPoints)';
            
            curEleProfile = nan(numPoints+2, 1);
            curEleProfile(1) = numPoints-1;
            curEleProfile(2) = distTxToRx/(numPoints-1);
            curEleProfile(3:end) = fetchZs(eleProfXs, eleProfYs);
            
            switch lower(libraryToUse)
                case 'cplusplus'
                    MedianBPLMap(idxRxY, idxRxX) ...
                        = calllib('ehata', 'ExtendedHata', ...
                        curEleProfile, fsMHz, ...
                        baseAntHeightInM, mobileAntHeightInM, ...
                        int8(region), NTIA_EHATA_RELIABILITY);
                case 'matlab'
                    MedianBPLMap(idxRxY, idxRxX) ...
                        = ExtendedHata_PropLoss( ...
                        fsMHz, baseAntHeightInM, mobileAntHeightInM, ...
                        region, curEleProfile);
                otherwise
                    error(['Unsupported library: ', model, '!'])
            end            
        end
    end
end

end
% EOF