% =========================================================================
%
% This software was developed by employees of the National Institute of 
% Standards and Technology (NIST), an agency of the Federal Government. 
% Pursuant to title 17 United States Code Section 105, works of NIST 
% employees are not subject to copyright protection in the United States 
% and are considered to be in the public domain. Permission to freely use, 
% copy, modify, and distribute this software and its documentation without 
% fee is hereby granted, provided that this notice and disclaimer of 
% warranty appears in all copies.
%
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, 
% EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, 
% ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY 
% IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, 
% AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION 
% WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE
% ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE FOR ANY DAMAGES, INCLUDING, 
% BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, 
% ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS 
% SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR 
% OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY 
% OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT 
% OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED 
% HEREUNDER.
%
% Distributions of NIST software should also include copyright and 
% licensing statements of any third-party software that are legally bundled
% with the code in compliance with the conditions of those licenses.
% 
% =========================================================================
%
% National Institute of Standards and Technology (NIST)
% Communications Technology Laboratory (CTL)
% Wireless Networks Division (673)
% By: Thao Nguyen 
% Date: 09/01/2016
%
% Description: Function to compute mixed land-sea path correction of 
% the extended Hata propagation model    
%
% Inputs: 
% - elev    : array containing elevation profile between Tx & Rx
%             where:
%             elev(1) = numPoints - 1 
%             (note, numPoints is the number of points between Tx & Rx)
%             elev(2) = distance between points (in meters). 
%             (thus, elev(1)-1)*elev(2)=distance between Tx & Rx)
%             elev(3) = Tx elevation (in meters)
%             elev(numPoints+2) = Rx elevation (in meters)
% Outputs:
% - Kmp     : correction factor for mixed path sea land (in decibels)
% 
% References: 
% [1] U.S. Department of Commerce, National Telecommunications and 
%     Information Administration, 3.5 GHz Exclusion Zone Analyses and 
%     Methodology (Jun. 18, 2015), available at 
%     http://www.its.bldrdoc.gov/publications/2805.aspx.
% [2] Y. Okumura, E. Ohmori, T. Kawano, and K. Fukuda, Field strength and
%     its variability in VHF and UHF land-mobile radio service, Rev. Elec. 
%     Commun. Lab., 16, 9-10, pp. 825-873, (Sept.-Oct. 1968).
% 
% History:
% - 2016/09/01: released version 1.0


function Kmp = ExtendedHata_MixedPathCorr(elev)

% Extract data from elevation profile 
numPoints = elev(1) + 1;             % number of points between Tx & Rx
pointRes_km = elev(2)/1e3;           % distance between points (km)
pointElev_m = elev(3:2+numPoints);   % elevation vector (m)
d_Tx_Rx_km = (numPoints-1)*pointRes_km;% distance between Tx & Rx (km)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute sea-land distance ratio 'beta' and 
% determine whether:
% - mobile station is adjacent to sea (case A) or
% - base station is adjacent to sea (case B) or
% - sea is in the middle between base station and mobile station (case C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find points that have elevation = 0 (i.e., sea level)
indSeaPoint = find(pointElev_m == 0);

% Return if path has no sea point
if isempty(indSeaPoint)
    Kmp = 0;
    return;
end

% Sea distance
ds_km = length(indSeaPoint)*pointRes_km;

% Distance ratio between sea distance and total distance (in percent)
beta = ds_km/d_Tx_Rx_km*100;

% Determine whether the sea is adjacent to the Tx, Rx, or in between
boundary_km = 1;
boudary_point = round(boundary_km / pointRes_km);

adjTx = find(indSeaPoint<=boudary_point); % find sea points near Tx
adjRx = find(indSeaPoint>=numPoints-boudary_point);%find sea points near Rx
if (isempty(adjTx) && ~isempty(adjRx))    % Only adjacent to the Rx
    curve = 'A';
elseif (~isempty(adjTx) && isempty(adjRx))% Only adjacent to the Tx
    curve = 'B';
else
    curve = 'C';                     % Sea can be anywhere between Tx & Rx 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute correction factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the curves from Figure 35 in Okumura's paper
betaVec = 0:10:100;   % distance ratio (percentage)
Kmp_dGreater60km_A = [0 3 5.3 7.5 9.2 10.5 11.8 12.7 13.5 14.2 14.8];
Kmp_dGreater60km_B = [0 1.5 2.5 4 5.5 7 9 10.8 12.5 14.2 14.8];

Kmp_dLessThan30km_A = [0 2 3.8 5 6.3 7.5 8.5 9.3 10 10.5 11];
Kmp_dLessThan30km_B = [0 0.5 1.5 2.5 3.5 4.8 6 7.8 9.5 10.5 11];

% Apply piecewise linear fits to find the mixed land-sea path
% correction factor
if (d_Tx_Rx_km > 60)           % Distance > 60 km
    
    if strcmp(curve, 'A')      % mobile station is adjacent to sea (case A)
       
        Kmp = piecelin(betaVec, Kmp_dGreater60km_A, beta);
        
    elseif strcmp(curve, 'B')   % base station is adjacent to sea (case B)
        
        Kmp = piecelin(betaVec, Kmp_dGreater60km_B, beta);
        
    elseif strcmp(curve, 'C')   % sea is in the middle between base station 
                                % and mobile station (case C) 
                                % => average of curve A and curve B
        Kmp = (piecelin(betaVec, Kmp_dGreater60km_A, beta) + ...
            piecelin(betaVec, Kmp_dGreater60km_B, beta)) / 2; 
        
    end
    
elseif (d_Tx_Rx_km < 30)        % Distance < 30 km
    
    if strcmp(curve, 'A')      % mobile station is adjacent to sea (case A)
       
        Kmp = piecelin(betaVec, Kmp_dLessThan30km_A, beta);
        
    elseif strcmp(curve, 'B')   % base station is adjacent to sea (case B)
        
        Kmp = piecelin(betaVec, Kmp_dLessThan30km_B, beta);
        
    elseif strcmp(curve, 'C')   % sea is in the middle between base station 
                                % and mobile station (case C) 
                                % => average of curve A and curve B
        Kmp = (piecelin(betaVec, Kmp_dLessThan30km_A, beta) + ...
            piecelin(betaVec, Kmp_dLessThan30km_B, beta)) / 2; 
        
    end
    
else                            % 30 km <= Distance <= 60 km
    
     if strcmp(curve, 'A')      % mobile station is adjacent to sea (case A)
       
         % Get correction factors for distance >60km and <30km
         KmpGT60km = piecelin(betaVec, Kmp_dGreater60km_A, beta); 
         KmpLT30km = piecelin(betaVec, Kmp_dLessThan30km_A, beta);
        
         % Interpolate correction factors
         Kmp = interp1([30, 60], [KmpLT30km, KmpGT60km], d_Tx_Rx_km);
        
    elseif strcmp(curve, 'B')   % base station is adjacent to sea (case B)
        
         % Get correction factors for distance >60km and <30km
         KmpGT60km = piecelin(betaVec, Kmp_dGreater60km_B, beta); 
         KmpLT30km = piecelin(betaVec, Kmp_dLessThan30km_B, beta);
        
         % Interpolate correction factors
         Kmp = interp1([30, 60], [KmpLT30km, KmpGT60km], d_Tx_Rx_km);
                
    elseif strcmp(curve, 'C')   % sea is in the middle between base station 
                                % and mobile station (case C) 
                                % => interpolate based on the distance
        
         % Get correction factors for distance >60km and <30km
         KmpGT60km = (piecelin(betaVec, Kmp_dGreater60km_A, beta) + ...
             piecelin(betaVec, Kmp_dGreater60km_B, beta)) / 2; 
         KmpLT30km = (piecelin(betaVec, Kmp_dLessThan30km_A, beta) + ...
             piecelin(betaVec, Kmp_dLessThan30km_B, beta)) / 2;
        
         % Interpolate correction factors
         Kmp = interp1([30, 60], [KmpLT30km, KmpGT60km], d_Tx_Rx_km);
        
    end
end



