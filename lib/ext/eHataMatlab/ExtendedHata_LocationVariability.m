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
% Description: Function to compute the estimate of the standard deviation
% of the location variability for urban and suburban environments of the
% extended Hata propagation model.
%
% Inputs:
% - freq_MHz    : frequency (in MHz), in the range of [1500, 3000] MHz
% - region      : region of the area ('DenseUrban', 'Urban', 'Suburban')
%
% Outputs:
% - sigma       : standard deviation of the location variability 
%                 (in decibels)
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

function sigma = ExtendedHata_LocationVariability(freq_MHz, region)

% Compute standard deviation of the location variability for urban region
% ([1] Eqn. (A-18c))
sigma = 4.0976291 - 1.2255656*log10(freq_MHz) + 0.68350345*(log10(freq_MHz))^2;

% Adjust std for suburban region ([1] Page 35)
if strcmp(region, 'Suburban')
    sigma = sigma + 2;
end;





