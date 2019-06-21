September 01, 2016

The Extended Hata (eHATA) propagation model code was implemented by employees 
of the National Institute of Standards and Technology (NIST), Communications 
Technology Laboratory (CTL). The code was written in MATLAB and the latest 
version of the code can be found at:

https://github.com/Thao-Nguyen

Version 1.0 release contains the following files and folders:

********************************** Files ***************************************

README                 				- This file.

ExtendedHata_PropLoss.m 			- The main function to compute the extended 
									  Hata propagation loss.

ExtendedHata_MedianBasicPropLoss.m 	- A function to compute median basic 
								      transmission loss.

ExtendedHata_EffHeightCorr.m 		- A function to compute terminals's 
									  "effective height" corrections.

ExtendedHata_RollingHillyCorr.m     - A function to compute "median" and 
									  "fine" corrections for rolling hilly terrain.

ExtendedHata_GeneralSlopeCorr.m 	- A function to compute general slope 
									  of terrain correction.

ExtendedHata_IsolatedRidgeCorr.m    - A function to compute isolated 
									  mountain (or isolated ridge) correction.

ExtendedHata_MixedPathCorr.m  		- A function to compute mixed 
									  land-sea path correction.

ExtendedHata_LocationVariability.m  - A function to compute the standard 
									  deviation estimate of the location  
									  variability for urban and suburban 
									  environments. 

piecelin.m 							- A function to get piecewise linear 
								      interpolation. It can be downloaded at:
									  https://www.mathworks.com							  

doc									- A folder contains a short tutorial  
									  explaining how to use the code.
							
tests								- A folder contains MATLAB scripts and   
									  a .mat file used to test each function. 
									  
									  Accuracy of the ExtendHata_MedianBasicPropLoss.m 
									  function, which is deterministic and does  
									  not depend on terrain information, has been  
									  validated by comparing with Paul McKenna's 
									  results.
										
									  However, other functions, which are not 
									  deterministic and depend on terrain  
									  information, has been tested using only a 
									  subset of terrain data. These functions 
									  might need more testing against reference 
									  data.

********************************************************************************

********************************************************************************

The folder containing the eHATA code needs to be added to the MATLAB path 
before use.

This can be accomplished by running the 'addpath' command from the 
matlab command prompt. For example:

>> addpath('C:\MATLAB\eHATA')

********************************************************************************

********************************************************************************

Legal Notice:
This software was developed by employees of the National Institute of Standards 
and Technology (NIST), an agency of the Federal Government. Pursuant to 
title 17 United States Code Section 105, works of NIST employees are not 
subject to copyright protection in the United States and are considered to 
be in the public domain. Permission to freely use, copy, modify, and distribute 
this software and its documentation without fee is hereby granted, provided that 
this notice and disclaimer of warranty appears in all copies.

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER 
EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY 
THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, 
AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY 
WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NASA BE LIABLE 
FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR 
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED 
WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR 
OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR 
OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE 
RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.

Distributions of NIST software should also include copyright and 
licensing statements of any third-party software that are legally bundled 
with the code in compliance with the conditions of those licenses. 

********************************************************************************

********************************************************************************

References: 
[1] U.S. Department of Commerce, National Telecommunications and 
    Information Administration, 3.5 GHz Exclusion Zone Analyses and 
    Methodology (Jun. 18, 2015), available at 
    http://www.its.bldrdoc.gov/publications/2805.aspx.
[2] Y. Okumura, E. Ohmori, T. Kawano, and K. Fukuda, Field strength and
    its variability in VHF and UHF land-mobile radio service, Rev. Elec. 
    Commun. Lab., 16, 9-10, pp. 825-873, (Sept.-Oct. 1968).
[3] M. Hata, Empirical formula for propagation loss in land mobile radio
    services, IEEE Transactions on Vehicular Technology, VT-29, 3,
    pp. 317-325 (Aug. 1980).
[4] Anita G. Longley, Radio Propagation in Urban Areas, United States 
    Department of Commerce, Office of Telecommunications, OT Report 
    78-144 (Apr.1978), available at 
    http://www.its.bldrdoc.gov/publications/2674.aspx.
	
********************************************************************************


