%DEFINEMODELITUOBSBYWOODLAND
%
% Yaguang Zhang, Purdue, 05/07/2021

% ITU-R P.833-9 obstruction by woodland model. Note that we do not
% differentiate the signal polarization here.
modelItuObsByWoodland.refParameterValues.freqInMHz ...
    = [105.9, 466.475, 949.0, 1852.2, 2117.5];
% Specific attenuation, gamma, in dB/m.
modelItuObsByWoodland.refParameterValues.gammaInDbPerM ...
    = [0.04, 0.12, 0.17, 0.30, 0.34];
% The maximum attenuation, Am, in dB.
modelItuObsByWoodland.refParameterValues.AmInDb ...
    = [9.4, 18.0, 26.5, 29.0, 34.1];

% Key points in figure 2 for estimating gamma. Each row is a point [fInMHz,
% gamma]. We will use the vertical polaization for the worst-case analyses.
modelItuObsByWoodland.refFittedGammaWrtFInMHz ...
    = [30, 2*10^(-2); 1000, 2*10^(-1); 30000, 6];

% The fomular to compute the excess loss caused by vegetation, where d is
% the length of path within woodland in meter.
modelItuObsByWoodland.excessLossFormula ...
    = @(d, gamma, Am) Am.*(1-exp(-d.*gamma./Am));

% Linearly interpolate and extrapolate the reference parameter values in
% log scale relative to the input frequency fInMHz based on Rec. ITU-R
% P.833-7 FIGURE 2 Specific attenuation due to woodland.
modelItuObsByWoodland.estimateGamma = @(fInMHz) 10^interp1( ...
    log10(modelItuObsByWoodland.refFittedGammaWrtFInMHz(:,1)), ...
    log10(modelItuObsByWoodland.refFittedGammaWrtFInMHz(:,2)), ...
    log10(fInMHz), 'linear','extrap'); 
% Rec. ITU-R P.833-7 Formula 2 with parameters from the Mulhouse (France)
% measurement campaign.
modelItuObsByWoodland.estimateAm = @(fInMHz) 1.15 * fInMHz^0.43;
% modelItuObsByWoodland.estimateAm = @(fInMHz) 10^interp1( ...
%     log10(modelItuObsByWoodland.refParameterValues.freqInMHz), ...
%      log10(modelItuObsByWoodland.refParameterValues.AmInDb), ...
%     log10(fInMHz), 'linear','extrap');

% Tests.
%  modelItuObsByWoodland.estimateGamma(28000)
% modelItuObsByWoodland.estimateAm(28000)

% EOF