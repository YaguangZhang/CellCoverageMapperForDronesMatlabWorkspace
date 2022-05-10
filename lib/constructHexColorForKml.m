function [colorHex] = constructHexColorForKml(color, alpha)
%CONSTRUCTHEXCOLORFORKML Construct a hex color variable for plotting in a
%.kml file. All input values, including elements of [r, g, b] color and
%alpha, are 0~255 integers.
%
% Yaguang Zhang, Purdue, 02/10/2021

if ~exist('alpha', 'var')
    alpha = 255;
end

[r, g, b] = deal(color(1), color(2), color(3));

[rhex, ghex, bhex, ahex ] = deal( ...
    dec2hex(r),dec2hex(g),dec2hex(b),dec2hex(alpha));
if length(rhex)==1,rhex=['0' rhex];end
if length(ghex)==1,ghex=['0' ghex];end
if length(bhex)==1,bhex=['0' bhex];end
if length(ahex)==1,ahex=['0' ahex];end

colorHex = [ahex rhex ghex bhex];
end
% EOF