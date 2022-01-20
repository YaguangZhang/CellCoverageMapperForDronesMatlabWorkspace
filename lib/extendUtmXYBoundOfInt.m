function [extendedUtmXYBound] ...
    = extendUtmXYBoundOfInt(utmXYBoundOfInt, dInM)
%EXTENDUTMXYBOUNDOFINT Extend a polygon, defined by the boundary vertices
%in UTM (x, y), by distance dInM.
%
% Yaguang Zhang, Purdue, 01/20/2022

extendedUtmXYBoundPoly ...
    = polybuffer(polyshape(utmXYBoundOfInt), dInM);
assert(extendedUtmXYBoundPoly.NumRegions==1, ...
    'The extended area of interest has more than one region!');

extendedUtmXYBound = extendedUtmXYBoundPoly.Vertices;

end
% EOF