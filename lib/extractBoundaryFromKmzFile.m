function [ utmXyBoundary, utmZone, kmzStruct, ...
    utmXYPolygons, lonLatPolygons] ...
    = extractBoundaryFromKmzFile(dirToKmzFile, flagSkipValidityCheck)
%EXTRACTBOUNDARYFROMKMZFILE Extract the boundary of the union of all the
%polygons stored in the input .kmz file.
%
% Input:
%   - dirToKmzFile
%     The absolute path to the .kmz file. Note that we only support
%     polygons without holes.
%
% Outputs:
%   - utmXyBoundary
%     The output polygon matrix for the output boundary, with each row
%     being one vertex in the form of UTM (x, y).
%   - utmZone
%     The UTM zone label for the coordinates, e.g., '16 T'.
%   - kmzStruct
%     A struct array for the raw items extracted from the .kmz file.
%   - utmXYPolygons, lonLatPolygons
%     Column polyshape arrays for "Polygon" items in the .kmz file, in the
%     UTM (x, y) and GPS (lon, lat) systems, respectively. Note that we
%     flipped the (lat, lon) to make it easier to plot.
%
% Yaguang Zhang, Purdue, 06/03/2020

if ~exist('flagSkipValidityCheck', 'var')
    flagSkipValidityCheck = false;
end

utmXyPoly = polyshape();
utmZone = '';

kmzStruct = kmz2structCxPlatform(dirToKmzFile);
polyonsLonLatStruct = kmzStruct(strcmp({kmzStruct.Geometry}, 'Polygon'));

numOfPolygons = length(polyonsLonLatStruct);
if nargout > 3
    [utmXYPolygons, lonLatPolygons] = deal( ...
        repmat(polyshape, numOfPolygons, 1));
end

for idxPolygon = 1:numOfPolygons
    curPolygonStruct = polyonsLonLatStruct(idxPolygon);
    
    % Make sure we have column vectors to work with.
    curLats = curPolygonStruct.Lat(:);
    % Remove the trailing NaN.
    curLats = curLats(1:(end-1));
    
    if ~isempty(curLats)
        curLons = curPolygonStruct.Lon(:);
        curLons = curLons(1:(end-1));
        
        assert(~any(isnan([curLats; curLons])), ...
            'We expect polygons without holes!');
        
        [curXs, curYs, curZones] = deg2utm(curLats, curLons);
        if isempty(utmZone)
            utmZone = curZones(1,:);
            assert(all(strcmp(cellstr(curZones),utmZone)), ...
                'Not all polygons are contained in the same UTM zone!');
        end
        
        curUtmXYPolygon = polyshape(curXs, curYs);
        if ~isempty(utmXyPoly.Vertices)
            utmXyPoly = union(utmXyPoly, curUtmXYPolygon);
        else
            utmXyPoly = curUtmXYPolygon;
        end
        
        if nargout > 3
            utmXYPolygons(idxPolygon) = curUtmXYPolygon;
            lonLatPolygons(idxPolygon) = polyshape(curLons, curLats);
        end
    end
end

utmXyBoundary = utmXyPoly.Vertices;

if ~flagSkipValidityCheck
    % Make sure the output vertices are clockwise and the polygon is
    % closed.
    if ~ispolycw(utmXyBoundary(:,1),utmXyBoundary(:,2))
        [utmXyBoundary(:,1),utmXyBoundary(:,2)] = poly2cw( ...
            utmXyBoundary(:,1),utmXyBoundary(:,2));
    end
    
    assert(all(ispolycw(utmXyBoundary(:,1),utmXyBoundary(:,2))), ...
        'The output vertices are not clockwise!')
end
utmXyBoundary(end+1, :) = utmXyBoundary(1, :);

end
% EOF