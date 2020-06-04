function [ utmXyBoundary, utmZone ] ...
    = extractBoundaryFromKmzFile(dirToKmzFile)
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
%     The output polygon matrix for the output boundaryq, with each row
%     being one vertex in the form of UTM (x, y).
%
% Yaguang Zhang, Purdue, 06/03/2020

utmXyBoundary = polyshape();
utmZone = '';

kmzStruct = kmz2struct(dirToKmzFile);
polyonsLonLatStruct = kmzStruct(strcmp({kmzStruct.Geometry}, 'Polygon'));

numOfPolygons = length(polyonsLonLatStruct);
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
        
        curPolygon = polyshape(curXs, curYs);
        if ~isempty(utmXyBoundary.Vertices)
            utmXyBoundary = union(utmXyBoundary, curPolygon);
        else
            utmXyBoundary = curPolygon;
        end
    end
end

utmXyBoundary = utmXyBoundary.Vertices;

% Make sure the output vertices are clockwise and the polygon is closed.
assert(ispolycw(utmXyBoundary(:,1),utmXyBoundary(:,2)), ...
    'The output vertices are not clockwise!')
utmXyBoundary(end+1, :) = utmXyBoundary(1, :);

end
% EOF