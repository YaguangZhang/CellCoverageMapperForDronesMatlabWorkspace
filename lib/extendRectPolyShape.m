function [ extendedPolyShape ] ...
    = extendRectPolyShape(rectPolyShape, amount, flagDebug)
%EXTENDRECTPOLYSHAPE Extend the input rectangle polyshape. We will move
%each vertex outward by the specified amount on both the x and y
%directions.
%
% Inputs:
%   - rectPolyShape
%     The input polyshape. We will only use its first four vertices, so it
%     can be either closed (with a redundent 5th vertex being the same as
%     its 1st one) or not.
%   - amount
%     The amount to shift for each vertex along the x (y) direction.
%   - flagDebug
%     Optional. Default to false. Set this to true to generate a debug
%     figure for comparing the input and output polyshapes.
%
% Output:
%   - extendedPolyShape
%     The extended polyshape. We will output the closed version (with a
%     redundent 5th vertex being the same as its 1st one) .
%
% Yaguang Zhang, Purdue, 06/04/2020

if ~exist('flagDebug', 'var')
    flagDebug = false;
end

% Make sure the vertices of the input polygon are organized as expected.
[numVertices, dimension] = size(rectPolyShape.Vertices);
assert((numVertices==4 | numVertices==5) & dimension==2, ...
    'The vertices of the input polygon are not organized as expected!');

% Use the centroid to determine to which direction we should shift each
% vertex.
centroid = mean(rectPolyShape.Vertices);
extendedPolyShapeVertices = nan(4,2);
for idxVertex = 1:4
    curVertex = rectPolyShape.Vertices(idxVertex,:);
    extendedPolyShapeVertices(idxVertex,:) = ...
        curVertex+sign(curVertex-centroid).*amount;
end
extendedPolyShape = polyshape( ...
    [extendedPolyShapeVertices(:,1); ...
    extendedPolyShapeVertices(1,1)], ...
    [extendedPolyShapeVertices(:,2); ...
    extendedPolyShapeVertices(1,2)]);

if flagDebug
    figure; hold on;
    hExtendedPoly = plot(extendedPolyShape);
    hInputPoly = plot(rectPolyShape, ...
        'FaceColor', 'white', 'FaceAlpha', 0.5);
    legend([hInputPoly, hExtendedPoly], 'Input Polygon', 'Extended')
    axis equal;
    grid on;
end
% EOF