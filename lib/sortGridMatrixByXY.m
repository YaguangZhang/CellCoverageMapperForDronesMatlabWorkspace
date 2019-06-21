function [ XsSorted, YsSorted, gridMatSorted ] ...
    = sortGridMatrixByXY(Xs, Ys, gridMat)
%SORTGRIDMATRIXBYXY Sort the data on a grid represented by gridMat
%according to the grid labels Xs and Ys.
%
% Inputs:
%   - Xs, Ys
%     Column vectors. The labels for the horizontal axis (the columns) and
%     the vertical axis (the rows) of the grid data, respectively.
%   - gridMat
%     The grid data stored as a matrix. The (i,j) element corresponds to
%     the location [Xs(j), Ys(i)].
%
% Outputs:
%   - XsSorted, YsSorted, gridMatSorted
%     The sorted results for the inputs so that the new X labels are
%     monotonically increasing while the new Y labels are monotonically
%     decreasing.
%
% Yaguang Zhang, Purdue, 06/11/2019

gridMatSorted = [Xs'; gridMat];
gridMatSorted = sortrows(gridMatSorted', 1)';
XsSorted = gridMatSorted(1, :)';
gridMatSorted = gridMatSorted(2:end, :);

gridMatSorted = [Ys gridMatSorted];
gridMatSorted = sortrows(gridMatSorted, 1, 'descend');
YsSorted = gridMatSorted(:, 1);
gridMatSorted = gridMatSorted(:, 2:end);

end
% EOF