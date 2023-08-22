%GETRATIOOFSHORTLINKS From variables simConfigs and simState, compuate the
%ratio of links are shorter than a given distance ratio.
%
% Note: Any link with a bigger length than
% simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M is ignored. 
% 
% Yaguang Zhang, Purdue, 03/25/2023

MIN_LENGTH_IN_KM_OF_LONG_LINKS = 1;

linkLengthsInM = pdist2( ...
    simState.CellAntsXyhEffective(:, 1:2), ... Tx locs.
    simState.mapGridXYPts);                  % Rx locs.
linkLengthsInM = linkLengthsInM(:);
linkLengthsInM = linkLengthsInM( ...
    linkLengthsInM<=(simConfigs.MAX_CELL_COVERAGE_RADIUS_IN_M));

ratioOfShortLinks = sum( ...
    linkLengthsInM<(MIN_LENGTH_IN_KM_OF_LONG_LINKS*1000)) ...
    /length(linkLengthsInM);

disp(['ratioOfShortLinks = ', num2str(ratioOfShortLinks*100), '%']);
% EOF