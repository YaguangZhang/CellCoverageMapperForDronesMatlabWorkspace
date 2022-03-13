function [fsplVs, latLonFsplMap] = genOptimisticFsplMap(simState)
%GENOPTIMISTICFSPLMAP Compute the best available (min) FSPL on a 2D plane.
%
% Yaguang Zhang, Purdue, 03/11/2022

[~, dists] = dsearchn(simState.CellAntsXyhEffective(:,1:2), ...
    simState.mapGridXYPts);
fsplVs = fspl(dists, simConfigs.CARRIER_WAVELENGTH_IN_M);
latLonFsplMap = [simState.mapGridLatLonPts, fsplVs];

end
% EOF