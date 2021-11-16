%CLEANFCCLOGVARS A helper snippet to clean the FCC log results.
%
% Yaguang Zhang, Purdue, 11/15/2021

osVersions(flagsEleToRemove) = [];
carriers(flagsEleToRemove) = [];
downSpeedsBps(flagsEleToRemove) = [];
upSpeedsBps(flagsEleToRemove) = [];
downEndLatLons(flagsEleToRemove, :) = [];
upEndLatLons(flagsEleToRemove, :) = [];
testIds(flagsEleToRemove) = [];

% EOF