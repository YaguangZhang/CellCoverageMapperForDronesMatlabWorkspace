function [projectionName] ...
    = validateProjectionByParameter(projectionParameters)
%VALIDATEPROJECTIONBYPARAMETER Validate the projection coordinate system
%via parameters.
%
% Input:
%   - projectionParameters
%     A struct/object with fields LatitudeOfNaturalOrigin,
%     LongitudeOfNaturalOrigin, ScaleFactorAtNaturalOrigin, FalseEasting,
%     and FalseNorthing.
% Output:
%   - projectionName
%     A string representing the projection coordinate system. We will try
%     to find a known projection stored in the file
%           parametersForKnownProjections.mat
%     and if none is find, "unknown" will be set as the output.
%
% Yaguang Zhang, Purdue, 02/01/2021

fieldsToCheck = {'LatitudeOfNaturalOrigin', 'LongitudeOfNaturalOrigin', ...
    'ScaleFactorAtNaturalOrigin', 'FalseEasting', 'FalseNorthing'};
maxErrorAllowedForFields = [10^(-10), 10^(-10), 10^(-8), 10^(-2), 10^(-2)];

projectionName = 'unknown';

knownProjections = load(fullfile( ...
    fileparts(mfilename('fullpath')), ...
    'parametersForKnownProjections.mat'));
knownProjectionNames = fieldnames(knownProjections);

for idxProj = 1:length(knownProjectionNames)
    curProjName = knownProjectionNames{idxProj};
    curProj = knownProjections.(curProjName);
    
    flagProjFound = true;
    for idxField = 1:length(fieldsToCheck)
        curField = fieldsToCheck{idxField};
        if curProj.(curField) - projectionParameters.(curField) ...
                > maxErrorAllowedForFields(idxField)
            flagProjFound = false;
        end
    end
    
    if flagProjFound
        projectionName = curProjName;
        break;
    end
end

end
% EOF