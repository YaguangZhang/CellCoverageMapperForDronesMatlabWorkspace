function [ ] = parsave(absPathToSave, varargin) %#ok<INUSD>
%PARSAVE A helper function to enable save in parfor loops.
%
% Inputs:
%   - absPathToSave
%     The absolute path for saving the variables.
%   - varargin
%     The variables to be saved.
%
% Note:
%   all the inputs should be independent and complete variables, that is,
%   they should have their own variable names and fully contained in
%   themselves (to save part of a matrix or part of a cell like
%   matToSave(1) or varToSave{1} are not supported).
%
% Example:
%   parsave(curFullPathToSaveResults, varToSave, moreVarToSave);
%
% Yaguang Zhang, Purdue, 09/10/2019

varLabels = strings(1, nargin-1);
for idxVar = 1:(nargin-1)
    curVarLabel = inputname(idxVar+1);
    eval([curVarLabel, ' = varargin{', num2str(idxVar), '};']);
    varLabels(idxVar) = curVarLabel;
end

eval(['save(absPathToSave, "', ...
    convertStringsToChars(join(varLabels, '", "')), '");']);

end
% EOF