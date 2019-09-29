function [ ] = parsave(absPathToSave, varargin) %#ok<INUSD>
%PARSAVE A helper function to enable save in parfor loops.
%
% Inputs:
%   - absPathToSave
%     The absolute path for saving the variables.
%   - varargin
%     The variables to be saved.
% Yaguang Zhang, Purdue, 09/10/2019

varLabels = strings(1, nargin-1);
for idxVar = 1:(nargin-1)
    curVarLabel = inputname(idxVar+1);
    eval([curVarLabel, ' = varargin{', num2str(idxVar), '};']);
    varLabels(idxVar) = curVarLabel;
end

eval(['save(absPathToSave, "', ...
    table2array(join(varLabels, '", "')), '");']);

end
% EOF