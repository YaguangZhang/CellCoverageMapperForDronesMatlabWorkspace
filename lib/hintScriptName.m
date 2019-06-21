function [ fileNameHintRuler ] ...
    = hintScriptName(curFileName)
% HINTSCRIPTNAME Display the script name stored in curFileName as a hint.
%
% Yaguang Zhang, Purdue, 06/19/2019

fileNameHintRuler = [' ', repmat('-', 1, length(curFileName)+2), ' '];
disp(fileNameHintRuler)
disp(['  ', curFileName, '  '])
disp(fileNameHintRuler)

end
% EOF