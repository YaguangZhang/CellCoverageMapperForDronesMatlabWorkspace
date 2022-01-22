% PREPARESIMULATIONENV Set environment for simulations.
%
% Yaguang Zhang, Purdue, 09/10/2019

% For logging computation time.
dataTimeStrStart = datestr(datetime('now'));
timerValueStart = tic;

% For separating hints from different sections of the program.
fileNameHintRuler = hintScriptName(curFileName);

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
setPath;

% The format to use for displaying datetime via datestr.
datetimeFormat = 'yyyy/mm/dd HH:MM:ss';

% EOF