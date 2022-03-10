function h = removeDrivesFromHistory()
% removeDrivesFromHistory()
%   Creates a GUI that allows the user to clean up the current folder history by 
%   removing selected paths. This is necessary to fix a common problem associated
%   with network locations on the current folder path history that are no longer 
%   available.  When an undefined function or variable error is elicited, Matlab 
%   will search the current folder history and when it encounters a path to an 
%   unavailable location, the program may stall for an unbearable amount of time 
%   before throwing the error.  See the links below for discussions on this matter 
%   in TMW forum.
%
%   Before making any changes, a backup of the current folder history is saved to a 
%   text file named pathHistoryBackup.txt and is saved in tempdir().
%
%   The GUI will list the paths that appear in the current folder history and 
%   contains the following button options:
%   * [Show current] displays the current folder path history in the command window.
%   * [Show backup] opens the backup text file showing the most recent backup.
%   * [Restore history] restores Matlab's current folder history with the backup.
%   * [Help] displays the documentation for this function.
%   * [Close] Closes the GUI without making any additional changes.
%   * [Remove] Updates the backup file with the current folder history and removes
%       all paths selected in the GUI. Select docmultiple paths by holding down the 
%       Ctrl key. 
%
% Discussion links
% <a href = "https://www.mathworks.com/matlabcentral/answers/395876#answer_315892">Undefined function error is very slow</a>
% <a href = "https://www.mathworks.com/matlabcentral/answers/364153">Waiting for Undefined function or variable error</a>
% <a href = "https://www.mathworks.com/matlabcentral/answers/412972">Matlab stalls on undefined variables</a>
% 
% Author: <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776-adam-danz">Adam Danz</a>
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/72519">removeDrivesFromHis​tory</a>

% Copyright (c) 2014-2021, Adam Danz
% All rights reserved

% Version history
% vs 1.0.0  190824  Initial release. 
% vs 1.0.1  190824  Fixed typo 'mfilename' -> mfilename
% vs 1.0.2  190825  Increased fontsize of dlgs; added greeting msg.
% vs 2.0.0  190828  Now users can select specific paths instead of drives; documentation improved.
% vs 2.1.0  190829  Shows current cd() at bottom of GUI; notification when removing curr cd().
% vs 2.1.1  191122  Fixed typo in the "Your current directory is..." message box.
% vs 2.2.0  210317  Fixed Matlab version check that failed in r21a; Close/reopen gui in same position; 
%                   Lists paths in cmd window as categories rather than cellstr.

%% Prepare for the GUI
persistent initialMessage

% This solution is not applicable to releases prior to r2017b (ver 9.3) [1]
if verLessThan('Matlab', '9.3')
    error('Requires Matlab r2017b or higher.')
end

% Input/output validity
narginchk(0,0)
nargoutchk(0,1)

% Issue a once-and-done question box that requires the user to confirm that they
% know what they are doing.
opts = struct('WindowStyle','non-modal','Interpreter','tex');
if isempty(initialMessage)
    msg = sprintf(['\\fontsize{12}Welcome!\n\nThis program will give you options to remove  '...
        'selected paths from your current folder history in order to avoid an obnoxious '...
        'delay that often preceeds the error message caused by an undefined function '...
        'or variable.  Press the "Help" button in the GUI for more information.']);
    mbh = msgbox(msg,mfilename,'help',opts);
    uiwait(mbh)
    initialMessage = true;
end

% Get history; see footnotes [1,2,3]
S = settings();
try
    dat.pathList = S.matlab.desktop.currentfolder.History.PersonalValue;
catch ME
    if strcmp(ME.identifier,'MATLAB:settings:config:UndefinedSettingValueForLevel')
        % This happens when the PersonalValue is empty.
        msb = msgbox(['\fontsize{12}',ME.message, ' (i.e. your current folder history is empty)'], mfilename, 'help',opts);
        uiwait(msb)
        dat.pathList = {};
        dat.drives = {};
    end
end

% Get drive letter from each path (for future development, maybe).
dat.pathList(cellfun(@isempty,dat.pathList)) = [];
% dat.drives = regexp(dat.pathList,'(?i)([a-z]):','tokens');
% dat.drives = cellfun(@(x)x{:},dat.drives);

%% Create the GUI
% First determine if GUI already exists.
search = findall(0,'type','Figure','tag',[mfilename,'GUI']);
if ~isempty(search)
    delete(search)
end
gui.colors = [
    0.148 0.148 0.148; ... %background, minshaft
    0.979 0.500 0.445; ... %Text, salmon
    0.660 0.660 0.660];    %bacgkround2, dark gray
gui.fh = figure('Name',[mfilename,'.m'],'Color',gui.colors(1,:),'menubar','none','toolbar','none',...
    'units','pixels','NumberTitle','off','HandleVisibility','off','tag',[mfilename,'GUI']);
gui.fh.Position(3:4) = [650,320];
movegui(gui.fh)
gui.title = uicontrol(gui.fh,'Style','text','Units','Normalize','Position',[0,.80,1,.15], ...
    'String',{'Select one or more locations to remove from history'},'FontSize',14', ...
    'FontWeight','Bold','BackgroundColor',gui.colors(1,:),'ForegroundColor',gui.colors(2,:));
gui.listbox = uicontrol(gui.fh,'Style','ListBox','Units','Normalize','Position',[.05,.09,.75,.68],...
    'String',dat.pathList,'FontSize',10,'BackgroundColor',gui.colors(3,:),'Min',1,'Max',1000);
gui.instruct = uicontrol(gui.fh,'Style','text','Units','Normalize','Position',[.05,.78,.75,.05], ...
    'String','Select 1 or more for removal','FontSize',8','BackgroundColor',gui.colors(1,:),...
    'ForegroundColor',gui.colors(2,:),'HorizontalAlignment','Left');
gui.currCDtxt = uicontrol(gui.fh,'Style','text','Units','Normalize','Position',[.05,.02,.75,.05], ...
    'String',sprintf('cd: %s',cd()),'FontSize',8','BackgroundColor',gui.colors(1,:),...
    'ForegroundColor',gui.colors(3,:),'FontName','monospaced','HorizontalAlignment','Left');
gui.HxTxt = uicontrol(gui.fh,'Style','text','Units','Normalize','Position',[.82,.78,.16,.05], ...
    'String','History','FontSize',8','BackgroundColor',gui.colors(1,:),...
    'ForegroundColor',gui.colors(2,:));
gui.showHxButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.71,.16,.06],...
    'String','Show current','BackgroundColor',gui.colors(3,:),'Fontsize', 10, 'callback', {@showCurrentHxFcn, dat.pathList});
gui.showBUHxButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.62,.16,.06],...
    'String','Show backup','BackgroundColor',gui.colors(3,:),'Fontsize', 10, 'callback', @showPathHxBUFcn);
gui.restoreButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.53,.16,.06],...
    'String','Restore history','BackgroundColor',gui.colors(3,:),'Fontsize', 10, 'callback', {@restoreHistoryFcn,gui});
gui.helpButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.29,.16,.06],...
    'String','Help','BackgroundColor',gui.colors(3,:),'Fontsize', 10, 'callback', @(hObj,event)doc(mfilename));
gui.CloseButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.20,.16,.06],...
    'String','Close','BackgroundColor',gui.colors(3,:),'Fontsize', 10, 'callback', @(hObj,event)close(hObj.Parent));
gui.goButton = uicontrol(gui.fh,'Style','PushButton','Units','Normalize','Position',[.82,.09,.16,.08],...
    'String','Remove','BackgroundColor',gui.colors(3,:),'Fontsize', 14, 'FontWeight','bold',...
    'callback', {@removeDrivesFromHistoryFcn,gui});

% Store drive lists in GUI
guidata(gui.fh,dat);
if nargout > 0
    h = gui.fh;
end

%% Callback function (in order of appearance in GUI)
function showCurrentHxFcn(~, ~,pathList)
% Responds to the "show current history" button.
T = table(pathList(:),'VariableNames',{'Current_History'});
if isempty(T)
    opts = struct('WindowStyle','non-modal','Interpreter','tex');
    msgbox('\fontsize{12}Matlab''s current folder history is currently empty.', mfilename, 'help',opts);
else
    % Display history in command window.
    fprintf('\n')
    disp(T)
    commandwindow();
end

function showPathHxBUFcn(~,~)
% Responds to the "Show history backup" button
td = tempdir();
ff = fullfile(td,'pathHistoryBackup.txt');
if exist(ff,'file') == 2
    winopen(ff);
else
    opts = struct('WindowStyle','non-modal','Interpreter','tex');
    errordlg(['\fontsize{12}Backup file not found. A backup file will be created '...
        'prior to making any changes to your current folder history.'], mfilename, opts)
end

function restoreHistoryFcn(~, ~, gui)
% Responds to the restore button.  It will search for a backup file and if
% it exists, it will restore the history to the previous backup.
% Search for backup file
td = tempdir();
ff = fullfile(td,'pathHistoryBackup.txt');
% If backfup file doesn't exist, throw error.
if exist(ff,'file') ~= 2
    opts = struct('WindowStyle','non-modal','Interpreter','tex');
    errordlg(['\fontsize{12}Backup file not found. A backup file will be created '...
        'prior to making any changes to your current folder history.'], mfilename, opts)
    return
end
% Read file
txt = strtrim(strsplit(fileread(ff),'\n'));
% Remove first 5 rows of text
txt(1:5) = [];
% Remove empty cells
txt(cellfun(@isempty,txt)) = [];
% assign to history
S = settings();
S.matlab.desktop.currentfolder.History.PersonalValue = txt;
% Close GUI and reopen it so the listbox is regenerated
closeAndReopenGUI(gui)

function removeDrivesFromHistoryFcn(~, ~, gui)
% Responds to pressing the 'Remove' button in the GUI created above;
% Get info from GUI
% Determine which locations are selected
driveOptions = gui.listbox.String;
opts = struct('WindowStyle','non-modal','Interpreter','tex');
if isempty(driveOptions)
    msgbox('\fontsize{12}Matlab''s current folder history is currently empty.', mfilename, 'help',opts);
    return
end
drives2rm = driveOptions(gui.listbox.Value);
% If no drives were selected, quit and warn.
if isempty(drives2rm)
    msgbox('\fontsize{12}Nothing was selected to be removed.',mfilename,'error',opts)
    return
end
% Backup current history to temp file (overwrite)
createBackupOfPathHxFcn(gui);
% Get gui data
dat = guidata(gui.fh);

% If the current directory is part of any drives flagged for removal, 
% change the cd to 1st user path
currcd = cd();
if any(strcmp(drives2rm, currcd))
    opts.Default = 'Continue'; 
    msg = sprintf(['Your current directory is one of the paths selected for removal. '...
        'The current directory will be reset to the directory listed first in the search '...
        'path.\n\nCurrent directory:  %s\nWill be set to:  %s'],regexprep(currcd,'\','\\\'),...
        regexprep(userpath,'\','\\\'));
    resp = questdlg(['\fontsize{12}',msg], mfilename, 'Continue','Abort',opts);
    if isempty(resp) || strcmpi(resp,'Abort')
        return
    else
        cd(userpath);
    end
end

% Remove selected drives from list
rmidx = ismember(dat.pathList, drives2rm);
dat.pathList(rmidx) = [];
% reassign history
s = settings();
s.matlab.desktop.currentfolder.History.PersonalValue = dat.pathList;

% Close GUI & reopen to refresh list
fprintf('\nThe following paths have been removed from the current folder history:\n')
disp(categorical(drives2rm))
if nargin < 4
    closeAndReopenGUI(gui)
end

function createBackupOfPathHxFcn(gui)
% The current path history will be written to file.  The file will be created or overwritten.
dat = guidata(gui.fh);
% Create filename and path
td = tempdir();
ff = fullfile(td,'pathHistoryBackup.txt'); %restoreHistoryFcn() & showPathHxBUFcn() relies on this filename and location
% Create header
header = {sprintf('Matlab current folder history\nCreated in %s.m on %s.\nBackup: %s',mfilename,datestr(now),ff)};
content = [header;{''};{''}; dat.pathList(:)]; % must start on line 6!
% Write to text file in temp dir
filePh = fopen(ff,'wt+');
fprintf(filePh,'%s\n',content{:});
fclose(filePh);
disp(['Saved backup of current folder history (','<a href="matlab: winopen(''',ff,''') ">', ...
    'pathHistoryBackup.txt','</a>'...
    ') to ','<a href="matlab: winopen(''',td,''') ">', 'temp directory','</a>'])

function closeAndReopenGUI(gui)
% Reset GUI by close/re-open in same position. 
guipos = gui.fh.Position;
close(gui.fh)
h = removeDrivesFromHistory();
h.Position = guipos;

%% notes
%[1] Documentation suggests that settings() came out in r2018a but it was also available on my r2017b release.
%      Futhermore, I don't think the findUnlicensedFunctions() was a problem before then anyway.
%[2] https://www.mathworks.com/matlabcentral/answers/395876-undefined-function-error-is-very-slow-to-occur#answer_315892
%[3] This function can be found in ...\toolbox\matlab\configtools\settings.m
%     It is equivalent to running: S = matlab.settings.internal.settings;
