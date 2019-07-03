% LOADEHATACPPLIB A short script to load the NTIA C++ eHata model.
%
% Note: this script will reset the current search path.
%
% Yaguang Zhang, Purdue, 06/19/2019

if ~libisloaded('ehata')
    if ~(strcmpi(computer, 'PCWIN64') || strcmpi(computer, 'GLNXA64'))
        warning(['The NTIA C++ eHata library was setup only ', ...
            'for x64 Windows/Linux machines!']);
    end
    
    try
        restoredefaultpath;
        
        if ispc
            % The .dll version.
            addpath(fullfile('lib', 'ext', 'eHataNtia'));
        elseif isunix
            % The .so version.
            addpath(fullfile('lib', 'ext', 'eHataNtiaLinux'));
        end
        
        loadlibrary('ehata');
    catch e
        switch e.identifier
            case 'MATLAB:loadlibrary:FileNotFound'
                disp(e.message);
                error(['Restarting Matlab and loading the lib first ', ...
                    'may resolve this issue.']);
            otherwise
                error(e.message);
        end
    end
end
% EOF