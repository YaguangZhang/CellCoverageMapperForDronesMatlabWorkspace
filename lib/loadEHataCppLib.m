% LOADEHATACPPLIB A short script to load the NTIA C++ eHata model.
%
% Note: this script will reset the current search path.
%
% Yaguang Zhang, Purdue, 06/19/2019

if ~libisloaded('ehata')
    if ~strcmpi(computer, 'PCWIN64')
        warning(['The NTIA C++ eHata library was setup only ', ...
            'for x64 Windows machines!']);
    end
    
    try
        restoredefaultpath;
        addpath(fullfile('lib', 'ext', 'eHataNtia'));
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