%GUARDMEMAVAILONLINUX A snippet to check the memory available on a Linux
%machine and make sure it is not too low.
%
% It will restart the parallel pool if availabe memory is less than a
% threshold.
%
% Yaguang Zhang, Purdue, 03/02/2022

% Only do the check if a parallel pool is present.
flagPoolRestarted = false;
curCluster = gcp('nocreate');
if ~isempty(curCluster)
    % Only do the check for Linux machines, e.g., the Purdue cluster
    % Frankie, where the simulation runs in the production mode.
    if isunix
        [checkMemStatus, availMemInGiBStr] ...
            = system(['awk ''/MemAvailable/ ', ...
            '{ printf "%.3f \n", $2/1024/1024 }'' /proc/meminfo']);
        assert(~checkMemStatus, 'Unable to check memory usage!');
        availMemInGiB = str2double(availMemInGiBStr);

        minAvailMemInGiBRequired = 32; % We want at least 32 GiB arailable.
        if availMemInGiB < minAvailMemInGiBRequired
            disp(' ');
            disp(['Warning: Alailable memory (', ...
                num2str(availMemInGiB), ...
                ' GiB) is less than the required amout (', ...
                num2str(minAvailMemInGiBRequired), ' GiB)!']);

            disp(' ');
            disp('Restarting parallel pool...');

            tic;
            delete(curCluster);
            curCluster = gcp;
            if exist('numOfWorkersInCurCluster', 'var')
                assert(numOfWorkersInCurCluster==curCluster.NumWorkers, ...
                    'New cluster has a different number of workers!');
            end
            flagPoolRestarted = true;
            toc;

            disp('Done!');
        end
    end
end

% EOF