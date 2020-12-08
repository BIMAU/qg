function [] = main(varargin)
% wrapper for experiment using the parpool
    switch nargin
      case 0
        pid      = 0;
        procs    = 1;

      case 2
        pid     = arg2value(varargin{1});
        procs   = arg2value(varargin{2});

      case 3
        pid     = arg2value(varargin{1});
        procs   = arg2value(varargin{2});
        threads = arg2value(varargin{3});

      otherwise
        error('Unexpected input');
    end

    fprintf('load training data...\n'); tic;
    fname_base = 'N128-N64_ff2_Re1.0e+04-Re1.0e+02_Tstart159_Tend187';
    trdata = load(['data/training/', fname_base, '.mat']);
    fprintf('load training data... done (%fs)\n', toc);
    
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        fprintf('deleting existing parallel pool\n');
        delete(gcp('nocreate'))
    end

    if exist('threads','var')
        if threads > 1
            poolobj = parpool(threads);        
        end
    else
        poolobj = parpool();
        threads = poolobj.NumWorkers;
    end

    fprintf('--------------------------------------------\n')
    fprintf('-----     main  \n');
    fprintf('-----    procs: %d\n', procs);
    fprintf('-----      pid: %d\n', pid);
    fprintf('-----  threads: %d\n', threads);    

    if threads > 1
        parfor t = 0:threads-1
            task = getCurrentTask();
            clst = getCurrentCluster();
            fprintf(['| parfor      ID     GID  GPROCS \n', ...
                     '|    %3d     %3d     %3d     %3d \n'], ...
                    t, task.ID-1, pid*threads+t, procs*threads);
            experiment(pid*threads+t, procs*threads, trdata);
        end
        delete(poolobj);
    else
        experiment(pid, procs, trdata);
    end
end
