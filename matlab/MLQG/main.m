function [] = main(varargin)
% wrapper for experiment using the parpool
    switch nargin
      case 0
        node     = 0;
        nodes    = 1;

      case 2
        node    = arg2value(varargin{1});
        nodes   = arg2value(varargin{2});

      case 3
        node    = arg2value(varargin{1});
        nodes   = arg2value(varargin{2});
        threads = arg2value(varargin{3});

      otherwise
        error('Unexpected input');
    end

    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(gcp('nocreate'))
    end

    if exist('threads','var')
        poolobj = parpool(threads);
    else
        poolobj = parpool();
    end

    threads = poolobj.NumWorkers;
    fprintf('--------------------------------------------\n')
    fprintf('-----         main  \n', nodes);
    fprintf('-----  number of nodes: %d\n', nodes);
    fprintf('-----  this node: %d\n', node);
    fprintf('-----  number of threads: %d\n', threads);

    parfor t = 0:threads-1
        experiment(node*threads+t, nodes*threads);
    end

    delete(poolobj);
end
