function [errs, nums, pids, metadata, predictions, truths] = gather_plotdata(varargin)

    switch nargin
      case 1
        dir     = varargin{1};
        
        % number of procs = number of files in dir:
        [~, fc] = system(['ls ', dir, ' -1 | wc -l']);
        procs   = str2num(fc);

      case 2
        dir   = varargin{1};
        procs = varargin{2};

      case 3
        dir      = varargin{1};
        procs    = varargin{2};
        tr_range = varargin{3};

      otherwise
        error('Unexpected input');
    end

    if procs > 1
        serial = false;
    else
        serial = true;
    end
    
        

    fileNames = cell(procs,1);

    if serial
        fileNames{1} = sprintf([dir, 'results.mat']);
    else
        for i = 1:procs
            fileNames{i} = sprintf([dir, 'results_%d.mat'], i-1);
            if ~exist(fileNames{i}, 'file')
                fprintf('%s does not exist\n', fileNames{i});
                fileNames{i} = 'failedproc';
            end
        end
    end

    initialize = true;

    for d = 1:procs
        if strcmp(fileNames{d}, 'failedproc')
            continue;
        end
        data = load(fileNames{d});

        if initialize
            trials  = size(data.num_predicted, 2);

            if exist('tr_range', 'var')
                trials = numel(tr_range);
            else
                tr_range = 1:trials;
            end

            n       = size(data.num_predicted, 1); % ensemble size
            errs    = cell(procs, trials);
            pids    = cell(procs, trials);

            predictions = cell(procs, trials);
            truths      = cell(procs, trials);

            nums    = nan(n, trials);
            ind_vis = zeros(procs, ceil(n / procs));
            initialize = false;
        end

        for i = data.my_inds
            for j = 1:trials
                errs{i, j} = data.errs{i, tr_range(j)};
                pids{i, j} = d;

                predictions{i, j} = data.predictions{i, j};
                truths{i, j}      = data.truths{i, j};

                num = data.num_predicted(i, tr_range(j));
                if num > 0
                    nums(i, j) = num;
                end
            end
        end
    end

    assert(i == size(data.num_predicted, 1), ...
           'failed assertion, probably wrong procs');

    % export additional labels for plotting
    metadata = struct();
    importlabs = {'xlab', ...
                  'ylab', ...
                  'hyp_range',...
                  'hyp',...
                  'exp_id',...
                  'exp_ind',...
                  'run_pars',...
                  'esn_pars',...
                 };

    for l = 1:numel(importlabs)
        lab = importlabs{l};
        if isfield(data, lab)
            metadata.(lab) = data.(lab);
        end
    end
end