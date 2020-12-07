function [errs, nums, pids, labels] = gather_plotdata(varargin)

    switch nargin
      case 3
        exp_name = varargin{1};
        exp_type = varargin{2};
        procs    = varargin{3};
      case 4
        exp_name = varargin{1};
        exp_type = varargin{2};
        procs    = varargin{3};
        tr_range = varargin{4};
      otherwise
        error('Unexpected input');
    end

    dir = ['data/experiments/', exp_name, '/', exp_type, '/'];

    fileNames = cell(procs,1);
    for i = 1:procs
        fileNames{i} = sprintf([dir, 'results_%d.mat'],i-1);
    end

    data    = load(fileNames{1});
    trials  = size(data.num_predicted, 2);


    if exist('tr_range', 'var')
        trials = numel(tr_range);
    else
        tr_range = 1:trials;
    end

    n       = size(data.num_predicted, 1); % ensemble size
    errs    = cell(procs, trials);
    pids    = cell(procs, trials);
    nums    = nan(n, trials);
    ind_vis = zeros(procs, ceil(n / procs));

    for d = 1:procs
        data = load(fileNames{d});
        for i = data.my_inds
            for j = 1:trials
                errs{i, j} = data.errs{i, tr_range(j)};
                pids{i, j} = d;
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
    labels = struct();
    importlabs = {'xlab', ...
                  'ylab', ...
                  'hyp_range'};

    for l = 1:numel(importlabs)
        lab = importlabs{l};
        if isfield(data, lab)
            labels.(lab) = data.(lab);
        end
    end
end