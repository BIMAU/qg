function [errs, nums] = plot_experiment(varargin)
% plot results of a parallel experiment

    switch nargin
      case 0
        exp_name = 'test';
        exp_type = 'parallel';
        procs    = 96;

      case 2
        exp_name = varargin{1};
        exp_type = varargin{2};
        procs    = 96;

      case 3
        exp_name = varargin{1};
        exp_type = varargin{2};
        procs    = varargin{3};

      otherwise
        error('Unexpected input');
    end

    dir = ['data/experiments/', exp_name, '/', exp_type, '/'];

    fileNames = cell(procs,1);
    for i = 1:procs
        fileNames{i} = sprintf([dir, 'results_%d.mat'],i-1);
    end

    data = load(fileNames{1});
    trials = size(data.num_predicted, 2);
    n = size(data.num_predicted, 1); % ensemble size
    errs = cell(procs, trials);
    nums  = nan(n, trials);
    ind_vis = zeros(procs, ceil(n / procs));

    for d = 1:procs
        data = load(fileNames{d});
        for i = data.my_inds
            for j = 1:trials
                errs{i, j} = data.errs{i, j};
                num = data.num_predicted(i, j);
                if num > 0
                    nums(i, j) = num;
                end
            end
        end
    end

    assert(i == size(data.num_predicted, 1), ...
           'failed assertion, probably wrong procs');

    my_boxplot(nums);

    % % alternative error bound:
    % nums  = nan(n, trials);
    % subplot(2,1,2)
    % run = 1;
    % for j = 1:trials
    %     for i = 1:n
    %         num = find(errs{i,j}>4.0,1);
    %         if ~isempty(num)
    %             nums(i,j) = num
    %         end
    %     end
    % end
    % my_boxplot(nums);
end