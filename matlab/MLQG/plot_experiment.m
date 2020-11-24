function [] = plot_experiment(varargin)
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
    
    dir      = ['data/experiments/', exp_name, '/', exp_type, '/'];

    fileNames = cell(procs,1);
    for i = 1:procs
        fileNames{i} = sprintf([dir, 'results_%d.mat'],i-1);
    end

    data = load(fileNames{1});
    assert(procs == size(data.num_predicted, 1));
    trials = size(data.num_predicted, 2);
    errs = cell(procs, trials);
    nums = zeros(procs, trials);

    for i = 1:procs
        data = load(fileNames{i});

        for j = 1:trials
            errs{i, j} = data.errs{i, j};
            nums(i, j) = data.num_predicted(i, j);
        end
    end

    my_boxplot(nums)
    axis off
end

