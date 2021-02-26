function [] = plot_experiment(varargin) 
    
    switch nargin
      case 2
        nums = varargin{1};
        mdat = varargin{2};
      case 0
        
        dir = ['data/experiments/',  'NR_4000-16000_ESN1_MDL0', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_500-16000_ESN1_MDL0', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_1000-4000_RF_8-16_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data`/experiments/', 'NR_1000-8000_FA_0.1-1_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_4000-8000_RF_1-16_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_500-8000_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_1000-8000_BS_1-4_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_1000-8000_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_1000-8000_SP_1000-5000_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_1000-32000_ESN1_MDL1', '/', 'parallel', '/'];
        dir = ['data/experiments/',  'NR_500-500_ESN1_MDL1', '/', 'serial', '/'];

        [errs,nums,pids,mdat] = gather_plotdata(dir, 1);
      otherwise
        error('Unexpected input')
    end

    % ignore cols with nans
    % plotJds = find(~isnan(sum(nums,1)));
    % nums = nums(:,plotJds);

    % ignore rows with nans
    % plotIds = find(~isnan(sum(nums,2)));
    % nums = nums(plotIds,:);

    % number of experiments
    exp_ind = mdat.exp_ind;
    Nexp    = numel(exp_ind);

    labels  = [];
    Nvalues = [];

    for i = 1:Nexp
        labels{i} = mdat.hyp_range(exp_ind{i}, :);
        Nvalues(i) = numel(unique(labels{i}));
    end

    [~, I] = sort(Nvalues, 'descend');
    xlab_index = I(1); % x label corresponds to parameter with largest number of values
    maxValues  = Nvalues(xlab_index);
    Ntotal     = size(nums,2);
    Nboxplots  = Ntotal/maxValues;
    assert(Nboxplots == round(Nboxplots))

    % plot results
    f = [];
    clrs = lines(Nboxplots);

    H = mdat.hyp_range(exp_ind{xlab_index}, :);
    range1 = 1:numel(H);
    if Nexp == 2
        M = reshape(1:numel(H), Nvalues(1), Nvalues(2));
        if H(2) ~= H(1)
            range1 = M(:);
            M = M';
            range2 = M(:);
        else
            range2 = M(:);
            M = M';
            range1 = M(:);
        end
    end

    for i = 1:Nboxplots
        subrange = range1((i-1)*maxValues+1:i*maxValues);
        f{i} = my_boxplot(nums(:, subrange), {clrs(i,:), clrs(i,:)});
        grid on;
        xticklabels([]);
        hold on
    end
    hold off

    xticklabels(mdat.hyp_range(exp_ind{xlab_index}, subrange));
    xtickangle(45);
    xlabel(mdat.xlab{xlab_index});
    ylabel(mdat.ylab);

    % for combined experiments and multiple boxplots we need a legend
    if Nexp == 2
        str = cell(Nboxplots,1);
        for i = 1:Nboxplots
            value = mdat.hyp_range(exp_ind{I(2)}, range2(i));
            str{i} = sprintf('%s: %1.1e', mdat.xlab{I(2)}, value);
        end
        legend([f{:}], str, 'location', 'north')
    end

    % create description
    descr = create_description(mdat);
    ylim([min(ylim), 1.1*max(ylim)])
    text(min(xlim), max(ylim), descr, ...
         'color', [0,0,0] , 'VerticalAlignment', 'top', ...
         'FontName', 'Monospaced', 'FontSize', 9);

    title(sprintf('Experiment: %d training sets, %d par combinations', ...
                  size(nums,1), size(nums,2)));

end