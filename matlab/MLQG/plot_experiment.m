dir = ['data/experiments/', 'NR_4000-16000_ESN1_MDL0', '/', 'parallel', '/'];
dir = ['data/experiments/', 'NR_500-16000_ESN1_MDL0', '/', 'parallel', '/'];

[errs,nums,pids,mdat] = gather_plotdata(dir, 2);

% ignore nans
plotIds = find(~isnan(sum(nums,2)));
plotJds = find(~isnan(sum(nums,1)));
nums = nums(plotIds,:);

% plot results
f = my_boxplot(nums);

% number of experiments
exp_ind = mdat.exp_ind;
Nexp    = numel(exp_ind);

labels  = [];
Nvalues = [];

for i = 1:Nexp
    labels{i} = mdat.hyp_range(exp_ind{i}, :);
    Nvalues(i) = numel(unique(labels{i}));
end

[~, I] = sort(Nvalues, 'descend')
xlab_index = I(1);

xticklabels(mdat.hyp_range(exp_ind{xlab_index},:));
xtickangle(45);
xlabel(mdat.xlab{xlab_index});
ylabel(mdat.ylab);
grid on;

% create description
descr = create_description(mdat)
text(0.7, 0.8*max(ylim), descr, 'color', [0,0,0] , 'FontSize', 10);