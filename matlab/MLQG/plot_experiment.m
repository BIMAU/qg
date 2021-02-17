dir = ['data/experiments/', 'NR_4000-16000_ESN1_MDL0', '/', 'parallel', '/'];
dir = ['data/experiments/', 'NR_500-16000_ESN1_MDL0', '/', 'parallel', '/'];
dir = ['data/experiments/', 'NR_500-8000_ESN1_MDL1', '/', 'parallel', '/'];
dir = ['data/experiments/', 'NR_4000-8000_RF_1-16_ESN1_MDL1', '/', 'parallel', '/'];

[errs,nums,pids,mdat] = gather_plotdata(dir, 2);

% ignore rows with nans
plotIds = find(~isnan(sum(nums,2)));
nums = nums(plotIds,:);

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
maxValues = Nvalues(xlab_index);
Ntotal    = size(nums,2);
Nboxplots = Ntotal/maxValues;

% plot results
f = [];
clrs = lines(Nboxplots);
for i = 1:Nboxplots
    range = (i:Nboxplots:Ntotal)
    f{i} = my_boxplot(nums(:,range), {clrs(i,:), clrs(i,:)}); 
    hold on
end
hold off

xticklabels(mdat.hyp_range(exp_ind{xlab_index}, range));
xtickangle(45);
xlabel(mdat.xlab{xlab_index});
ylabel(mdat.ylab);
grid on;

% for combined experiments and multiple boxplots we need a legend
if Nexp == 2
    str = cell(Nboxplots,1);
    for i = 1:Nboxplots
        value = mdat.hyp_range(exp_ind{I(2)}, i);
        str{i} = sprintf('%2d', value);
    end
    legend([f{:}], str, 'location', 'northwest')    
end    

% create description
descr = create_description(mdat)
text(0.6*max(xlim), max(ylim), descr, ...
     'color', [0,0,0] , 'VerticalAlignment', 'top', ...
     'FontName', 'Monospaced', 'FontSize', 9);