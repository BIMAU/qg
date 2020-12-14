function [] = my_boxplot(varargin)

    switch nargin
      case 1
        array    = varargin{1};
        colors   = {'k', 'k'};
        
      case 2
        array    = varargin{1};
        colors   = varargin{2};
        
      otherwise
        error('Unexpected input');
    end
    
    

    x_index   = 1:size(array,2);
    alpha     = 0.2;
    ylimMax   = 0.0;
    plot_mean    = true;
    plot_scatter = true;
    Q = zeros(numel(x_index), 3);

    for idx = x_index
        arr = array(:,idx);        
        % scatterplot of all data
        
        if plot_scatter
            scatter(repmat(idx,1,size(array,1)), arr, ...
                    '+', 'markeredgealpha', alpha, ...
                    'sizedata', 10, ...
                    'markerfacecolor', colors{1}, ...
                    'markeredgecolor', colors{1}, ...
                    'markerfacealpha', alpha, ...
                    'linewidth', 0.01); hold on
        end

        % plot median
        q1 = quantile(arr,0.25);
        q2 = quantile(arr,0.5);
        q3 = quantile(arr,0.75);
        Q(idx, :) = [q1,q2,q3];

        plot(idx, q2, '.','markersize', 15, 'linewidth',1, 'color', colors{1});
        hold on

        % plot quantiles
        plot(repmat(idx,1,2), ...
             [q1, q3], ...
             '.-','markersize', 12,'linewidth',1, 'color', colors{1});

        if plot_mean
            mn = mean(arr(~isnan(arr)));
            st = std(arr(~isnan(arr)));
            plot(idx, mn, ...
                 '*', 'markersize', 8, 'color', colors{2});
        end

        ylimMax = max(ylimMax, 1.25*quantile(arr, 0.75));
    end
    ht = plot(x_index, Q, 'color', colors{2});
    uistack(ht, 'bottom');

    xlim([min(x_index)-0.5, max(x_index)+0.5]);
    xticks([x_index]);
    ylim([0,ylimMax]);

    hold off
end