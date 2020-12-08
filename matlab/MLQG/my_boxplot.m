function [] = my_boxplot(varargin)

    switch nargin
      case 1
        array    = varargin{1};
        x_index  = 1:size(array,2);

      case 2
        x_index  = varargin{1};
        array    = varargin{2};

      otherwise
        error('Unexpected input');
    end

    alpha     = 0.2;
    ylimMax   = 0.0;
    plot_mean    = false;
    plot_scatter = true;
    Q = zeros(numel(x_index), 3);

    for idx = x_index
        % scatterplot of all data
        if plot_scatter
            scatter(repmat(idx,1,size(array,1)), array(:,idx), ...
                    'k+', 'markeredgealpha', alpha, ...
                    'sizedata', 10, ...
                    'markerfacecolor', 'k', ...
                    'markerfacealpha', alpha, ...
                    'linewidth', 0.01); hold on
        end

        % plot median
        q1 = quantile(array(:,idx),0.25);
        q2 = quantile(array(:,idx),0.5);
        q3 = quantile(array(:,idx),0.75);
        Q(idx, :) = [q1,q2,q3];

        plot(idx, q2, 'k.','markersize', 15, 'linewidth',1);
        hold on

        % plot quantiles
        plot(repmat(idx,1,2), ...
             [q1, q3], ...
             'k.-','markersize', 12,'linewidth',1);

        if plot_mean
            plot(repmat(idx,1,2), ...
                 mean(array(:,idx)), 'k*', 'markersize', 12)
        end

        ylimMax = max(ylimMax, 1.2*quantile(array(:,idx), 0.75));
    end
    ht = plot(x_index, Q');
    uistack(ht, 'bottom');

    xlim([min(x_index)-0.5, max(x_index)+0.5]);
    xticks([x_index]);
    ylim([0,ylimMax]);

    hold off
end