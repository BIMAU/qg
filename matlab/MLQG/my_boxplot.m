function [] = my_boxplot(varargin)
    
    switch nargin
      case 1
        array    = varargin{1};
        x_index  = 1:size(array,2);

      case 2
        x_index  = varargin{1}; % assuming input is a string
        array    = varargin{2};

      otherwise
        error('Unexpected input');
    end
    
    alpha = 0.1;
    ylimMax = 0.0;
    for idx = x_index
        scatter(repmat(idx,1,size(array,1)), array(:,idx), ...
                'k+', 'markeredgealpha', alpha, ... 
                'sizedata', 10, ...
                'markerfacecolor', 'k', ...
                'markerfacealpha', alpha, ...
                'linewidth', 0.01); hold on
        plot(idx, quantile(array(:,idx),0.5), 'k.','markersize', 15, 'linewidth',1); 
        hold on
        plot(repmat(idx,1,2), ...
             [quantile(array(:,idx),0.25), quantile(array(:,idx),0.75)], ...
             'k.-','markersize', 12,'linewidth',1); 
        
        ylimMax = max(ylimMax, 2*quantile(array(:,idx),0.75));
    end    
    xlim([min(x_index)-0.5, max(x_index)+0.5]);
    ylim([0,ylimMax]);
    hold off    
end