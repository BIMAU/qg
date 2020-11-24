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
    for idx = x_index
        scatter(repmat(idx,1,size(array,1)), array(:,idx), ...
                'k+', 'markeredgealpha', alpha, ... 
                'sizedata', 10, ...
                'markerfacecolor', 'k', ...
                'markerfacealpha', alpha, ...
                'linewidth', 0.01); hold on
        plot(idx, median(array(:,idx)), 'k.','markersize', 15, 'linewidth',1); 
        hold on
        plot(repmat(idx,1,2), ...
             [quantile(array(:,idx),0.25), quantile(array(:,idx),0.75)], ...
             'k.-','markersize', 12,'linewidth',1); 

    end
    hold off
end