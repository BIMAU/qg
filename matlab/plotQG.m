function [] = plotQG(n,m,XX,state,contours)
% Plotting QG model:
% n,m: horizontal grid size
% XX: 1,2 index of field you want to plot

    if nargin < 5
        contours = true;
    end

    % It's QG so nun = 2
    nun = 2;
    plotfield = reshape(state(XX:nun:end),n,m);
    
    % grid
    for i=1:n
        x(i) = (i-1)/(n-1);
    end

    for j=1:m
        y(j) = (j-1)/(m-1); 
    end

    % scaling
    maxp = max(max(abs(plotfield))); 
    % plot

    imagesc(x,y,plotfield');
    if contours
        hold on
        contour(x,y,plotfield',10,'k'); hold off
    end

    % colorbar
    set(gca,'ydir','normal')
    
    if (XX == 1)
        caxis([-maxp, maxp]);
    end

end