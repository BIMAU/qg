
% Henk QG params:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365*day;

fname = 'N256vortT300days.avi';

writerObj = VideoWriter(fname, 'Motion JPEG AVI');
writerObj.FrameRate = 22;
writerObj.Quality = 90;
open(writerObj);
set(0,'DefaultFigureWindowStyle','normal')
fhandle = figure('units','pixels','position',[100,100,900,800]);
set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
set(gca,'color','w','fontsize',15);

% relative vorticity (1/day):
scaling = 3600*24/tdim
crange = [-0.2,0.2];

N = size(states,2);

for i = max(1,N-300):N
    plotQG(nx,ny,1,scaling*states(:,i),false)
    caxis(crange);
    titleString = sprintf('rel. vorticity (day^{-1}), t = %3.0fd', (times(i)-times(1)) / day);
    title(titleString);
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);
    