fname = 'N256eddyforcing2.avi';
writerObj = VideoWriter(fname, 'Motion JPEG AVI');
writerObj.FrameRate = 22;
writerObj.Quality = 90;
open(writerObj);
set(0,'DefaultFigureWindowStyle','normal')
fhandle = figure('units','pixels','position',[100,100,1400,600]);
set(gca,'position',[0.05 0.1 .92 0.85],'units','normalized');
set(gca,'color','w','fontsize',15);

% Initialize QG
nun = 2;
nx = 256;
ny = 256;
qg = QG(nx, ny, 1);

% QG params:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365*day;

Re   = 4e4;   % Reynolds number
ampl = 0.5;   % Forcing amplitude
Tend = 100;   % End time, nondim. timescale is in years
dt   = 0.01;  % Time step 
th   = 1.0;   % Theta
t0   = 0;     % time

% Set parameters in QG
qg.set_par(18,  0.0);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qg.set_par(11, ampl);  % stirring amplitude
qg.set_par(5,    Re);  % Reynolds number

% Also create a coarse variant
ff = 2; % coarsening factor
nxc = 256 / ff;
nyc = 256 / ff;
qg_c = QG(nxc, nyc, 1);

% Set parameters in QG
qg_c.set_par(18,  0.0);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qg_c.set_par(11, ampl);  % stirring amplitude
qg_c.set_par(5, Re/10);  % Reynolds number

% load data for the fine setup 
% data = load('N256_Re4.0e+04_Tstart135_Tend235_F0.5.mat');

assert(nx == data.nx);

scaling = 3600*24/tdim;
crange = [-0.2,0.2];
Frange = [-10,10];

N = size(data.states,2);

for i = N-300:N
    xf = data.states(:,i);
    
    Flast     = qg.rhs(xf);
    Flast_av  = average(Flast,nx,ny,nun,ff); 
    Fcoarse   = qg_c.rhs(average(xf,nx,ny,nun,ff));

    subplot(1,2,1)
    plotQG(nx,ny,1,scaling*xf,false)
    caxis(crange);
    titleString = sprintf('Fine vorticity (day^{-1}), t = %3.1fy', ...
                          data.times(i) / year);
    title(titleString);
    
    subplot(1,2,2)
    plotQG(nxc, nyc, 1, scaling*(Flast_av-Fcoarse), false);
    colorbar
    caxis(Frange);
    title('Eddy forcing for coarse model (Re/10)')

    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);



