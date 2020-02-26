% Initialize QG
nx = 128;
ny = 128;
n  = nx * ny * 2;
qg = QG(nx, ny, 1);

% Henk QG params:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365*day;

Re   = 10000;
wind = 0;

qg.set_par(11, wind);  % wind stress (stirring) amplitude
qg.set_par(5,    Re);  % Reynolds number
qg.set_par(2,   0.0);  % no rotation

rhs = @ (x) qg.rhs(x);

xgrid = ((1:nx)-1)*2*pi/nx;
ygrid = ((1:ny)-1)*2*pi/ny;

z0 = sin(4*xgrid)'*sin(4*ygrid) + ...
     0.4*cos(3*xgrid)'*cos(3*ygrid) + ...
     0.3*cos(5.0*xgrid)'*cos(5.0*ygrid) + ...
     0.02*sin(xgrid) + 0.02*cos(ygrid);


z0 = sin(5*xgrid)'*sin(5*ygrid) + ...
     0.5*cos(4*xgrid+0.1)'*cos(4*ygrid+0.1) + ...
     0.2*sin(3*xgrid+0.3)'*sin(3*ygrid+0.3);

z0 = (rand(nx,ny)-0.5)+(sin(4*xgrid)'*sin(4*ygrid));

x0 = zeros(n,1);
x0(1:2:end) = 0.2*z0(:)/(3600*24/tdim); % nondimensional and
                                        % realistic vorticity

dt = 0.01;

th = 1.0;          % theta
s  = 1.0/(dt*th);
B  = qg.mass(n);

F = @(x) qg.rhs(x) ;

x  = x0;
F0 = F(x);

kDes = 3.3;
t = 0;
states = [];
times = [];
storeTime = 0;
Tend = Re;     % e-folding timescale largest mode (?)

tic
while t < Tend
    fprintf(' t = %2.2e years,  \n',  t / year);
    fprintf('dt = %2.2e days \n Newton: \n', dt / day);
    fprintf('start Newton \n')

    for k = 1:10
        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
        J   = qg.jacobian(x, s);
        dx  = J \ rhs;
        x   = x + dx;
        if norm(dx,2) < 1e-3
            fprintf('||dx|| = %2.5e\n', norm(dx));
            break;
        end
    end   
    t  = t + dt;
    dt = kDes / k * dt;
    s  = 1.0 / (dt*th);
         
    x0 = x;
    F0 = F(x);

    if t > storeTime || t > Tend
        states = [states, x];
        times  = [times,  t];
        
        subplot(1,2,1);
        plotQG(nx,ny,2,x);
        titleString = sprintf('t = %f year', t / year);
        title(titleString);
        
        subplot(1,2,2);
        plotQG(nx,ny,1,3600*24/tdim*x,false);
        titleString = sprintf('t = %f year', t / year);
        title(titleString);
        
        fnamebase = ['N',num2str(nx), '_Re', num2str(Re), '_Tend', ...
                 num2str(Tend), '_wind', num2str(wind)];
        
        exportfig([fnamebase,'.eps'],10,[30,10]);

        fprintf('saving data to %s\n', [fnamebase,'.mat']);
        
        save([fnamebase,'.mat'], 'states', 'times', 'nx', 'ny', 'Re', ...
             't', 'dt', 'wind');
        

        storeTime = t + 0.1*year;
    end
end
toc
toc / 3600

