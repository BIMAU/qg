% Initialize QG
nx = 512;
ny = 512;
n  = nx * ny * 2;
qg = QG(nx, ny, 1);

% qg.set_par(11, 0.1);  % wind stress
                        % qg.set_par(5, 45);    % Reynolds number

% load testing data
% testdata = load('testdata.mat');
% xr = testdata.xr(1:n);

%% Test 1: analytical jacobian vs numerical jacobian

% J  = -qg.jacobian(xr, 0.0);
% 
% F0 = qg.rhs(xr);
% 
% pert = 1e-6;
% err  = zeros(n,1);
% D    = zeros(n,n);
% Jn   = zeros(n,n);
% 
% v = zeros(n,1);
% for i = 1:n
    %     v(i) = pert;
    %     Jn(:,i) = (qg.rhs(xr+v) - qg.rhs(xr-v)) / (2*pert);
    %     v(i) = 0;
    %     
    %     D(:,i)  = abs(Jn(:,i)  - J(:,i));
    %     err(i)  = max(D(:,i));
    % end
    %

%%  --> check dit in relatieve zin...
% Jn = sparse(Jn);
% D(abs(D) < 1) = 0.0;
% D = sparse(D);
% 
% assert(max(err) < 1e-4);

%vsm(J)
%vsm(Jn-J)

%% Test 2: perform a few backward Euler time steps

% Edeling/Verkley params:
nu  = 1.79957939712e-05;
mu  = 0.00176358780918;
L   = 6.371e6;
T   = 1/(7.292e-5);
nud = nu * L^2 / T
mud = mu / T

% Henk QG params:
Ldim = 1e6;
Udim = 1.6e-2;
tdim = Ldim / Udim; % in seconds

day = 3600 * 24 / tdim;

Ah = 1.8e-05;
Ldim*Udim/Ah
Re   = 8000;
wind = 0;

% qg.set_par(1, 1);    % alpha_tau
qg.set_par(11, wind);  % enable  wind stress
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

%z0 = 1*(rand(nx,ny)-0.5);

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
days = 600;
tic
while t < days*day
    fprintf(' t = %2.2e days,  \n',  t / day);
    fprintf('dt = %2.2e days \n Newton: \n', dt / day);
    fprintf('start Newton\n')

    for k = 1:10

        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
        J   = qg.jacobian(x, s);
        dx  = J \ rhs;
        x   = x + dx;

        fprintf('||dx|| = %2.5e\n', norm(dx));
        if norm(dx,2) < 1e-5
            break;
        end
    end   
    t  = t + dt;
    dt = kDes / k * dt;
    s  = 1.0 / (dt*th);
         
    %assert(norm(dx,2) < 1e-7);
    
    x0 = x;
    F0 = F(x);

    if t > storeTime || t > days*day
        states = [states, x];
        times  = [times,  t];

% $$$         figure(1)
% $$$         plotQG(nx,ny,2,x)
% $$$         titleString = sprintf('t = %f days', t / day);
% $$$         title(titleString);
% $$$         %        exportfig(['psi_Re',num2str(Re),'.eps'])

        figure(2)
        plotQG(nx,ny,1,3600*24/tdim*x,false)
        titleString = sprintf('t = %f days', t / day);
        title(titleString);
        %        exportfig(['zeta_Re',num2str(Re),'.eps'])

        storeTime = t + day;
        drawnow
    end
end
toc

save(['N',num2str(nx),'_Re',num2str(Re),'_days',num2str(days),'.mat']...
     ,'states','times','nx','ny','Re','t','dt','wind');