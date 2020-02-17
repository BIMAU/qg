% Initialize QG
nx = 128;
ny = 128;
n  = nx * ny * 2;
qg = QG(nx, ny, 0);

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
Re = 7000;
wind = 20;
%qg.set_par(1, 1);   % alpha_tau
qg.set_par(11, wind);  % enable  wind stress
qg.set_par(5, Re);   % Reynolds number
qg.set_par(2, 0.0);  % no rotation

rhs = @ (x) qg.rhs(x);

xgrid = linspace(0,2*pi,nx);
ygrid = linspace(0,2*pi,ny);

%z0 = sin(4*xgrid)'*sin(4*ygrid) + ...
%     0.4*cos(3*xgrid)'*cos(3*ygrid) + ...
%     0.3*cos(5.0*xgrid)'*cos(5.0*ygrid) + ...
%     0.02*sin(xgrid) + 0.02*cos(ygrid);
 
z0 = 0.01*(rand(nx,ny)-0.5);

x0 = zeros(n,1);
x0(1:2:end) = z0(:);

dt = 0.001;

th = 1.0;          % theta
s  = 1.0/(dt*th);
B  = qg.mass(n);

F = @(x) qg.rhs(x) ;

x  = x0;
F0 = F(x);

kDes = 3.5;
figure(1)
t = 0;
states = [];
times = [];
storeTime = 0;
tic
while t < 365*day
    fprintf('dt = %2.2e, Newton: \n', dt);
    %[L,U] = ilu(J, struct('type','ilutp','droptol',1e-5));
    fprintf('start Newton\n')
    for k = 1:10
        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
        J   = qg.jacobian(x, s);
        dx  = J \ rhs;

        %qg.jacob(x,s);
        %qg.compute_precon();
        %dx = qg.solve(-rhs);
        x  = x + dx;
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

    if t > storeTime
        states = [states, x];
        times  = [times, t];    
        plotQG(nx,ny,2,x)
        titleString = sprintf('t = %f days', t / day);
        title(titleString);
        exportfig(['out_D_Re',num2str(Re),'_wind',num2str(wind),'.eps'])
        storeTime = t + day;
    end
end
toc

save(['D_N128Re',num2str(Re),'_wind',[num2str(wind)],'.mat'],'states','times');