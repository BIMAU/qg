% load existing states array
%fprintf('loading existing fullmodel data ...\n'); tic;
%fmdata = load('data/fullmodel/N256_Re4.0e+04_Tstart142_Tend151_F0.5_Stir0_Rot1.mat');
%fprintf('loading existing fullmodel data ... done (%fs)\n', toc)

load 'data/fullmodel/N128_Re1.0e+04_Tstart151_Tend179_F0.5_Stir0_Rot1.mat'


% [~,~,R1,~] = gridTransfers(fmdata.nx, 'periodic');
% R1         = kron(R1, [1,0;0,1]);
% [~,~,R2,~] = gridTransfers(fmdata.nx/2, 'periodic');
% R2         = kron(R2, [1,0;0,1]);
% x_init     = R1*fmdata.states;
% x_init     = x_init(:,end);
x_init = states(:,end);

%%

% Initialize QG
nx = 128;
ny = 128;
n  = nx * ny * 2;
qg = QG(nx, ny, 1);

% QG params:
Ldim = 1e6;
Udim = 3.171e-2;
tdim = Ldim / Udim; % in seconds
day  = 3600 * 24 / tdim;
year = 365*day;

% Job parameters
restartFlag      = true;  % restart from existing states array
adaptiveTimeStep = false;
storeTimeIncr = 100*day;   % save to mat and eps
rotation = true;          % enable or disable rotation (beta)

Re   = 1e4;       % Reynolds number
ampl = 0.5;       % Forcing amplitude
stir = 0;         % stirring type: 0 = cos(5x), 1 = sin(16x)
Tend = 10000*day; % End time, nondim. timescale is in years
dt   = 1.0*day;   % Time step
th   = 1.0;       % Theta
t0   = 0;         % time

% Set parameters in QG
qg.set_par(18, stir);  % stirring type: 0 = cos(5x), 1 = sin(16x)
qg.set_par(11, ampl);  % stirring amplitude
qg.set_par(5,    Re);  % Reynolds number
if ~rotation
    qg.set_par(2,   0.0);  % no rotation
end

% Several initial vorticity profile choices can be made here:
xgrid = ((1:nx)-1)*2*pi/nx;
ygrid = ((1:ny)-1)*2*pi/ny;
rng(7);

% z0 = sin(4*xgrid)'*sin(4*ygrid) + ...
%      0.4*cos(3*xgrid)'*cos(3*ygrid) + ...
%      0.3*cos(5.0*xgrid)'*cos(5.0*ygrid) + ...
%      0.02*sin(xgrid) + 0.02*cos(ygrid);
%
%
% z0 = sin(5*xgrid)'*sin(5*ygrid) + ...
%      0.5*cos(4*xgrid+0.1)'*cos(4*ygrid+0.1) + ...
%      0.2*sin(3*xgrid+0.3)'*sin(3*ygrid+0.3);
%
% z0 = 1*(rand(nx,ny)-0.5)+(sin(4*xgrid)'*sin(4*ygrid));

z0 = 0;
for i = 1:4
    z0 = z0 + 1/i*(sin((4+0.1*randn())*xgrid+0.1*randn())'*sin((4+0.1* ...
                                                      randn())*ygrid+0.1*randn()));
end

x0 = zeros(n,1);
x0(1:2:end) = 0.2*z0(:) / (3600*24/tdim); % nondimensional and
                                          % realistic vorticity

% we may also restart from existing states array
if restartFlag
    if exist('x_init') && exist('fmdata')
        x0 = x_init;
        assert(numel(x0) == n);
        t0 = fmdata.times(end);
        Tend = fmdata.times(end) + Tend;
    else
        fprintf('Cannot restart, no states array available\n');
        return;
    end
end

s  = 1.0/(dt*th);   % Contribution on diagonal of Jacobian
B  = qg.mass(n);    % Mass matrix

F = @(x) qg.rhs(x); % Right hand side

x  = x0;
F0 = F(x);

kDes   = 3.3;    % optimal number of Newton iterations
states = [];
times  = [];
storeTime = 0;
t = t0;
%%
while t < Tend
    fprintf(' t = %2.2e years,  \n',  t / year);
    fprintf('dt = %2.2e days \n Newton: \n', dt / day);
    fprintf('start Newton \n');

    for k = 1:10
        rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
        J   = qg.jacobian(x, s);
        dx  = J \ rhs;
        x   = x + dx;
        if norm(dx,2) < 1e-3
            fprintf('||dx|| = %2.5e, k = %d\n', norm(dx),k);
            rhs = B*(x-x0)/(dt*th) + F(x) + (1-th)/th * F0;
            fprintf('||F|| = %2.5e \n', norm(rhs, 2));
            break;
        end
    end
    t  = t + dt;
    if adaptiveTimeStep
        dt = kDes / k * dt;
        s  = 1.0 / (dt*th);
    end

    x0 = x;
    F0 = F(x);

    states = [states, x];
    times  = [times,  t];

    if t > storeTime || t > Tend

        subplot(2,2,1);
        plotQG(nx,ny,2,x);
        titleString = sprintf('psi, t = %4.0fd', (t-times(1)) / day);
        title(titleString);

        subplot(2,2,2);
        plotQG(nx,ny,1,3600*24/tdim*x,false);
        titleString = sprintf('vorticity');
        title(titleString);

        [u,v] = qg.compute_uv(x);
        subplot(2,2,3);
        u = reshape(u,nx,ny);
        v = reshape(v,nx,ny);
        imagesc((u.^2+v.^2)')
        titleString = sprintf('u^2 + v^2');
        title(titleString);

        subplot(2,2,4)
        plotQGspectrum(qg, nx, ny, x, 5);

        ReStr = sprintf('_Re%1.1e',Re);
        fnamebase = [ 'N', num2str(nx), ReStr, '_Tstart', ...
                      num2str(round(t0)), '_Tend', ...
                      num2str(round(Tend)), '_F', ...
                      num2str(ampl), '_Stir', num2str(stir), ...
                      '_Rot',num2str(rotation)];


        exportfig(['data/eps/',fnamebase,'.eps'],10,[40,10]);

        fprintf('saving data to %s\n', [fnamebase,'.mat']);

        % this is slow
        save(['data/fullmodel/',fnamebase,'.mat'], ...
             'states', 'times', 'nx', 'ny', 'Re', ...
             't', 'dt', 'ampl', 'stir', '-v7.3');

        storeTime = t + storeTimeIncr;
    end
end